//! Structures for representing sequence to graph alignments.
//!
//! An [`Alignment`] object represents the alignment of a query sequence to a target path in a graph.
//! It corresponds to a single line in a GAF file or a row in table `Alignments` in a [`crate::GAFBase`].
//!
//! The GAF format is a text-based format for representing sequence alignments to a graph.
//! See [the specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) for an overview.
//! Some details are better documented in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html#10).

use crate::utils;

use std::fmt::Display;
use std::io::{Read, Write};
use std::ops::Range;
use std::str;

use zstd::stream::Encoder as ZstdEncoder;
use zstd::stream::Decoder as ZstdDecoder;

use gbwt::{GBWT, Orientation, Pos};
use gbwt::support::{self, ByteCode, ByteCodeIter, RLE, Run, RLEIter};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

// FIXME examples
/// An alignment between a query sequence and a target path in a graph.
///
/// This object corresponds either to a line in a GAF file or to a row in table `Alignments` in [`crate::GAFBase`].
/// When the alignment is built from a GAF line, the target path is stored explicitly.
/// For alignments stored in a database, only the GBWT starting position is stored.
/// See [`TargetPath`] for details.
#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    /// Name of the query sequence.
    pub name: String,
    /// Length of the query sequence.
    pub seq_len: usize,
    /// Aligned interval of the query sequence.
    pub seq_interval: Range<usize>,
    /// Target path in the orientation of the query sequence.
    pub path: TargetPath,
    /// Length of the target path.
    pub path_len: usize,
    /// Aligned interval of the target path.
    pub path_interval: Range<usize>,
    /// Number of matches in the alignment.
    pub matches: usize,
    /// Number of mismatches and gaps in the alignment.
    pub edits: usize,
    /// Mapping quality.
    pub mapq: Option<usize>,
    /// Alignment score.
    pub score: Option<isize>,
    /// Base quality values for the query sequence.
    pub base_quality: Option<Vec<u8>>,
    /// Difference string, or an empty vector if not present.
    pub difference: Vec<Difference>,
    /// Information about the paired alignment.
    pub pair: Option<PairedRead>,
    /// Optional typed fields that have not been interpreted.
    pub optional: Vec<TypedField>,
}

/// Construction from a GAF line and matching the alignment with a GBWT path.
impl Alignment {
    // Number of mandatory fields in a GAF line.
    const MANDATORY_FIELDS: usize = 12;

    // Placeholder value for a missing mapping quality.
    const MISSING_MAPQ: usize = 255;

    // The field is empty and the value is missing; typically used with unaligned sequences.
    const MISSING_VALUE: [u8; 1] = [b'*'];

    // Parses a string field from a GAF field.
    fn parse_string(field: &[u8], field_name: &str) -> Result<String, String> {
        String::from_utf8(field.to_vec()).map_err(|err| {
            format!("Invalid {}: {}", field_name, err)
        })
    }

    // Parses an unsigned integer from a GAF field.
    // Returns `0` if the value is missing.
    fn parse_usize(field: &[u8], field_name: &str) -> Result<usize, String> {
        if field == Self::MISSING_VALUE {
            return Ok(0);
        }
        let number = str::from_utf8(field).map_err(|err| {
            format!("Invalid {}: {}", field_name, err)
        })?;
        number.parse().map_err(|err| {
            format!("Invalid {}: {}", field_name, err)
        })
    }

    // Parses an interval from two GAF fields.
    fn parse_interval(start: &[u8], end: &[u8]) -> Result<Range<usize>, String> {
        let start = Self::parse_usize(start, "interval start")?;
        let end = Self::parse_usize(end, "interval end")?;
        Ok(start..end)
    }

    // Parses an orientation from a GAF field.
    // Returns [`Orientation::Forward`] if the value is missing.
    fn parse_orientation(field: &[u8], field_name: &str) -> Result<Orientation, String> {
        if field == Self::MISSING_VALUE {
            return Ok(Orientation::Forward);
        }
        if field.len() != 1 {
            return Err(format!("Invalid {}: {}", field_name, String::from_utf8_lossy(field)));
        }
        match field[0] {
            b'+' => Ok(Orientation::Forward),
            b'-' => Ok(Orientation::Reverse),
            _ => Err(format!("Invalid {}: {}", field_name, String::from_utf8_lossy(field))),
        }
    }

    // Parses an oriented path from a GAF field.
    // Returns an empty path if the value is missing.
    fn parse_path(field: &[u8]) -> Result<Vec<usize>, String> {
        let mut result = Vec::new();
        if field == Self::MISSING_VALUE {
            return Ok(result);
        }

        let mut start = 0;
        while start < field.len() {
            let orientation = match field[start] {
                b'>' => Orientation::Forward,
                b'<' => Orientation::Reverse,
                _ => return Err(format!("Invalid segment orientation: {}", String::from_utf8_lossy(field))),
            };
            start += 1;
            let end = field[start..].iter().position(|&c| c == b'>' || c == b'<').map_or(field.len(), |x| start + x);
            let node = str::from_utf8(&field[start..end]).map_err(|err| {
                format!("Invalid segment name: {}", err)
            })?.parse().map_err(|_| {
                String::from("Only numerical segment names are supported")
            })?;
            result.push(support::encode_node(node, orientation));
            start = end;
        }
        Ok(result)
    }

    // Parses a pair name from the value of a typed field.
    fn parse_pair(value: Vec<u8>, is_next: bool, output: &mut Option<PairedRead>) -> Result<(), String> {
        if output.is_some() {
            return Err(String::from("Multiple pair fields"));
        }
        let name = String::from_utf8(value).map_err(|err| {
            format!("Invalid pair name: {}", err)
        })?;
        *output = Some(PairedRead {
            name,
            is_next,
            is_proper: false,
        });
        Ok(())
    }

    /// Parses an alignment from a GAF line.
    ///
    /// Returns an error if the line cannot be parsed.
    /// The line may end with up to one endline character, which is ignored.
    /// The returned alignment stores the target path explicitly.
    /// Parsing is based on bytes rather than characters to avoid unnecessary UTF-8 validation.
    ///
    /// If a difference string is present, some numerical fields will be recalculated from it.
    /// These include interval ends on both the query and the target, as well as the number of matches and edits.
    /// This behavior is justified, because some aligners may not calculate these values correctly.
    pub fn from_gaf(line: &[u8]) -> Result<Self, String> {
        // Check for an endline character which may be present.
        let line = if line.last() == Some(&b'\n') {
            &line[..line.len() - 1]
        } else {
            line
        };

        // Split the line into fields.
        let fields = line.split(|&c| c == b'\t').collect::<Vec<_>>();
        if fields.len() < Self::MANDATORY_FIELDS {
            let line = String::from_utf8_lossy(line);
            let message = format!("GAF line with fewer than {} fields: {}", Self::MANDATORY_FIELDS, line);
            return Err(message);
        }

        // Query sequence.
        let name = Self::parse_string(fields[0], "query sequence name")?;
        let seq_len = Self::parse_usize(fields[1], "query sequence length")?;
        let mut seq_interval = Self::parse_interval(fields[2], fields[3])?;

        // Target path.
        let orientation = Self::parse_orientation(fields[4], "target orientation")?;
        let mut path = Self::parse_path(fields[5]).map_err(|err| {
            format!("Invalid target path: {}", err)
        })?;
        if orientation == Orientation::Reverse {
            support::reverse_path_in_place(&mut path);
        }
        let path = TargetPath::Path(path);
        let path_len = Self::parse_usize(fields[6], "target path length")?;
        // TODO: This is a hack. VG sometimes gets the target interval end wrong.
        // If we can't parse the interval, we will later try to parse the start and infer the length from the difference string.
        let path_interval = Self::parse_interval(fields[7], fields[8]);

        // Alignment statistics.
        let mut matches = Self::parse_usize(fields[9], "matches")?;
        let alignment_len = Self::parse_usize(fields[10], "alignment length")?;
        let mut edits = if matches <= alignment_len { alignment_len - matches } else { 0 };
        let mapq = Self::parse_usize(fields[11], "mapping quality")?;
        let mapq = if mapq == Self::MISSING_MAPQ { None } else { Some(mapq) }; // TODO: Too large values?

        // Optional fields.
        let mut score = None;
        let mut base_quality = None;
        let mut difference = Vec::new();
        let mut pair = None;
        let mut properly_paired = None;
        let mut optional = Vec::new();
        for field in fields[Self::MANDATORY_FIELDS..].iter() {
            let parsed = TypedField::parse(field)?;
            match parsed {
                TypedField::Int([b'A', b'S'], value) => {
                    if score.is_some() {
                        return Err(String::from("Multiple alignment score fields"));
                    }
                    score = Some(value);
                },
                TypedField::String([b'b', b'q'], value) => {
                    if base_quality.is_some() {
                        return Err(String::from("Multiple base quality fields"));
                    }
                    base_quality = Some(value.to_vec());
                },
                TypedField::String([b'c', b's'], value) => {
                    if !difference.is_empty() {
                        return Err(String::from("Multiple difference fields"));
                    }
                    difference = Difference::parse_normalized(&value)?;
                },
                TypedField::Bool([b'p', b'd'], value) => {
                    if properly_paired.is_some() {
                        return Err(String::from("Multiple properly paired fields"));
                    }
                    properly_paired = Some(value);
                },
                TypedField::String([b'f', b'n'], value) => {
                    Self::parse_pair(value, true, &mut pair)?;
                },
                TypedField::String([b'f', b'p'], value) => {
                    Self::parse_pair(value, false, &mut pair)?;
                },
                _ => { optional.push(parsed); },
            }
        }
        if pair.is_some() && properly_paired.is_some() {
            pair.as_mut().unwrap().is_proper = properly_paired.unwrap();
        }

        // TODO: Part of the interval end hack.
        let mut path_interval = match path_interval {
            Ok(interval) => interval,
            Err(_) => {
                if !difference.is_empty() {
                    let target_start = Self::parse_usize(fields[7], "target interval start")?;
                    target_start..target_start
                } else {
                    return Err(String::from("Target interval end cannot be parsed or inferred from the difference string"));
                }
            },
        };

        // If we have a difference string, recalculate the redundant numerical fields.
        if !difference.is_empty() {
            let (query_len, target_len, num_matches, num_edits) = Difference::stats(&difference);
            seq_interval.end = seq_interval.start + query_len;
            path_interval.end = path_interval.start + target_len;
            matches = num_matches;
            edits = num_edits;
        }

        // Now we have the final path interval. Flip its orientation if necessary.
        if orientation == Orientation::Reverse {
            let start = if path_interval.end < path_len { path_len - path_interval.end } else { 0 };
            let end = if path_interval.start < path_len { path_len - path_interval.start } else { 0 };
            path_interval = start..end;
        }

        Ok(Alignment {
            name, seq_len, seq_interval,
            path, path_len, path_interval,
            matches, edits, mapq, score,
            base_quality, difference, pair,
            optional,
        })
    }

    // TODO: back to a GAF line
}

//-----------------------------------------------------------------------------

/// Encoding / decoding the alignment in the format used in GAF-base.
impl Alignment {
    // Normalizes the interval and encodes it as (left flank, length, right flank).
    fn encode_coordinates(interval: Range<usize>, len: usize, include_redundant: bool, encoder: &mut ByteCode) {
        let start = if interval.start <= interval.end { interval.start } else { interval.end };
        let end = if interval.end <= len { interval.end } else { len };
        encoder.write(start);
        if include_redundant {
            encoder.write(end - start);
        }
        encoder.write(len - end);
    }

    // Encodes a signed integer. Small absolute values are represented as small numbers.
    fn encode_signed(value: isize, encoder: &mut ByteCode) {
        let value = if value < 0 { (-2 * value - 1) as usize } else { 2 * value as usize };
        encoder.write(value);
    }

    // Encodes numerical fields as a blob.
    //
    // This includes the following fields:
    //
    // * Query sequence: `seq_len`, `seq_interval.start`, (`seq_interval.end`).
    // * Target path: `path_len`, `path_interval.start`, (`path_interval.end`).
    // * Alignment statistics: (`matches`), (`edits`), `mapq`, `score`.
    //
    // The coordinates are normalized so that `start <= end <= len`.
    // If a difference string is present, the numbers in parentheses are not stored, as they can be derived from it.
    fn encode_numbers_into(&self, encoder: &mut ByteCode) {
        let include_redundant = self.difference.is_empty();

        // Coordinates.
        Self::encode_coordinates(self.seq_interval.clone(), self.seq_len, include_redundant, encoder);
        Self::encode_coordinates(self.path_interval.clone(), self.path_len, include_redundant, encoder);

        // Alignment statistics.
        if include_redundant {
            encoder.write(self.matches);
            encoder.write(self.edits);
        }
        if let Some(mapq) = self.mapq {
            encoder.write(mapq);
        }
        if let Some(score) = self.score {
            Self::encode_signed(score, encoder);
        };
    }

    // Encodes the difference string into the given encoder.
    //
    // [`Difference::End`] values are not included in the encoding, but one is always added at the end.
    fn encode_difference_into(&self, encoder: &mut RLE) {
        if !self.difference.is_empty() {
            for diff in self.difference.iter() {
                match diff {
                    Difference::Match(len) => encoder.write(Run::new(0, *len)),
                    Difference::Mismatch(base) => encoder.write(Run::new(1, utils::encode_base(*base))),
                    Difference::Insertion(seq) => {
                        let len = utils::encoded_length(seq.len());
                        encoder.write(Run::new(2, len));
                        let encoded = utils::encode_sequence(seq);
                        for byte in encoded {
                            encoder.write_byte(byte);
                        }
                    },
                    Difference::Deletion(len) => encoder.write(Run::new(3, *len)),
                    Difference::End => {},
                }
            }
        }
        encoder.write(Run::new(4, 1));
    }

    // TODO: encode optional fields

    // Decodes an interval and a length from the iterator.
    fn decode_coordinates(iter: &mut ByteCodeIter, include_redundant: bool) -> Option<(Range<usize>, usize)> {
        let left_flank = iter.next()?;
        let interval_len = if include_redundant { iter.next()? } else { 0 };
        let right_flank = iter.next()?;
        Some((left_flank..left_flank + interval_len, left_flank + interval_len + right_flank))
    }

    // Decodes a signed integer from the iterator.
    fn decode_signed(iter: &mut ByteCodeIter) -> Option<isize> {
        let value = iter.next()?;
        if value % 2 == 0 {
            Some((value / 2) as isize)
        } else {
            Some(-(value as isize + 1) / 2)
        }
    }

    // Decodes a difference string from the iterator until an End value.
    // Also needs the underlying slice for decoding insertions.
    // Returns an error if the if there is no End value.
    fn decode_difference_from(encoded: &[u8], decoder: &mut RLEIter) -> Result<Vec<Difference>, String> {
        let mut result: Vec<Difference> = Vec::new();
        while let Some(run) = decoder.next() {
            match run.value {
                0 => result.push(Difference::Match(run.len)),
                1 => result.push(Difference::Mismatch(utils::decode_base(run.len))),
                2 => {
                    let offset = decoder.offset();
                    for _ in 0..run.len {
                        let _ = decoder.byte().ok_or(String::from("Missing insertion base"))?;
                    }
                    let encoded = &encoded[offset..offset + run.len];
                    let seq = utils::decode_sequence(encoded);
                    result.push(Difference::Insertion(seq));
                },
                3 => result.push(Difference::Deletion(run.len)),
                4 => {
                    return Ok(result);
                },
                _ => return Err(format!("Invalid difference string operation: {}", run.value)),
            }
        }
        Err(String::from("Encoded difference string ended without an End value"))
    }
}

//-----------------------------------------------------------------------------

// FIXME test these
// Operations on the Alignment object.

impl Alignment {
    /// Returns `true` if the read is unaligned.
    pub fn is_unaligned(&self) -> bool {
        self.seq_interval.is_empty()
    }

    /// Returns the minimum GBWT node identifier in the target path, or [`None`] if there is no path.
    pub fn min_node(&self) -> Option<usize> {
        match self.path {
            TargetPath::Path(ref path) => path.iter().copied().min(),
            TargetPath::StartPosition(_) => None,
        }
    }

    /// Returns the maximum GBWT node identifier in the target path, or [`None`] if there is no path.
    pub fn max_node(&self) -> Option<usize> {
        match self.path {
            TargetPath::Path(ref path) => path.iter().copied().max(),
            TargetPath::StartPosition(_) => None,
        }
    }

    // FIXME: This should also work with GAF-base.
    // FIXME: error handling
    /// Sets the target path from the GBWT index if it is not already present.
    pub fn set_target_path(&mut self, index: &GBWT) {
        let mut pos = match self.path {
            TargetPath::Path(_) => return,
            TargetPath::StartPosition(pos) => Some(pos),
        };
        let mut path = Vec::new();
        while let Some(p) = pos {
            path.push(p.node);
            pos = index.forward(p);
        }
        self.path = TargetPath::Path(path);
    }
}

//-----------------------------------------------------------------------------

/// Target path in [`Alignment`].
///
/// The path is stored either explicitly as a sequence of oriented nodes or as its starting position in the GBWT index.
/// This corresponds to table `Nodes` in [`crate::GAFBase`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TargetPath {
    /// Path as a sequence of oriented nodes (GBWT node identifiers; GAF-base node handles).
    Path(Vec<usize>),
    /// Starting position in the GBWT.
    StartPosition(Pos),
}

/// Information about the paired alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairedRead {
    /// Name of the pair.
    pub name: String,
    /// Is the pair the next the next fragment?
    pub is_next: bool,
    /// Are the alignments properly paired?
    pub is_proper: bool,
}

//-----------------------------------------------------------------------------

/// An operation in a difference string describing an alignment between a query sequence and a target sequence.
///
/// This implementation supports the following operations:
///
/// * `=`: A match given as the matching sequence.
/// * `:`: A match given as the match length.
/// * `*`: A mismatch given as the target base and the query base.
/// * `+`: An insertion given as the inserted sequence.
/// * `-`: A deletion given as the deleted sequence.
///
/// The operations do not store target bases, as the query sequence can be reconstructed without that information.
/// Operation `~` (intron length and splice signal) is not supported yet.
/// Parsing is based on bytes rather than characters to avoid unnecessary UTF-8 validation.
///
/// # Examples
///
/// ```
/// use gbz_base::alignment::Difference;
///
/// let with_gaps = b":48-CAT:44+GATTACA:51";
/// let ops = Difference::parse(with_gaps);
/// assert!(ops.is_ok());
/// let ops = ops.unwrap();
/// assert_eq!(ops.len(), 5);
/// assert_eq!(ops[0], Difference::Match(48));
/// assert_eq!(ops[1], Difference::Deletion(3));
/// assert_eq!(ops[2], Difference::Match(44));
/// assert_eq!(ops[3], Difference::Insertion(b"GATTACA".to_vec()));
/// assert_eq!(ops[4], Difference::Match(51));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Difference {
    /// A match of the given length.
    Match(usize),
    /// Mismatch represented as the query base.
    Mismatch(u8),
    /// Insertion to the reference represented as the inserted sequence.
    Insertion(Vec<u8>),
    /// Deletion from the reference represented as deletion length.
    Deletion(usize),
    /// Marker for the end of the sequence.
    End,
}

impl Difference {
    /// Number of supported operation types.
    pub const NUM_TYPES: usize = 5;

    // TODO: This does not support `~` (intron length and splice signal) yet.
    const OPS: &'static [u8] = b"=:*+-";

    fn base_to_upper(c: u8) -> u8 {
        if c.is_ascii_lowercase() {
            c - b'a' + b'A'
        } else {
            c
        }
    }

    fn seq_to_upper(seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|&c| Self::base_to_upper(c)).collect()
    }

    fn matching_sequence(value: &[u8]) -> Option<Self> {
        Some(Self::Match(value.len()))
    }

    fn match_length(value: &[u8]) -> Option<Self> {
        let len = str::from_utf8(value).ok()?;
        let len = len.parse::<usize>().ok()?;
        Some(Self::Match(len))
    }

    fn mismatch(value: &[u8]) -> Option<Self> {
        if value.len() != 2 {
            return None;
        }
        Some(Self::Mismatch(Self::base_to_upper(value[1])))
    }

    fn insertion(value: &[u8]) -> Option<Self> {
        Some(Self::Insertion(Self::seq_to_upper(value)))
    }

    fn deletion(value: &[u8]) -> Option<Self> {
        Some(Self::Deletion(value.len()))
    }

    /// Parses a difference string and returns it as a vector of operations.
    ///
    /// Returns an error if the difference string is invalid.
    pub fn parse(difference_string: &[u8]) -> Result<Vec<Self>, String> {
        let mut result: Vec<Self> = Vec::new();
        if difference_string.is_empty() {
            return Ok(result);
        }
        if !Self::OPS.contains(&difference_string[0]) {
            return Err(format!("Invalid difference string operation: {}", difference_string[0] as char));
        }

        let mut start = 0;
        while start < difference_string.len() {
            let mut end = start + 1;
            while end < difference_string.len() && !Self::OPS.contains(&difference_string[end]) {
                end += 1;
            }
            let value = &difference_string[start + 1..end];
            let op = match difference_string[start] {
                b'=' => Self::matching_sequence(value),
                b':' => Self::match_length(value),
                b'*' => Self::mismatch(value),
                b'+' => Self::insertion(value),
                b'-' => Self::deletion(value),
                _ => return Err(format!("Invalid difference string operation: {}", difference_string[start] as char)),
            }.ok_or(format!("Invalid difference string field: {}", String::from_utf8_lossy(&difference_string[start..end])))?;
            result.push(op);
            start = end;
        }

        Ok(result)
    }

    /// Parses a difference string and returns it as a normalized vector of operations.
    ///
    /// The operations are merged and empty operations are removed.
    /// Returns an error if the difference string is invalid.
    pub fn parse_normalized(difference_string: &[u8]) -> Result<Vec<Self>, String> {
        let ops = Self::parse(difference_string)?;
        Ok(Self::normalize(ops))
    }

    /// Calculates various statistics from a sequence of operations.
    ///
    /// The return value is (query length, target length, matches, edits).
    pub fn stats(difference_string: &[Self]) -> (usize, usize, usize, usize) {
        let mut query_len = 0;
        let mut target_len = 0;
        let mut matches = 0;
        let mut edits = 0;
        for diff in difference_string.iter() {
            match diff {
                Self::Match(len) => {
                    query_len += len;
                    target_len += len;
                    matches += len;
                },
                Self::Mismatch(_) => {
                    query_len += 1;
                    target_len += 1;
                    edits += 1;
                }
                Self::Insertion(seq) => {
                    query_len += seq.len();
                    edits += seq.len();
                }
                Self::Deletion(len) => {
                    target_len += len;
                    edits += len;
                },
                Self::End => {},
            }
        }
        (query_len, target_len, matches, edits)
    }

    /// Returns the length of the operation.
    pub fn len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(len) => *len,
            Self::End => 0,
        }
    }

    /// Returns `true` if the operation is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the length of the operation in the target sequence.
    pub fn target_len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(_) => 0,
            Self::Deletion(len) => *len,
            Self::End => 0,
        }
    }

    /// Returns the length of the operation in the query sequence.
    pub fn query_len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(_) => 0,
            Self::End => 0,
        }
    }

    /// Merges the given operation into this operation if they can be merged.
    ///
    /// Returns `true` if the operations were merged.
    pub fn try_merge(&mut self, op: &Self) -> bool {
        match (self, op) {
            (Self::Match(len1), Self::Match(len2)) => {
                *len1 += len2;
                true
            },
            (Self::Insertion(seq1), Self::Insertion(seq2)) => {
                seq1.extend_from_slice(seq2);
                true
            },
            (Self::Deletion(len1), Self::Deletion(len2)) => {
                *len1 += len2;
                true
            },
            _ => false,
        }
    }

    /// Normalizes the sequence of operations.
    ///
    /// This merges adjacent matches and insertions and removes empty operations.
    pub fn normalize(ops: Vec<Self>) -> Vec<Self> {
        let mut result = ops;
        let mut tail = 0;
        for i in 0..result.len() {
            if result[i].is_empty() {
                continue;
            }
            if tail > 0 {
                // We need to be careful around the borrow checker.
                let (left, right) = result.split_at_mut(i);
                if left[tail - 1].try_merge(&right[0]) {
                    continue;
                }
            }
            result.swap(tail, i);
            tail += 1;
        }

        result.truncate(tail);
        result
    }
}

//-----------------------------------------------------------------------------

/// A typed optional field used in formats such as SAM, GFA, and GAF.
///
/// The field corresponds to a TAG:TYPE:VALUE string.
/// Supported types include A (single character), Z (string), i (integer), f (float), and b (boolean).
/// The field is stored as a tag and a value.
/// Parsing is based on bytes rather than characters to avoid unnecessary UTF-8 validation.
///
/// # Examples
///
/// ```
/// use gbz_base::alignment::TypedField;
///
/// let alignment_score = "AS:i:160";
/// let field = TypedField::parse(alignment_score.as_bytes());
/// assert_eq!(field, Ok(TypedField::Int([b'A', b'S'], 160)));
/// assert_eq!(field.unwrap().to_string(), alignment_score);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub enum TypedField {
    /// A single character.
    Char([u8; 2], u8),
    /// A string.
    String([u8; 2], Vec<u8>),
    /// An integer.
    Int([u8; 2], isize),
    /// A float.
    Float([u8; 2], f64),
    /// A boolean value.
    Bool([u8; 2], bool),
}

impl TypedField {
    /// Parses the field from a TAG:TYPE:VALUE string.
    ///
    /// Returns [`None`] if the field cannot be parsed or the type is unsupported.
    pub fn parse(field: &[u8]) -> Result<Self, String> {
        if field.len() < 5 || field[2] != b':' || field[4] != b':' {
            return Err(format!("Invalid typed field: {}", String::from_utf8_lossy(field)));
        }
        let tag = [field[0], field[1]];
        match field[3] {
            b'A' => {
                if field.len() != 6 {
                    return Err(format!("Invalid char field {}", String::from_utf8_lossy(field)));
                }
                Ok(TypedField::Char(tag, field[5]))
            },
            b'Z' => Ok(TypedField::String(tag, field[5..].to_vec())),
            b'i' => {
                let value = String::from_utf8_lossy(&field[5..]);
                let value = value.parse::<isize>().map_err(|err| {
                    format!("Invalid int field {}: {}", value, err)
                })?;
                Ok(TypedField::Int(tag, value))
            },
            b'f' => {
                let value = String::from_utf8_lossy(&field[5..]);
                let value = value.parse::<f64>().map_err(|err| {
                    format!("Invalid float field {}: {}", value, err)
                })?;
                Ok(TypedField::Float(tag, value))
            },
            b'b' => {
                if field.len() != 6 {
                    return Err(format!("Invalid bool field {}", String::from_utf8_lossy(field)));
                }
                match field[5] {
                    b'0' => Ok(TypedField::Bool(tag, false)),
                    b'1' => Ok(TypedField::Bool(tag, true)),
                    _ => Err(format!("Invalid bool field {}", String::from_utf8_lossy(field))),
                }
            },
            _ => Err(format!("Unsupported field type: {}", field[3] as char)),
        }
    }

    /// Returns the tag of the field.
    pub fn tag(&self) -> [u8; 2] {
        match self {
            TypedField::Char(tag, _) => *tag,
            TypedField::String(tag, _) => *tag,
            TypedField::Int(tag, _) => *tag,
            TypedField::Float(tag, _) => *tag,
            TypedField::Bool(tag, _) => *tag,
        }
    }
}

impl Display for TypedField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TypedField::Char(tag, value) => {
                write!(f, "{}{}:A:{}", tag[0] as char, tag[1] as char, *value as char)
            },
            TypedField::String(tag, value) => {
                let value = String::from_utf8_lossy(value);
                write!(f, "{}{}:Z:{}", tag[0] as char, tag[1] as char, value)
            },
            TypedField::Int(tag, value) => {
                write!(f, "{}{}:i:{}", tag[0] as char, tag[1] as char, value)
            },
            TypedField::Float(tag, value) => {
                // TODO: Precision?
                write!(f, "{}{}:f:{}", tag[0] as char, tag[1] as char, value)
            },
            TypedField::Bool(tag, value) => {
                write!(f, "{}{}:b:{}", tag[0] as char, tag[1] as char, if *value { '1' } else { '0' })
            },
        }
    }
}

//-----------------------------------------------------------------------------

/// An ad hoc bitvector for flags in [`AlignmentBlock`].
#[derive(Debug, Clone)]
pub struct Flags {
    bits: Vec<u8>,
}

impl Flags {
    /// Number of flags per alignment.
    pub const NUM_FLAGS: usize = 4;

    // Flags expressed as a bit offset.
    pub const FLAG_PAIR_IS_NEXT: usize = 0;
    pub const FLAG_PAIR_IS_PROPER: usize = 1;
    pub const FLAG_HAS_MAPQ: usize = 2;
    pub const FLAG_HAS_SCORE: usize = 3;

    /// Creates a new flags object with the given number of alignments.
    /// The flags are initialized to zero (false).
    pub fn new(num_alignments: usize) -> Self {
        let bits = vec![0; (num_alignments * Self::NUM_FLAGS + 7) / 8];
        Self { bits }
    }

    /// Sets the given flag.
    ///
    /// # Arguments
    ///
    /// * `index`: The index of the alignment.
    /// * `flag`: The flag to set.
    /// * `value`: The value to set the flag to.
    pub fn set(&mut self, index: usize, flag: usize, value: bool) {
        let bit = index * Self::NUM_FLAGS + flag;
        if value {
            self.bits[bit / 8] |= 1 << (bit % 8);
        } else {
            self.bits[bit / 8] &= !(1 << (bit % 8));
        }
    }

    /// Gets the value of the given flag.
    ///
    /// # Arguments
    ///
    /// * `index`: The index of the alignment.
    /// * `flag`: The flag to get.
    pub fn get(&self, index: usize, flag: usize) -> bool {
        let bit = index * Self::NUM_FLAGS + flag;
        (self.bits[bit / 8] >> (bit % 8)) & 0x01 != 0
    }

    /// Returns the number of bytes used to store the flags.
    pub fn bytes(&self) -> usize {
        self.bits.len()
    }
}

impl From<Vec<u8>> for Flags {
    fn from(bits: Vec<u8>) -> Self {
        Self { bits }
    }
}

impl AsRef<[u8]> for Flags {
    fn as_ref(&self) -> &[u8] {
        &self.bits
    }
}

//-----------------------------------------------------------------------------

// FIXME: document, example, tests
// An encoded block of alignments.
#[derive(Debug, Clone)]
pub struct AlignmentBlock {
    pub min_node: Option<usize>,
    pub max_node: Option<usize>,
    pub alignments: usize,
    pub gbwt_starts: Vec<u8>,
    pub names: Vec<u8>,
    pub quality_strings: Vec<u8>,
    pub difference_strings: Vec<u8>,
    pub flags: Flags,
    pub numbers: Vec<u8>,
}

impl AlignmentBlock {
    // TODO: Make this a parameter?
    /// Compression level for Zstandard.
    pub const COMPRESSION_LEVEL: i32 = 7;

    // FIXME: Error handling. Especially if there is a mix of aligned and unaligned reads.
    pub fn new(alignments: &[Alignment], index: &GBWT, first_id: usize) -> Self {
        let min_node = alignments.iter().map(|aln| aln.min_node()).min().flatten();
        let max_node = alignments.iter().map(|aln| aln.max_node()).max().flatten();

        let mut gbwt_starts = ByteCode::new();
        let mut names: Vec<u8> = Vec::new();
        let mut quality_strings: Vec<u8> = Vec::new();
        let mut difference_strings = RLE::with_sigma(Difference::NUM_TYPES);
        let mut flags = Flags::new(alignments.len());
        let mut numbers = ByteCode::new();

        let base_node = min_node.unwrap_or(0);
        for (i, aln) in alignments.iter().enumerate() {
            // GBWT start as (node - min_node, offset).
            if let Some(start) = index.start(support::encode_path(first_id + i, Orientation::Forward)) {
                gbwt_starts.write(start.node - base_node);
                gbwt_starts.write(start.offset);
            } else {
                assert!(gbwt_starts.is_empty(), "Line {}: Unaligned read in a block with aligned reads", first_id + i);
            }

            // Read/pair names as 0-terminated strings.
            names.extend_from_slice(aln.name.as_bytes());
            names.push(0);
            if let Some(pair) = &aln.pair {
                names.extend_from_slice(pair.name.as_bytes());
            }
            names.push(0);

            // Quality strings as 0-terminated strings.
            if let Some(quality) = &aln.base_quality {
                quality_strings.extend_from_slice(quality);
            }
            quality_strings.push(0);

            // Difference strings using a custom encoder.
            aln.encode_difference_into(&mut difference_strings);

            // Flags as an ad hoc bitvector.
            if let Some(pair) = &aln.pair {
                if pair.is_next {
                    flags.set(i, Flags::FLAG_PAIR_IS_NEXT, true);
                }
                if pair.is_proper {
                    flags.set(i, Flags::FLAG_PAIR_IS_PROPER, true);
                }
            }
            if aln.mapq.is_some() {
                flags.set(i, Flags::FLAG_HAS_MAPQ, true);
            }
            if aln.score.is_some() {
                flags.set(i, Flags::FLAG_HAS_SCORE, true);
            }

            // Numbers with a custom encoder.
            aln.encode_numbers_into(&mut numbers);
        }

        // FIXME: error handling
        // Compress the names and quality strings.
        let mut name_encoder = ZstdEncoder::new(Vec::new(), Self::COMPRESSION_LEVEL).unwrap();
        name_encoder.write_all(&names).unwrap();
        names = name_encoder.finish().unwrap();
        let mut quality_encoder = ZstdEncoder::new(Vec::new(), Self::COMPRESSION_LEVEL).unwrap();
        quality_encoder.write_all(&quality_strings).unwrap();
        quality_strings = quality_encoder.finish().unwrap();

        Self {
            min_node,
            max_node,
            alignments: alignments.len(),
            gbwt_starts: Vec::from(gbwt_starts),
            names,
            quality_strings,
            difference_strings: Vec::from(difference_strings),
            flags,
            numbers: Vec::from(numbers),
        }
    }

    pub fn len(&self) -> usize {
        self.alignments
    }

    pub fn is_empty(&self) -> bool {
        self.alignments == 0
    }

    fn decode_gbwt_starts(&self) -> Result<Vec<Pos>, String> {
        let mut result = Vec::new();
        if self.min_node.is_none() {
            return Ok(result);
        }
        let base_node = self.min_node.unwrap();
        let mut decoder = ByteCodeIter::new(&self.gbwt_starts[..]);
        for i in 0..self.len() {
            let start = decoder.next().ok_or(
                format!("Missing GBWT start for alignment {}", i)
            )?;
            let offset = decoder.next().ok_or(
                format!("Missing GBWT offset for alignment {}", i)
            )?;
            result.push(Pos::new(start + base_node, offset));
        }
        Ok(result)
    }

    fn zstd_decompress(&self, data: &[u8]) -> Result<Vec<u8>, String> {
        let mut decoder = ZstdDecoder::new(data).map_err(|err| format!("Zstd decompression error: {}", err))?;
        let mut buffer = Vec::new();
        decoder.read_to_end(&mut buffer).map_err(|err| format!("Zstd decompression error: {}", err))?;
        Ok(buffer)
    }

    fn decode_pair(&self, i: usize, names: &[&[u8]]) -> Option<PairedRead> {
        if names[2 * i + 1].is_empty() {
            return None;
        }
        let name = String::from_utf8_lossy(names[2 * i]).to_string();
        let is_next = self.flags.get(i, Flags::FLAG_PAIR_IS_NEXT);
        let is_proper = self.flags.get(i, Flags::FLAG_PAIR_IS_PROPER);
        Some(PairedRead { name, is_next, is_proper })
    } 

    // FIXME: What if we have GBZ-base instead of GBWT?
    // FIXME: Maybe don't collect the slices?
    pub fn decode(&self) -> Result<Vec<Alignment>, String> {
        let mut result = Vec::with_capacity(self.len());

        // Decompress the GBWT starting positions.
        let gbwt_starts = self.decode_gbwt_starts()?;

        // Decompress the names.
        let name_buffer = self.zstd_decompress(&self.names[..])?;
        let names = name_buffer.split(|&c| c == 0).collect::<Vec<_>>();
        if names.len() != 2 * self.len() + 1 {
            // We have an empty slice at the end, as the buffer is 0-terminated.
            return Err(format!("Split the names into {} slices in a block of {} alignments", names.len(), self.len()));
        }

        // Decompress the quality strings.
        let quality_buffer = self.zstd_decompress(&self.quality_strings[..])?;
        let quality_strings = quality_buffer.split(|&c| c == 0).collect::<Vec<_>>();
        if quality_strings.len() != self.len() + 1 {
            // We have an empty slice at the end, as the buffer is 0-terminated.
            return Err(format!("Split the quality strings into {} slices in a block of {} alignments", quality_strings.len(), self.len()));
        }

        // Decompress the difference strings on the fly.
        let mut difference_decoder = RLEIter::with_sigma(&self.difference_strings[..], Difference::NUM_TYPES);

        // We decode the numbers on the fly, as there are too many special cases.
        let mut number_decoder = ByteCodeIter::new(&self.numbers[..]);

        // TODO: Refactor this. Maybe a function for decoding the numbers?
        // Build the alignments.
        for i in 0..self.len() {
            let name = String::from_utf8_lossy(names[2 * i]).to_string();

            let path = if gbwt_starts.is_empty() {
                TargetPath::Path(Vec::new())
            } else {
                TargetPath::StartPosition(gbwt_starts[i])
            };

            let base_quality = if quality_strings[i].is_empty() {
                None
            } else {
                Some(Vec::from(quality_strings[i]))
            };

            let difference = Alignment::decode_difference_from(&self.difference_strings, &mut difference_decoder)?;
            let include_redundant = difference.is_empty();

            // Now we can decode the numbers.
            let (mut seq_interval, mut seq_len) = Alignment::decode_coordinates(&mut number_decoder, include_redundant).ok_or(
                format!("Missing query sequence coordinates for alignment {}", i)
            )?;
            let (mut path_interval, mut path_len) = Alignment::decode_coordinates(&mut number_decoder, include_redundant).ok_or(
                format!("Missing target path coordinates for alignment {}", i)
            )?;
            let mut matches = if include_redundant {
                number_decoder.next().ok_or(format!("Missing number of matches for alignment {}", i))?
            } else { 0 };
            let mut edits = if include_redundant {
                number_decoder.next().ok_or(format!("Missing number of edits for alignment {}", i))?
            } else { 0 };
            let mapq = if self.flags.get(i, Flags::FLAG_HAS_MAPQ) {
                Some(number_decoder.next().ok_or(format!("Missing mapping quality for alignment {}", i))?)
            } else { None };
            let score = if self.flags.get(i, Flags::FLAG_HAS_SCORE) {
                Some(Alignment::decode_signed(&mut number_decoder).ok_or(format!("Missing alignment score for alignment {}", i))?)
            } else { None };

            // And update the numbers if we have a difference string.
            if !difference.is_empty() {
                let (query_len, target_len, num_matches, num_edits) = Difference::stats(&difference);
                seq_interval.end = seq_interval.start + query_len;
                seq_len += query_len;
                path_interval.end = path_interval.start + target_len;
                path_len += target_len;
                matches = num_matches;
                edits = num_edits;
            }

            let pair = self.decode_pair(i, &names);
            let optional = Vec::new();

            result.push(Alignment {
                name,
                seq_len, seq_interval,
                path,
                path_len, path_interval,
                matches, edits,
                mapq, score,
                base_quality,
                difference,
                pair,
                optional,
            });
        }

        Ok(result)
    }
}

//-----------------------------------------------------------------------------
