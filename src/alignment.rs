//! Structures for representing sequence to graph alignments.
//!
//! An [`Alignment`] object represents the alignment of a query sequence to a target path in a graph.
//! It corresponds to a single line in a GAF file.
//!
//! An [`AlignmentBlock`] object represents an ordered collection of alignments, using column-based compression.
//! It corresponds to a row in table `Alignments` in a [`crate::GAFBase`].
//!
//! The GAF format is a text-based format for representing sequence alignments to a graph.
//! See [the specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) for an overview.
//! Some details are better documented in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html#10).
//!
//! # Assumptions
//!
//! Some details of the GAF format are not well documented, and this implementation makes some assumptions.
//!
//! * Softclips are included either in both the query interval and the difference string or in neither.
//! * If a difference string is present, it overrides query/path interval ends, the number of matches, and alignment block length.
//! * Alignment block length is the sum of matches, mismatches, insertions, and deletions.
//!
//! # Supported optional fields
//!
//! * `AS:i`: Alignment score.
//! * `bq:Z`: Base quality values for the query sequence.
//! * `cs:Z`: Difference string for the alignment.
//! * `pd:b`: Properly paired flag.
//! * `fn:Z`: Name of the next read in the pair.
//! * `fp:Z`: Name of the previous read in the pair.
//!
//! Only matches, mismatches, insertions, and deletions are supported in the difference string.
//! Matches are represented as a match length rather than as the matching sequence.
//!
//! The order of optional fields is not preserved.
//! Compression with [`AlignmentBlock`] discards unsupported optional fields.

use crate::{Subgraph, Mapping, Difference};
use crate::formats::{self, TypedField};
use crate::utils;

use std::io::{Read, Write};
use std::ops::Range;
use std::sync::Arc;
use std::{cmp, str};

use zstd::stream::Encoder as ZstdEncoder;
use zstd::stream::Decoder as ZstdDecoder;

use gbwt::{GBWT, Orientation, Pos};
use gbwt::support::{self, ByteCode, ByteCodeIter, RLE, Run, RLEIter};

#[cfg(test)]
mod tests;

pub mod mapping;

//-----------------------------------------------------------------------------

// TODO: add validate() against a graph, a subgraph, or a read set
/// An alignment between a query sequence and a target path in a graph.
///
/// This object corresponds either to a line in a GAF file or to a row in table `Alignments` in [`crate::GAFBase`].
/// When the alignment is built from a GAF line, the target path is stored explicitly.
/// For alignments stored in a database, only the GBWT starting position is stored.
/// See [`TargetPath`] for details.
///
/// A GAF line can be converted to an `Alignment` object with [`Alignment::from_gaf`].
/// An `Alignment` object can be serialized as a GAF line with [`Alignment::to_gaf`].
/// The conversion can be lossy (see [`crate::alignment`] for details).
///
/// # Examples
///
/// ```
/// use gbz_base::Alignment;
/// use gbwt::Orientation;
/// use gbwt::support;
///
/// // Unnamed empty sequence.
/// let empty = Alignment::new();
/// assert!(empty.name.is_empty());
/// assert_eq!(empty.seq_len, 0);
/// assert!(empty.is_unaligned());
///
/// // Construct a GAF line.
/// let name = "query";
/// let seq_len = 7;
/// let seq_interval = 0..6;
/// let path = ">1<2>3";
/// let path_len = 8;
/// let path_interval = 2..8;
/// let matches = 5;
/// let edits = 1;
/// let mapq = 60;
/// let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
///     name,
///     seq_len, seq_interval.start, seq_interval.end,
///     '+', path,
///     path_len, path_interval.start, path_interval.end,
///     matches, matches + edits, mapq
/// );
///
/// // Create an alignment from the line.
/// let aln = Alignment::from_gaf(line.as_bytes()).unwrap();
/// assert_eq!(aln.name, name);
/// assert_eq!(aln.seq_len, seq_len);
/// assert_eq!(aln.seq_interval, seq_interval);
/// assert!(aln.has_target_path());
/// let path = vec!(
///     support::encode_node(1, Orientation::Forward),
///     support::encode_node(2, Orientation::Reverse),
///     support::encode_node(3, Orientation::Forward)
/// );
/// assert_eq!(aln.target_path(), Some(path.as_ref()));
/// assert_eq!(aln.path_len, path_len);
/// assert_eq!(aln.path_interval, path_interval);
/// assert_eq!(aln.matches, matches);
/// assert_eq!(aln.edits, edits);
/// assert_eq!(aln.mapq, Some(mapq));
///
/// // All optional fields are missing.
/// assert!(aln.score.is_none());
/// assert!(aln.base_quality.is_empty());
/// assert!(aln.difference.is_empty());
/// assert!(aln.pair.is_none());
/// assert!(aln.optional.is_empty());
///
/// // Convert back to GAF.
/// let _query_sequence =   b"GATTACA"; // Not used yet.
/// let target_sequence = b"GAGATCAC".to_vec();
/// let from_aln = aln.to_gaf(&target_sequence);
/// assert_eq!(&from_aln, line.as_bytes());
/// ```
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
    pub base_quality: Vec<u8>,
    /// Difference string, or an empty vector if not present.
    pub difference: Vec<Difference>,
    /// Information about the paired alignment.
    pub pair: Option<PairedRead>,
    /// Optional typed fields that have not been interpreted.
    pub optional: Vec<TypedField>,
}

impl Default for Alignment {
    fn default() -> Self {
        Alignment {
            name: String::new(),
            seq_len: 0,
            seq_interval: 0..0,
            path: TargetPath::Path(Vec::new()),
            path_len: 0,
            path_interval: 0..0,
            matches: 0,
            edits: 0,
            mapq: None,
            score: None,
            base_quality: Vec::new(),
            difference: Vec::new(),
            pair: None,
            optional: Vec::new(),
        }
    }
}

//-----------------------------------------------------------------------------

/// Construction; conversions between GAF lines and `Alignment` objects.
impl Alignment {
    // Number of mandatory fields in a GAF line.
    const MANDATORY_FIELDS: usize = 12;

    // Placeholder value for a missing mapping quality.
    const MISSING_MAPQ: usize = 255;

    // The field is empty and the value is missing; typically used with unaligned sequences.
    const MISSING_VALUE: [u8; 1] = [b'*'];

    /// Creates an empty alignment.
    pub fn new() -> Self {
        Alignment::default()
    }

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
        let mut base_quality = Vec::new();
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
                    if !base_quality.is_empty() {
                        return Err(String::from("Multiple base quality fields"));
                    }
                    base_quality = value;
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

    // TODO: Should we take the sequence for the full target path instead?
    /// Converts the alignment to a GAF line.
    ///
    /// If the target path is stored as a GBWT starting position, it will be missing (`*`).
    /// The path can be set with [`Alignment::set_target_path`] or extracted from a GBWT index with [`Alignment::extract_target_path`].
    /// The returned line does not end with an endline character.
    ///
    /// A target sequence is necessary for reconstructing the difference string, if present.
    /// It must correspond only to `path_interval` of the target path.
    pub fn to_gaf(&self, target_sequence: &[u8]) -> Vec<u8> {
        let mut result = Vec::new();

        result.extend_from_slice(self.name.as_bytes());

        // Coordinates in the query sequence.
        result.push(b'\t');
        utils::append_usize(&mut result, self.seq_len);
        result.push(b'\t');
        utils::append_usize(&mut result, self.seq_interval.start);
        result.push(b'\t');
        utils::append_usize(&mut result, self.seq_interval.end);

        // Target path; always in the forward orientation.
        result.push(b'\t');
        result.push(b'+');
        result.push(b'\t');
        match &self.path {
            TargetPath::Path(path) => {
                if path.is_empty() {
                    result.extend_from_slice(&Self::MISSING_VALUE);
                } else {
                    formats::append_walk(&mut result, path);
                }
            },
            TargetPath::StartPosition(_) => {
                result.extend_from_slice(&Self::MISSING_VALUE);
            },
        }

        // Coordinates in the target path.
        result.push(b'\t');
        utils::append_usize(&mut result, self.path_len);
        result.push(b'\t');
        utils::append_usize(&mut result, self.path_interval.start);
        result.push(b'\t');
        utils::append_usize(&mut result, self.path_interval.end);

        // Alignment statistics.
        result.push(b'\t');
        utils::append_usize(&mut result, self.matches);
        result.push(b'\t');
        utils::append_usize(&mut result, self.matches + self.edits);
        result.push(b'\t');
        let mapq = self.mapq.unwrap_or(Self::MISSING_MAPQ);
        utils::append_usize(&mut result, mapq);

        // Known optional fields.
        if let Some(score) = self.score {
            let field = TypedField::Int([b'A', b'S'], score);
            field.append_to(&mut result, true);
        }
        if !self.base_quality.is_empty() {
            TypedField::append_string(&mut result, [b'b', b'q'], &self.base_quality, true);
        }
        if !self.difference.is_empty() {
            // TODO: Can we write the difference string directly?
            let field = TypedField::String([b'c', b's'], Difference::to_bytes(&self.difference, target_sequence));
            field.append_to(&mut result, true);
        }
        if let Some(pair) = &self.pair {
            if pair.is_next {
                TypedField::append_string(&mut result, [b'f', b'n'], pair.name.as_bytes(), true);
            } else {
                TypedField::append_string(&mut result, [b'f', b'p'], pair.name.as_bytes(), true);
            }
            let field = TypedField::Bool([b'p', b'd'], pair.is_proper);
            field.append_to(&mut result, true);
        }

        // Other optional fields.
        for field in self.optional.iter() {
            field.append_to(&mut result, true);
        }

        result
    }
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
    // * Query sequence: `seq_len`, `seq_interval.start`, `seq_interval.end`.
    // * Target path: `path_len`, `path_interval.start`, `path_interval.end`.
    // * Alignment statistics: `matches`, `edits`, `mapq`, `score`.
    //
    // The coordinates are normalized so that `start <= end <= len`.
    // If `perfect_alignment` is set, this is assumed to be a perfect alignment of the expected length.
    // In that case, numbers that can be derived from the length are not stored.
    // If a difference string is present, numbers that can be derived from it are not stored.
    //
    // Mapping quality and alignment score are stored if they are present.
    // Their presence should be stored separately.
    fn encode_numbers_into(&self, encoder: &mut ByteCode, perfect_alignment: bool) {
        let include_redundant = self.difference.is_empty();

        // Coordinates.
        if perfect_alignment {
            // We can also derive target interval length from the sequence length.
            Self::encode_coordinates(self.path_interval.clone(), self.path_len, false, encoder);
        } else {
            Self::encode_coordinates(self.seq_interval.clone(), self.seq_len, include_redundant, encoder);
            Self::encode_coordinates(self.path_interval.clone(), self.path_len, include_redundant, encoder);
        }

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

// TODO: Add an operation for reconstructing the query sequence.
/// Operations on the Alignment object.
impl Alignment {
    /// Returns `true` if the read is unaligned.
    ///
    /// An aligned read has a non-empty query interval aligned to a non-empty target interval.
    /// NOTE: An empty sequence is by definition unaligned.
    pub fn is_unaligned(&self) -> bool {
        self.seq_interval.is_empty() || self.path_interval.is_empty()
    }

    /// Returns `true` if this is a perfect alignment of the entire read.
    pub fn is_perfect(&self) -> bool {
        self.seq_len > 0 &&
            self.seq_interval.start == 0 &&
            self.seq_interval.end == self.seq_len &&
            self.edits == 0
    }

    /// Returns the minimum GBWT node identifier in the target path, or [`None`] if there is no path.
    pub fn min_handle(&self) -> Option<usize> {
        match self.path {
            TargetPath::Path(ref path) => path.iter().copied().min(),
            TargetPath::StartPosition(_) => None,
        }
    }

    /// Returns the maximum GBWT node identifier in the target path, or [`None`] if there is no path.
    pub fn max_handle(&self) -> Option<usize> {
        match self.path {
            TargetPath::Path(ref path) => path.iter().copied().max(),
            TargetPath::StartPosition(_) => None,
        }
    }

    /// Returns `true` if the target path is stored explicitly.
    pub fn has_target_path(&self) -> bool {
        matches!(self.path, TargetPath::Path(_))
    }

    /// Returns `true` if the target path is stored explicitly and is non-empty.
    pub fn has_non_empty_target_path(&self) -> bool {
        match &self.path {
            TargetPath::Path(path) => !path.is_empty(),
            TargetPath::StartPosition(_) => false,
        }
    }

    /// Returns the target path if it is stored explicitly.
    pub fn target_path(&self) -> Option<&[usize]> {
        match &self.path {
            TargetPath::Path(path) => Some(path),
            TargetPath::StartPosition(_) => None,
        }
    }

    // TODO: Should this update an existing target path?
    /// Sets the target path from the GBWT index if it is not already present.
    pub fn extract_target_path(&mut self, index: &GBWT) {
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

    /// Sets the given path as the target path.
    pub fn set_target_path(&mut self, path: Vec<usize>) {
        self.path = TargetPath::Path(path);
    }

    // FIXME: Make this work with unaligned reads represented as insertions.
    /// Returns an iterator over the alignment as a sequence of mappings.
    ///
    /// Returns [`None`] if a valid iterator cannot be built.
    /// The iterator requires a difference string and an explicitly stored target path.
    /// It may stop early if the alignment is invalid.
    ///
    /// The iterator needs a function that provides the sequence length for each node.
    /// This function may be based on [`gbwt::GBZ`], [`Subgraph`], or [`crate::ReadSet`].
    pub fn iter<'a>(&'a self, sequence_len: Arc<dyn Fn(usize) -> Option<usize> + 'a>) -> Option<AlignmentIter<'a>> {
        if !self.has_target_path() {
            return None;
        }
        let target_path = self.target_path().unwrap();
        if target_path.is_empty() != self.is_unaligned() {
            return None;
        }
        if self.difference.is_empty() != self.is_unaligned() {
            return None;
        }

        let mut iter = AlignmentIter {
            parent: self,
            sequence_len,
            seq_offset: self.seq_interval.start,
            path_offset: self.path_interval.start,
            path_vec_offset: 0,
            path_node_offset: self.path_interval.start,
            diff_vec_offset: 0,
            diff_op_offset: 0,
        };
        if self.is_unaligned() {
            return Some(iter);
        }

        // In some edge cases, there may be unused nodes at the start of the target path.
        let mut handle = *target_path.get(iter.path_vec_offset)?;
        let mut node_len = (*iter.sequence_len)(handle)?;
        while iter.path_node_offset >= node_len {
            iter.path_node_offset -= node_len;
            iter.path_vec_offset += 1;
            handle = *target_path.get(iter.path_vec_offset)?;
            node_len = (*iter.sequence_len)(handle)?;
        }

        Some(iter)
    }

    // Creates an empty fragment of this alignment.
    fn empty_fragment(&self, fragment_id: usize) -> Alignment {
        let mut aln = Alignment {
            name: self.name.clone(),
            seq_len: self.seq_len,
            seq_interval: 0..0,
            path: TargetPath::Path(Vec::new()),
            path_len: 0,
            path_interval: 0..0,
            matches: self.matches,
            edits: self.edits,
            mapq: self.mapq,
            score: self.score,
            base_quality: self.base_quality.clone(),
            difference: Vec::new(),
            pair: self.pair.clone(),
            optional: self.optional.clone(),
        };
        aln.optional.push(TypedField::Int([b'f', b'i'], fragment_id as isize));

        aln
    }

    // Extends the given fragment with the next mapping of the same original alignment.
    fn extend(&mut self, mapping: Mapping) -> Result<(), String> {
        // Start by updating the difference string.
        if self.seq_interval.is_empty() && self.path_interval.is_empty() {
            // If we started with a gap, `is_unaligned()` would be true.
            if !self.difference.is_empty() {
                return Err(String::from("Cannot extend an unaligned fragment with a non-empty difference string"));
            }
            self.difference.push(mapping.edit().clone());
        } else {
            if self.difference.is_empty() {
                return Err(String::from("Cannot extend an aligned fragment without a difference string"));
            }
            if !self.difference.last_mut().unwrap().try_merge(mapping.edit()) {
                self.difference.push(mapping.edit().clone());
            }
        }

        // Update the query sequence.
        if mapping.seq_interval().end > self.seq_len {
            return Err(String::from("Cannot extend a fragment beyond the query sequence length"));
        }
        if self.seq_interval.is_empty() {
            // Either the first mapping or we only have deletions so far.
            self.seq_interval = mapping.seq_interval().clone();
        } else {
            if mapping.seq_interval().start != self.seq_interval.end {
                return Err(String::from("Cannot append a non-contiguous query interval"));
            }
            self.seq_interval.end = mapping.seq_interval().end;
        }

        // TODO: enum?
        // Determine how we are extending the alignment.
        let target_path = self.target_path().ok_or(
            "Cannot extend a fragment without an explicit target path"
        )?;
        let last_node = target_path.last().copied();
        let path_left = self.path_len.saturating_sub(self.path_interval.end);
        let reverse_offset = mapping.node_len().saturating_sub(mapping.node_interval().start);
        let first_mapping = last_node.is_none();
        let continues_in_same_node = Some(mapping.handle()) == last_node && reverse_offset == path_left;
        let starts_a_new_node = path_left == 0 && mapping.is_at_start();

        // Update the target path.
        if first_mapping {
            if let TargetPath::Path(path) = &mut self.path {
                path.push(mapping.handle());
            } else {
                unreachable!();
            }
            self.path_len += mapping.node_len();
            self.path_interval = mapping.node_interval().clone();
        } else if starts_a_new_node {
            if let TargetPath::Path(path) = &mut self.path {
                path.push(mapping.handle());
            } else {
                unreachable!();
            }
            self.path_len += mapping.node_len();
            self.path_interval.end += mapping.target_len();
        } else if continues_in_same_node {
            self.path_interval.end += mapping.target_len();
        } else {
            return Err(String::from("Cannot append a non-contiguous target interval"));
        }

        Ok(())
    }

    /// Clips the alignment into fragments that are fully contained in the given subgraph.
    ///
    /// Returns an empty vector if the read is unaligned or lacks an explicit non-empty target path or a difference string.
    /// There will be one alignment fragment for every maximal subpath in the subgraph.
    /// The fragments have the same name, pair, optional fields, and statistics as the original alignment.
    /// Only the aligned intervals and difference strings depend on the fragment.
    /// Fragments of the same alignment are identified by a fragment index stored as an optional field `fi:i`.
    ///
    /// Clipping requires a function that returns the sequence length for the node with the given handle.
    /// This function may be based on [`gbwt::GBZ`], [`Subgraph`], or [`crate::ReadSet`].
    ///
    /// # Errors
    ///
    /// Returns an error if an [`AlignmentIter`] cannot be created.
    /// May return an error if the alignment is invalid and the iterator returns non-consecutive mappings.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbwt::{GBZ, Orientation};
    /// use gbwt::support;
    /// use gbz_base::{Alignment, Difference, Subgraph};
    /// use gbz_base::utils;
    /// use simple_sds::serialize;
    /// use std::sync::Arc;
    ///
    /// let gbz_filename = utils::get_test_data("micb-kir3dl1.gbz");
    /// let graph: GBZ = serialize::load_from(&gbz_filename).unwrap();
    ///
    /// // Create a perfect alignment with node lengths of 68 bp, 1 bp, and 6 bp.
    /// let mut aln = Alignment::new();
    /// aln.seq_len = 50;
    /// aln.seq_interval = 0..50;
    /// let target_path = vec![
    ///     support::encode_node(13, Orientation::Forward),
    ///     support::encode_node(14, Orientation::Forward),
    ///     support::encode_node(16, Orientation::Forward),
    /// ];
    /// aln.set_target_path(target_path.clone());
    /// aln.path_len = 75;
    /// aln.path_interval = 20..70;
    /// aln.matches = 50;
    /// aln.difference.push(Difference::Match(50));
    ///
    /// // Create a subgraph without the middle node.
    /// let mut subgraph = Subgraph::new();
    /// let _ = subgraph.add_node_from_gbz(&graph, 13);
    /// let _ = subgraph.add_node_from_gbz(&graph, 16);
    ///
    /// // Clip the alignment to the subgraph.
    /// let sequence_len = Arc::new(|handle| {
    ///     let (node_id, _) = support::decode_node(handle);
    ///     graph.sequence_len(node_id)
    /// });
    /// let result = aln.clip(&subgraph, sequence_len.clone());
    /// assert!(result.is_ok());
    /// let clipped = result.unwrap();
    ///
    /// // We should have two fragments.
    /// assert_eq!(clipped.len(), 2);
    /// assert_eq!(clipped[0].seq_interval, 0..48);
    /// assert_eq!(clipped[0].target_path().unwrap(), &target_path[0..1]);
    /// assert_eq!(clipped[0].path_interval, 20..68);
    /// assert_eq!(clipped[0].difference, vec![Difference::Match(48)]);
    /// assert_eq!(clipped[1].seq_interval, 49..50);
    /// assert_eq!(clipped[1].target_path().unwrap(), &target_path[2..3]);
    /// assert_eq!(clipped[1].path_interval, 0..1);
    /// assert_eq!(clipped[1].difference, vec![Difference::Match(1)]);
    /// ```
    pub fn clip<'a>(&self, subgraph: &Subgraph, sequence_len: Arc<impl Fn(usize) -> Option<usize> + 'a>) -> Result<Vec<Alignment>, String> {
        let mut result = Vec::new();
        if self.is_unaligned() || !self.has_non_empty_target_path() || self.difference.is_empty() {
            return Ok(result);
        }

        let mut aln: Option<Alignment> = None; // The alignment we are currently building.
        let iter = self.iter(sequence_len);
        if iter.is_none() {
            eprintln!("Warning: Cannot build an alignment iterator for {}", self.name);
            return Ok(result);
        }
        let iter = iter.unwrap();
        for mapping in iter {
            if !subgraph.has_handle(mapping.handle()) {
                if let Some(prev) = aln {
                    result.push(prev);
                    aln = None;
                }
                continue;
            }

            if let Some(aln) = &mut aln {
                aln.extend(mapping)?;
            } else {
                let mut curr = self.empty_fragment(result.len() + 1);
                curr.extend(mapping)?;
                aln = Some(curr);
            }
        }
        if let Some(aln) = aln {
            result.push(aln);
        }

        Ok(result)
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

/// An iterator over an [`Alignment`] as a sequence of [`Mapping`] objects.
///
/// This iterator assumes that the alignment is valid, and it may stop early if it is not.
/// It requires that the alignment stores the target path explicitly and has a difference string.
/// Insertions at node boundaries are assigned to the left.
///
/// # Examples
///
/// ```
/// use gbz_base::{Alignment, Mapping, Difference};
/// use std::sync::Arc;
///
/// // Construct an alignment.
/// let mut aln = Alignment::new();
/// aln.seq_len = 15;
/// aln.seq_interval = 0..15;
/// aln.set_target_path(vec![1, 2, 3]);
/// aln.path_len = 15;
/// aln.path_interval = 0..15;
/// aln.matches = 14;
/// aln.edits = 1;
/// aln.difference = vec![
///     Difference::Match(6),
///     Difference::Mismatch(b'A'),
///     Difference::Match(8),
/// ];
///
/// // We use this in place of a graph for node lengths.
/// let node_lengths = vec![0, 4, 6, 5];
/// let sequence_len = Arc::new(|handle| node_lengths.get(handle).copied());
///
/// // This is what we expect from the iterator.
/// let truth = vec![
///     Mapping::new(0, 1, node_lengths[1], 0, Difference::Match(4)),
///     Mapping::new(4, 2, node_lengths[2], 0, Difference::Match(2)),
///     Mapping::new(6, 2, node_lengths[2], 2, Difference::Mismatch(b'A')),
///     Mapping::new(7, 2, node_lengths[2], 3, Difference::Match(3)),
///     Mapping::new(10, 3, node_lengths[3], 0, Difference::Match(5)),
/// ];
///
/// // Iterator creation may fail if the alignment is invalid.
/// let iter = aln.iter(sequence_len);
/// assert!(iter.is_some());
/// let iter = iter.unwrap();
/// assert!(iter.eq(truth.iter().cloned()));
/// ```
#[derive(Clone)]
pub struct AlignmentIter<'a> {
    parent: &'a Alignment,
    sequence_len: Arc<dyn Fn(usize) -> Option<usize> + 'a>,
    // Position in the query sequence.
    seq_offset: usize,
    // Position in the target sequence.
    path_offset: usize,
    // Position in the target path.
    path_vec_offset: usize,
    // Position within the current node in the target path.
    path_node_offset: usize,
    // Position in the difference string.
    diff_vec_offset: usize,
    // Target position in the current difference operation.
    diff_op_offset: usize,
}

impl<'a> Iterator for AlignmentIter<'a> {
    type Item = Mapping;

    fn next(&mut self) -> Option<Self::Item> {
        if self.diff_vec_offset >= self.parent.difference.len() {
            return None;
        }

        // First find the edit that forms the basis of the mapping.
        let full_edit = self.parent.difference.get(self.diff_vec_offset)?;
        let edit_left = full_edit.target_len() - self.diff_op_offset;

        // Then find the target node. We may need to advance to the next node we are
        // at the end of the node and we have an edit with non-zero target length.
        // We prefer mapping insertions to the end of a node rather than to the start
        // of the next node, as the next node may not exist.
        let target_path = self.parent.target_path()?;
        let mut handle = *target_path.get(self.path_vec_offset)?;
        let mut node_len = (*self.sequence_len)(handle)?;
        if self.path_node_offset >= node_len && edit_left > 0 {
            self.path_vec_offset += 1;
            self.path_node_offset = 0;
            handle = *target_path.get(self.path_vec_offset)?;
            node_len = (*self.sequence_len)(handle)?;
        }
        let node_left = node_len - self.path_node_offset;

        // Take the offsets for the mapping before we update them.
        let seq_offset = self.seq_offset;
        let node_offset = self.path_node_offset;

        // Now determine the actual edit.
        let edit_len = cmp::min(edit_left, node_left);
        let edit = match full_edit {
            Difference::Match(_) => Difference::Match(edit_len),
            Difference::Mismatch(base) => Difference::Mismatch(*base),
            Difference::Insertion(seq) => {
                // We can always take the full insertion within the current node.
                Difference::Insertion(seq.clone())
            },
            Difference::Deletion(_) => Difference::Deletion(edit_len),
            Difference::End => {
                self.diff_vec_offset = self.parent.difference.len();
                return None;
            },
        };

        // And then advance the iterator.
        self.seq_offset += edit.query_len();
        self.path_offset += edit.target_len();
        // We do not advance the node, as the next edit may be an insertion.
        self.path_node_offset += edit.target_len();
        self.diff_op_offset += edit.target_len();
        if self.diff_op_offset >= full_edit.target_len() {
            self.diff_vec_offset += 1;
            self.diff_op_offset = 0;
        }

        Some(Mapping::new(seq_offset, handle, node_len, node_offset, edit))
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
    pub const NUM_FLAGS: usize = 5;

    // Flags expressed as a bit offset.
    // TODO: As enum?

    /// Flag: Pair is the next fragment in the read pair.
    pub const FLAG_PAIR_IS_NEXT: usize = 0;

    /// Flag: The the alignment is properly paired.
    pub const FLAG_PAIR_IS_PROPER: usize = 1;

    /// Flag: A mapping quality score is present.
    pub const FLAG_HAS_MAPQ: usize = 2;

    /// Flag: An alignment score is present.
    pub const FLAG_HAS_SCORE: usize = 3;

    /// Flag: This is a perfect alignment of the entire read, and read length is as expected.
    pub const FLAG_PERFECT_ALIGNMENT: usize = 4;

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

// TODO: expected alignment score for perfect alignments? Or just a scoring model?
/// An encoded block of [`Alignment`] objects.
///
/// This is a compressed representation of multiple sequences aligned to the same graph.
/// With short reads, a properly chosen block size enables both column-based compression and random access to the alignments.
/// Long reads do not benefit from column-based compression, as the amount of metadata is insignificant.
/// Reasonable block sizes could be 1000 alignments for short reads and 10 for long reads.
///
/// For best results, node identifiers should approximate a topological order in the graph.
/// The range of node identifiers in a path is typically proportional to the length of the path.
/// If the alignments are sorted by (min id, max id), a block will then consist of alignments that are close to each other in the graph.
///
/// # Notes
///
/// * Target paths must be stored separately (e.g. using a GBWT-based encoding).
/// * A block must contain either aligned reads or unaligned reads, but not both.
///
/// # Examples
///
/// ```
/// use gbz_base::{Alignment, AlignmentBlock};
/// use gbz_base::utils;
/// use gbwt::GBWT;
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// // We need a GBWT index of the paths for compressing alignments.
/// let gbwt_filename = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
/// let index: GBWT = serialize::load_from(&gbwt_filename).unwrap();
///
/// // Read the some lines from the GAF file.
/// let gaf_filename = utils::get_test_data("micb-kir3dl1_HG003.gaf");
/// let mut gaf_file = utils::open_file(&gaf_filename).unwrap();
/// let mut alignments = Vec::new();
/// for _ in 0..10 {
///     let mut buf: Vec<u8> = Vec::new();
///     let _ = gaf_file.read_until(b'\n', &mut buf).unwrap();
///     let mut aln = Alignment::from_gaf(&buf).unwrap();
///     // A block cannot have a mix of aligned and unaligned reads.
///     assert!(!aln.is_unaligned());
///     // We do not store unsupported optional fields.
///     aln.optional.clear();
///     alignments.push(aln);
/// }
///
/// // Compress the block.
/// let mut first_id = 0;
/// let block = AlignmentBlock::new(&alignments, &index, first_id).unwrap();
/// assert_eq!(block.len(), alignments.len());
/// first_id += alignments.len(); // Next block would start there.
///
/// // Decompress the block and extract the paths from the GBWT.
/// let mut decompressed = block.decode().unwrap();
/// assert_eq!(decompressed.len(), alignments.len());
/// for (i, aln) in decompressed.iter_mut().enumerate() {
///     aln.extract_target_path(&index);
///     assert_eq!(*aln, alignments[i]);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct AlignmentBlock {
    /// Minimum GBWT node identifier in the target paths, or [`None`] if this is a block of unaligned reads.
    pub min_handle: Option<usize>,
    /// Maximum GBWT node identifier in the target paths, or [`None`] if this is a block of unaligned reads.
    pub max_handle: Option<usize>,
    /// Number of alignments in the block.
    pub alignments: usize,
    /// Expected read length in the block, or [`None`] if the lengths vary.
    pub read_length: Option<usize>,
    /// GBWT starting positions for the target paths.
    pub gbwt_starts: Vec<u8>,
    /// Read and pair names.
    pub names: Vec<u8>,
    /// Quality strings.
    pub quality_strings: Vec<u8>,
    /// Difference strings.
    pub difference_strings: Vec<u8>,
    /// Binary flags for each alignment.
    pub flags: Flags,
    /// Encoded numerical information that cannot be derived from the other fields.
    pub numbers: Vec<u8>,
}

impl AlignmentBlock {
    // TODO: Make this a parameter?
    /// Compression level for Zstandard.
    pub const COMPRESSION_LEVEL: i32 = 7;

    // TODO: Allow outliers instead of requiring the same length from every read.
    fn expected_read_length(alignments: &[Alignment]) -> Option<usize> {
        let mut read_length = None;
        for aln in alignments.iter() {
            if let Some(len) = read_length {
                if aln.seq_len != len {
                    read_length = None;
                    break;
                }
            } else {
                if aln.seq_len == 0 {
                    break;
                }
                read_length = Some(aln.seq_len);
            }
        }
        read_length
    }

    fn compress_gbwt_starts(alignments: &[Alignment], index: &GBWT, first_id: usize, min_handle: Option<usize>) -> Result<Vec<u8>, String> {
        let base_node = min_handle.unwrap_or(0);
        let mut encoder = ByteCode::new();
        for i in 0..alignments.len() {
            // GBWT start as (node - min_handle, offset).
            let gbwt_sequence_id = if index.is_bidirectional() {
                support::encode_path(first_id + i, Orientation::Forward)
            } else { first_id + i };
            if let Some(start) = index.start(gbwt_sequence_id) {
                encoder.write(start.node - base_node);
                encoder.write(start.offset);
            } else if !encoder.is_empty() {
                return Err(format!("Line {}: Unaligned read in a block with aligned reads", first_id + i));
            }
        }
        Ok(Vec::from(encoder))
    }

    // TODO: Somewhere else?
    fn zstd_compress(data: &[u8]) -> Result<Vec<u8>, String> {
        let mut encoder = ZstdEncoder::new(Vec::new(), Self::COMPRESSION_LEVEL).unwrap();
        encoder.write_all(data).map_err(|err| format!("Zstd compression error: {}", err))?;
        encoder.finish().map_err(|err| format!("Zstd compression error: {}", err))
    }

    fn compress_names_pairs(alignments: &[Alignment]) -> Result<Vec<u8>, String> {
        let mut names: Vec<u8> = Vec::new();
        for aln in alignments.iter() {
            // Read name as a 0-terminated string.
            names.extend_from_slice(aln.name.as_bytes());
            names.push(0);
            // Pair name as a 0-terminated string, if present.
            if let Some(pair) = &aln.pair {
                names.extend_from_slice(pair.name.as_bytes());
            }
            names.push(0);
        }
        Self::zstd_compress(&names)
    }

    fn compress_quality_strings(alignments: &[Alignment]) -> Result<Vec<u8>, String> {
        let mut quality_strings: Vec<u8> = Vec::new();
        for aln in alignments.iter() {
            // Quality strings as 0-terminated strings.
            if !aln.base_quality.is_empty() {
                quality_strings.extend_from_slice(&aln.base_quality);
            }
            quality_strings.push(0);
        }
        Self::zstd_compress(&quality_strings)
    }

    /// Creates a new alignment block from the given read alignments and GBWT index.
    ///
    /// If the reads are aligned, they correspond to paths `first_id` to `first_id + alignments.len() - 1` in the GBWT index.
    /// The GBWT index may be bidirectional or unidirectional.
    ///
    /// # Arguments
    ///
    /// * `alignments`: The alignments to include in the block.
    /// * `index`: The GBWT index to use for encoding the GBWT starts.
    /// * `first_id`: Path identifier of the first alignment in the block.
    ///
    /// # Errors
    ///
    /// Returns an error if the block contains a mix of aligned and unaligned reads.
    /// Returns an error if compression fails.
    pub fn new(alignments: &[Alignment], index: &GBWT, first_id: usize) -> Result<Self, String> {
        let min_handle = alignments.iter().map(|aln| aln.min_handle()).min().flatten();
        let max_handle = alignments.iter().map(|aln| aln.max_handle()).max().flatten();
        let read_length = Self::expected_read_length(alignments);
        let gbwt_starts = Self::compress_gbwt_starts(alignments, index, first_id, min_handle)?;
        let names = Self::compress_names_pairs(alignments)?;
        let quality_strings = Self::compress_quality_strings(alignments)?;

        // We encode the remaining fields together, as they depend on each other.
        let mut difference_strings = RLE::with_sigma(Difference::NUM_TYPES);
        let mut flags = Flags::new(alignments.len());
        let mut numbers = ByteCode::new();

        for (i, aln) in alignments.iter().enumerate() {
            let perfect_alignment = aln.is_perfect() && Some(aln.seq_len) == read_length;
            if !perfect_alignment {
                aln.encode_difference_into(&mut difference_strings);
            }

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
            flags.set(i, Flags::FLAG_PERFECT_ALIGNMENT, perfect_alignment);

            aln.encode_numbers_into(&mut numbers, perfect_alignment);
        }

        Ok(Self {
            min_handle,
            max_handle,
            alignments: alignments.len(),
            read_length,
            gbwt_starts,
            names,
            quality_strings,
            difference_strings: Vec::from(difference_strings),
            flags,
            numbers: Vec::from(numbers),
        })
    }

    pub fn len(&self) -> usize {
        self.alignments
    }

    pub fn is_empty(&self) -> bool {
        self.alignments == 0
    }

    fn decompress_gbwt_starts(&self, result: &mut [Alignment]) -> Result<(), String> {
        if self.min_handle.is_none() {
            return Ok(());
        }
        let base_node = self.min_handle.unwrap();
        let mut decoder: ByteCodeIter<'_> = ByteCodeIter::new(&self.gbwt_starts[..]);
        for (i, aln) in result.iter_mut().enumerate() {
            let start = decoder.next().ok_or(
                format!("Missing GBWT start for alignment {}", i)
            )?;
            let offset = decoder.next().ok_or(
                format!("Missing GBWT offset for alignment {}", i)
            )?;
            aln.path = TargetPath::StartPosition(Pos::new(start + base_node, offset));
        }
        Ok(())
    }

    // TODO: Somewhere else?
    fn zstd_decompress(data: &[u8]) -> Result<Vec<u8>, String> {
        let mut decoder = ZstdDecoder::new(data).map_err(|err| format!("Zstd decompression error: {}", err))?;
        let mut buffer = Vec::new();
        decoder.read_to_end(&mut buffer).map_err(|err| format!("Zstd decompression error: {}", err))?;
        Ok(buffer)
    }

    fn decompress_names_pairs(&self, result: &mut [Alignment]) -> Result<(), String> {
        let buffer = Self::zstd_decompress(&self.names[..])?;
        let mut iter = buffer.split(|&c| c == 0);
        for (i, aln) in result.iter_mut().enumerate() {
            let name = iter.next().ok_or(format!("Missing name for alignment {}", i))?;
            aln.name = String::from_utf8_lossy(name).to_string();
            let pair_name = iter.next().ok_or(format!("Missing pair name for alignment {}", i))?;
            if !pair_name.is_empty() {
                aln.pair = Some(PairedRead {
                    name: String::from_utf8_lossy(pair_name).to_string(),
                    is_next: self.flags.get(i, Flags::FLAG_PAIR_IS_NEXT),
                    is_proper: self.flags.get(i, Flags::FLAG_PAIR_IS_PROPER),
                });
            }
        }
        // There should be an empty slice at the end, as the buffer is 0-terminated.
        // TODO: Should we check this?
        Ok(())
    }

    fn decompress_quality_strings(&self, result: &mut [Alignment]) -> Result<(), String> {
        let buffer = Self::zstd_decompress(&self.quality_strings[..])?;
        let mut iter = buffer.split(|&c| c == 0);
        for (i, aln) in result.iter_mut().enumerate() {
            let quality = iter.next().ok_or(format!("Missing quality string for alignment {}", i))?;
            aln.base_quality = quality.to_vec();
        }
        // There should be an empty slice at the end, as the buffer is 0-terminated.
        // TODO: Should we check this?
        Ok(())
    }

    fn decompress_difference_strings(&self, result: &mut [Alignment]) -> Result<(), String> {
        let mut decoder = RLEIter::with_sigma(&self.difference_strings[..], Difference::NUM_TYPES);
        for (i, aln) in result.iter_mut().enumerate() {
            if self.flags.get(i, Flags::FLAG_PERFECT_ALIGNMENT) {
                let len = self.read_length.ok_or(format!("Missing read length for perfect alignment {}", i))?;
                aln.difference = vec![Difference::Match(len)];
            } else {
                let res = Alignment::decode_difference_from(&self.difference_strings, &mut decoder);
                aln.difference = res.map_err(|err| format!("Missing difference string for alignment {}: {}", i, err))?;
            }
        }
        Ok(())
    }

    fn decompress_numbers(&self, result: &mut [Alignment]) -> Result<(), String> {
        let mut decoder = ByteCodeIter::new(&self.numbers[..]);
        for (i, aln) in result.iter_mut().enumerate() {
            // If this is a perfect alignment, we created a difference string even if the original alignment did not have one.
            let include_redundant = aln.difference.is_empty();
            if !self.flags.get(i, Flags::FLAG_PERFECT_ALIGNMENT) {
                (aln.seq_interval, aln.seq_len) = Alignment::decode_coordinates(&mut decoder, include_redundant)
                    .ok_or(format!("Missing query sequence coordinates for alignment {}", i))?;
            }
            (aln.path_interval, aln.path_len) = Alignment::decode_coordinates(&mut decoder, include_redundant)
                .ok_or(format!("Missing target path coordinates for alignment {}", i))?;
            if include_redundant {
                aln.matches = decoder.next().ok_or(format!("Missing number of matches for alignment {}", i))?;
                aln.edits = decoder.next().ok_or(format!("Missing number of edits for alignment {}", i))?;
            }
            if self.flags.get(i, Flags::FLAG_HAS_MAPQ) {
                aln.mapq = Some(decoder.next().ok_or(format!("Missing mapping quality for alignment {}", i))?);
            }
            if self.flags.get(i, Flags::FLAG_HAS_SCORE) {
                aln.score = Some(Alignment::decode_signed(&mut decoder).ok_or(
                    format!("Missing alignment score for alignment {}", i)
                )?);
            }

            // Derive redundant numbers from the difference string.
            // This also handles numbers derived from a perfect alignment.
            if !include_redundant {
                let (query_len, target_len, num_matches, num_edits) = Difference::stats(&aln.difference);
                aln.seq_interval.end = aln.seq_interval.start + query_len;
                aln.seq_len += query_len;
                aln.path_interval.end = aln.path_interval.start + target_len;
                aln.path_len += target_len;
                aln.matches = num_matches;
                aln.edits = num_edits;
            }
        }

        Ok(())
    }

    pub fn decode(&self) -> Result<Vec<Alignment>, String> {
        let mut result = vec![Alignment::default(); self.len()];
        self.decompress_gbwt_starts(&mut result)?;
        self.decompress_names_pairs(&mut result)?;
        self.decompress_quality_strings(&mut result)?;
        self.decompress_difference_strings(&mut result)?;
        // We decompress the numbers last, as their presence may depend on the other fields.
        self.decompress_numbers(&mut result)?;
        Ok(result)
    }
}

//-----------------------------------------------------------------------------
