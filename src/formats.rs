//! Support for reading and writing various file formats.
//!
//! ### GFA (writing)
//!
//! The GFA format is a text-based format for representing sequence graphs.
//! See [the specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for details.
//! The following functions support line-by-line writing of GFA version 1.1:
//!
//! * [`write_gfa_header`]: Write the header line.
//! * [`write_gfa_segment`]: Write a segment line for a node.
//! * [`write_gfa_link`]: Write a link line for an edge.
//! * [`write_gfa_walk`]: Write a walk line for a path.
//!
//! A walk line contains metadata, which is stored in a [`WalkMetadata`] object.
//! The object contains a structured path name, the end position of the path, an optional weight, and an optional CIGAR string.
//!
//! ### JSON (writing)
//!
//! The JSON format is a text-based format for representing structured data.
//! The support for it is based on building a [`JSONValue`] object recursively and then writing it using the [`Display`] trait.
//! There is also a helper function [`json_path`] for building a JSON object for a path with metadata.
//!
//! ### GAF (reading)
//!
//! The GAF format is a text-based format for representing sequence alignments to a graph.
//! See [the specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) for an overview.
//! Some details are better documented in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html#10).
//!
//! FIXME: details

use crate::db;

use std::fmt::Display;
use std::io::{self, Write};
use std::ops::Range;
use std::path::PathBuf;
use std::str;

use gbwt::{GBWT, Metadata, Orientation, Pos, FullPathName};
use gbwt::support::{self, ByteCode, ByteCodeIter, RLE, Run, RLEIter};

use simple_sds::ops::Select;
use simple_sds::raw_vector::{RawVector, AccessRaw};
use simple_sds::sparse_vector::{SparseVector, SparseBuilder};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

// TODO: Move somewhere else?
/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

//-----------------------------------------------------------------------------

/// Metadata for a walk line in a GFA file.
///
/// In addition to the standard fields, the metadata may contain an optional weight and an optional CIGAR string.
/// Weights represent the number of duplicate paths collapsed into a single line.
/// They are stored as tag `WT` of type `i`.
/// The CIGAR string is typically relative to the reference path in the same graph component.
/// It is stored as tag `CG` of type `Z`.
pub struct WalkMetadata {
    // Structured name with a sample name, contig name, haplotype / phase number, and starting offset.
    name: FullPathName,

    // Past-the-end offset of the path.
    end: usize,

    // Optional weight for the path, representing the number of identical paths collapsed into a single line.
    weight: Option<usize>,

    // Optional CIGAR string for the path, typically relative to the reference path in the same graph component.
    cigar: Option<String>,
}

impl WalkMetadata {
    /// Creates new metadata for an interval of a path.
    pub fn path_interval(path_name: &FullPathName, interval: Range<usize>) -> Self {
        let mut name = path_name.clone();
        let end = name.fragment + interval.end;
        name.fragment += interval.start;
        WalkMetadata { name, end, weight: None, cigar: None }
    }

    /// Creates new metadata for a haplotype path using GBWT metadata.
    ///
    /// Returns [`None`] if the path does not exist.
    ///
    /// # Arguments
    ///
    /// * `metadata`: The GBWT metadata.
    /// * `path_id`: The path identifier.
    /// * `len`: The length of the path in base pairs.
    pub fn haplotype(metadata: &Metadata, path_id: usize, len: usize) -> Option<Self> {
        let name = FullPathName::from_metadata(metadata, path_id)?;
        Some(WalkMetadata { name, end: len, weight: None, cigar: None })
    }

    /// Creates new metadata for a haplotype of unknown origin.
    ///
    /// # Arguments
    ///
    /// * `haplotype`: The haplotype / phase number.
    /// * `contig`: The contig name.
    /// * `len`: The length of the path in base pairs.
    pub fn anonymous(haplotype: usize, contig: &str, len: usize) -> Self {
        let path_name = FullPathName::haplotype("unknown", contig, haplotype, 0);
        WalkMetadata {
            name: path_name,
            end: len,
            weight: None,
            cigar: None,
        }
    }

    /// Adds a weight to the metadata.
    pub fn add_weight(&mut self, weight: Option<usize>) {
        self.weight = weight;
    }

    /// Adds a CIGAR string to the metadata.
    pub fn add_cigar(&mut self, cigar: Option<String>) {
        self.cigar = cigar;
    }
}

//-----------------------------------------------------------------------------

// TODO: These should be shared with gbunzip.

/// Writes the GFA header line.
///
/// The header line may contain a list of reference sample names.
/// Following the convention set by vg, the reference sample names are stored as a string in the `RS` tag of type `Z`.
/// If there are multiple reference samples, their names are separated by a single space.
pub fn write_gfa_header<T: Write>(reference_samples: Option<&str>, output: &mut T) -> io::Result<()> {
    let header = if let Some(sample_names) = reference_samples {
        format!("H\tVN:Z:1.1\tRS:Z:{}\n", sample_names)
    } else {
        "H\tVN:Z:1.1\n".to_string()
    };
    output.write_all(header.as_bytes())?;
    Ok(())
}

/// Writes a GFA segment line corresponding to a node with an integer identifier.
pub fn write_gfa_node<T: Write>(node_id: usize, sequence: &[u8], output: &mut T) -> io::Result<()> {
    write_gfa_segment(node_id.to_string().as_bytes(), sequence, output)
}

/// Writes a GFA segment line corresponding to a segment with a string name.
pub fn write_gfa_segment<T: Write>(name: &[u8], sequence: &[u8], output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();

    buffer.extend_from_slice(b"S\t");
    buffer.extend_from_slice(name);
    buffer.push(b'\t');
    buffer.extend_from_slice(sequence);
    buffer.push(b'\n');

    output.write_all(&buffer)?;
    Ok(())
}

/// Writes a GFA link line corresponding to an edge between two oriented nodes.
pub fn write_gfa_edge<T: Write>(from: (usize, Orientation), to: (usize, Orientation), output: &mut T) -> io::Result<()> {
    write_gfa_link(
        (from.0.to_string().as_bytes(), from.1),
        (to.0.to_string().as_bytes(), to.1),
        output
    )
}

/// Writes a GFA link line corresponding to a link between two oriented segments.
pub fn write_gfa_link<T: Write>(from: (&[u8], Orientation), to: (&[u8], Orientation), output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();

    buffer.extend_from_slice(b"L\t");
    buffer.extend_from_slice(from.0);
    match from.1 {
        Orientation::Forward => buffer.extend_from_slice(b"\t+\t"),
        Orientation::Reverse => buffer.extend_from_slice(b"\t-\t"),
    }
    buffer.extend_from_slice(to.0);
    match to.1 {
        Orientation::Forward => buffer.extend_from_slice(b"\t+\t0M\n"),
        Orientation::Reverse => buffer.extend_from_slice(b"\t-\t0M\n"),
    }

    output.write_all(&buffer)?;
    Ok(())
}

/// Writes a GFA walk line corresponding to a path.
pub fn write_gfa_walk<T: Write>(path: &[usize], metadata: &WalkMetadata, output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();

    buffer.extend_from_slice(b"W\t");
    buffer.extend_from_slice(metadata.name.sample.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.name.haplotype.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.name.contig.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.name.fragment.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.end.to_string().as_bytes());
    buffer.push(b'\t');
    for handle in path.iter() {
        match support::node_orientation(*handle) {
            Orientation::Forward => buffer.push(b'>'),
            Orientation::Reverse => buffer.push(b'<'),
        }
        buffer.extend_from_slice(support::node_id(*handle).to_string().as_bytes());
    }
    if let Some(weight) = metadata.weight {
        buffer.extend_from_slice(b"\tWT:i:");
        buffer.extend_from_slice(weight.to_string().as_bytes());
    }
    if let Some(cigar) = &metadata.cigar {
        buffer.extend_from_slice(b"\tCG:Z:");
        buffer.extend_from_slice(cigar.as_bytes());
    }
    buffer.push(b'\n');

    output.write_all(&buffer)?;
    Ok(())
}

//-----------------------------------------------------------------------------

/// A structured JSON value.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum JSONValue {
    /// A boolean value.
    Boolean(bool),

    /// A string value.
    String(String),

    /// A number value.
    Number(usize),

    /// An array of JSON values.
    Array(Vec<JSONValue>),

    /// A JSON object storing a list of JSON values with string names.
    Object(Vec<(String, JSONValue)>),
}

impl Display for JSONValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            JSONValue::Boolean(b) => write!(f, "{}", b),
            JSONValue::String(s) => write!(f, "\"{}\"", s),
            JSONValue::Number(n) => write!(f, "{}", n),
            JSONValue::Array(v) => {
                write!(f, "[")?;
                let mut first = true;
                for value in v.iter() {
                    if first {
                        first = false;
                    } else {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", value)?;
                }
                write!(f, "]")
            },
            JSONValue::Object(v) => {
                write!(f, "{{")?;
                let mut first = true;
                for (key, value) in v.iter() {
                    if first {
                        first = false;
                    } else {
                        write!(f, ", ")?;
                    }
                    write!(f, "\"{}\": {}", key, value)?;
                }
                write!(f, "}}")
            },
        }
    }
}

/// Creates a JSON object for a path with metadata.
///
/// The object contains the following fields:
///
/// * `name`: A path name compatible with vg (see [`FullPathName::path_fragment_name`]).
/// * `weight`: An optional number of duplicate paths collapsed into a single line.
/// * `cigar`: An optional CIGAR string for the path.
/// * `path`: An array of objects with fields `id` (string) and `is_reverse` (boolean) for each node visit in the path.
pub fn json_path(path: &[usize], metadata: &WalkMetadata) -> JSONValue {
    let mut values: Vec<(String, JSONValue)> = Vec::new();
    values.push(("name".to_string(), JSONValue::String(metadata.name.path_fragment_name(metadata.end))));
    if let Some(weight) = metadata.weight {
        values.push(("weight".to_string(), JSONValue::Number(weight)));
    }
    if let Some(cigar) = &metadata.cigar {
        values.push(("cigar".to_string(), JSONValue::String(cigar.clone())));
    }
    values.push(("path".to_string(), JSONValue::Array(path.iter().map(
        |x| JSONValue::Object(vec![
            ("id".to_string(), JSONValue::String(support::node_id(*x).to_string())),
            ("is_reverse".to_string(), JSONValue::Boolean(support::node_orientation(*x) == Orientation::Reverse)),
        ])
    ).collect())));

    JSONValue::Object(values)
}

//-----------------------------------------------------------------------------

// TODO: Parse pairing information from tags.
// FIXME examples
/// An alignment between a query sequence and a target path in a graph.
///
/// This object corresponds either to a line in a GAF file or to a row in table `Alignments` in [`crate::GAFBase`].
/// When the alignment is built from a GAF line, the name of the query sequence and the target path are stored explicitly.
/// For alignments stored in a database, these are stored relative to other tables.
/// See [`SequenceName`] and [`TargetPath`] for details and [`Self::set_relative_information`] for conversion.
#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    /// Name or identifier of the query sequence.
    pub name: SequenceName,
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
    /// Difference string.
    pub difference: Option<Vec<Difference>>,
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
        if field == &Self::MISSING_VALUE {
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
        if field == &Self::MISSING_VALUE {
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
        if field == &Self::MISSING_VALUE {
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
                format!("Only numerical segment names are supported")
            })?;
            result.push(support::encode_node(node, orientation));
            start = end;
        }
        Ok(result)
    }

    // TODO: Permissive mode where we ignore errors in optional fields? Strict mode that returns an error?
    /// Parses an alignment from a GAF line.
    ///
    /// Returns an error if the line cannot be parsed.
    /// The line may end with up to one endline character, which is ignored.
    /// The returned alignment stores the query sequence name and the target path explicitly.
    /// Parsing is based on bytes rather than characters to avoid unnecessary UTF-8 validation.
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
        let name = SequenceName::Name(Self::parse_string(fields[0], "query sequence name")?);
        let seq_len = Self::parse_usize(fields[1], "query sequence length")?;
        let seq_interval = Self::parse_interval(fields[2], fields[3])?;

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
        let matches = Self::parse_usize(fields[9], "matches")?;
        let alignment_len = Self::parse_usize(fields[10], "alignment length")?;
        let edits = if matches <= alignment_len { alignment_len - matches } else { 0 }; // TODO: Error in strict mode?
        let mapq = Self::parse_usize(fields[11], "mapping quality")?;
        let mapq = if mapq == Self::MISSING_MAPQ { None } else { Some(mapq) }; // TODO: Too large values?

        // Optional fields.
        let mut score = None;
        let mut base_quality = None;
        let mut difference = None;
        let mut target_len: Option<usize> = None; // TODO: Part of the interval end hack.
        let mut optional = Vec::new();
        for field in fields[Self::MANDATORY_FIELDS..].iter() {
            let parsed = TypedField::parse(field)?;
            match parsed {
                TypedField::Int([b'A', b'S'], value) => {
                    if let Some(_) = score {
                        return Err(String::from("Multiple alignment score fields"));
                    }
                    score = Some(value);
                },
                TypedField::String([b'b', b'q'], value) => {
                    if let Some(_) = base_quality {
                        return Err(String::from("Multiple base quality fields"));
                    }
                    base_quality = Some(value.to_vec());
                },
                TypedField::String([b'c', b's'], value) => {
                    if let Some(_) = difference {
                        return Err(String::from("Multiple difference fields"));
                    }
                    let ops = Difference::parse_normalized(&value)?;
                    target_len = Some(ops.iter().map(|op| op.target_len()).sum()); // TODO: Part of the interval end hack.
                    difference = Some(ops);
                },
                _ => { optional.push(parsed); },
            }
        }

        // TODO: Part of the interval end hack.
        let mut path_interval = match path_interval {
            Ok(interval) => interval,
            Err(_) => {
                if let Some(target_len) = target_len {
                    let target_start = Self::parse_usize(fields[7], "target interval start")?;
                    target_start..target_start + target_len
                } else {
                    return Err(String::from("Target interval end cannot be parsed or inferred from the difference string"));
                }
            },
        };

        // TODO: Move to the target path section once the interval end hack becomes unnecessary.
        // TODO: Error in strict mode, sanitize in permissive mode?
        if orientation == Orientation::Reverse {
            let start = if path_interval.end < path_len { path_len - path_interval.end } else { 0 };
            let end = if path_interval.start < path_len { path_len - path_interval.start } else { 0 };
            path_interval = start..end;
        }

        Ok(Alignment {
            name, seq_len, seq_interval,
            path, path_len, path_interval,
            matches, edits, mapq, score,
            base_quality, difference, optional,
        })
    }

    /// Returns the list of all path identifiers in GBWT metadata, sorted by sample identifier.
    ///
    /// Path identifiers are returned as [`u32`], because a GBWT index cannot have more than [`u32::MAX`] paths.
    /// Also returns a multiset [`SparseVector`] that maps sample identifiers to intervals in the returned list.
    /// The interval for sample identifier `i` is `select(i)..select(i + 1)`.
    /// Note that there is a sentinel to make the general formula work for the last sample.
    ///
    /// Returns an error if something goes wrong.
    /// [`crate::GAFBase::check_gbwt_metadata`] has more detailed error messages for insufficient metadata.
    ///
    /// # Panics
    ///
    /// May panic during [`SparseVector`] construction if the metadata is in an inconsistent state.
    pub fn paths_by_sample(metadata: &Metadata) -> Result<(Vec<u32>, SparseVector), String> {
        let mut sample_path: Vec<(u32, u32)> = metadata.path_iter().enumerate().map(|(id, name)| {
            (name.sample() as u32, id as u32)
        }).collect();
        sample_path.sort_unstable();

        let mut builder = SparseBuilder::multiset(metadata.paths() + 1, metadata.samples() + 1);
        let mut sample_id = 0;
        for (offset, (sample, _)) in sample_path.iter().enumerate() {
            while sample_id <= *sample {
                builder.set(offset);
                sample_id += 1;
            }
        }
        while (sample_id as usize) <= metadata.samples() {
            builder.set(metadata.paths());
            sample_id += 1;
        }
        let path_ids: Vec<u32> = sample_path.iter().map(|(_, path)| *path).collect();
        let path_index = SparseVector::try_from(builder).map_err(|err| {
            format!("Failed to build a sample-to-path index: {}", err)
        })?;

        Ok((path_ids, path_index))
    }

    // Replaces query sequence name with the corresponding sample identifier in the GBWT metadata.
    fn set_sequence_id(&mut self, metadata: &Metadata) -> Result<(), String> {
        let name = match &self.name {
            SequenceName::Name(n) => n,
            SequenceName::Identifier(_) => return Ok(()),
        };
        let id = metadata.sample_id(name).ok_or(
            format!("Sequence name {} not found in the GBWT metadata", name)
        )?;
        self.name = SequenceName::Identifier(id);
        Ok(())
    }

    // TODO: We could simplify this greatly by storing GAF line number as fragment id in the path name.
    // Replaces target path with its starting position in the given GBWT index.
    fn set_target_start(&mut self,
        index: &GBWT,
        paths_by_sample: &(Vec<u32>, SparseVector),
        used_paths: Option<&mut RawVector>
    ) -> Result<(), String> {
        if let TargetPath::StartPosition(_) = self.path {
            return Ok(());
        }
        let (path_ids, path_index) = paths_by_sample;

        // Determine the interval of possible path identifiers.
        let sample_id = match self.name {
            SequenceName::Name(_) => return Err(String::from("Sequence name not replaced with an identifier")),
            SequenceName::Identifier(id) => id,
        };
        let mut iter = path_index.select_iter(sample_id);
        let start = iter.next().unwrap_or((0, path_ids.len())).1;
        let end = iter.next().unwrap_or((0, path_ids.len())).1;

        // Check the paths in the interval.
        for offset in start..end {
            // Is it still available?
            let path_id = path_ids[offset] as usize;
            if let Some(used) = used_paths.as_ref() {
                if used.bit(path_id) {
                    continue;
                }
            }

            // Does it match the query sequence?
            let iter = index.sequence(support::encode_path(path_id, Orientation::Forward));
            if iter.is_none() {
                continue;
            }
            let gbwt_path = iter.unwrap();
            let query_path = match self.path {
                TargetPath::Path(ref path) => path.iter().copied(),
                _ => unreachable!(),
            };
            if !gbwt_path.eq(query_path) {
                continue;
            }

            // Mark the path as used and set the start position.
            if let Some(used) = used_paths {
                used.set_bit(path_id, true);
            }
            let start_pos = db::path_start(index, path_id, Orientation::Forward);
            self.path = TargetPath::StartPosition(start_pos);
            return Ok(());
        }

        Err(format!("Could not match query sequence {} to an unused GBWT path (interval: {}..{})", self.name, start, end))
    }

    /// Replaces the explicitly stored data with data relative to the given GBWT index.
    ///
    /// Replaces query sequence name with an identifier and the target path with a starting position.
    /// No effect if the data is already relative to a GBWT index.
    /// Returns an error if the replacement fails for any reason.
    /// [`crate::GAFBase::check_gbwt_metadata`] has more detailed error messages for insufficient metadata.
    ///
    /// # Arguments
    ///
    /// * `index`: The GBWT index.
    /// * `paths_by_sample`: Path identifiers indexed by sample identifier (from [`Self::paths_by_sample`]).
    /// * `used_paths`: Optional vector to mark paths that have been used and should not be used again.
    ///
    /// The length of `used_paths` must be at least `index.metadata().paths()`.
    ///
    /// # Panics
    ///
    /// May panic if `paths_by_sample` or `used_paths` are inconsistent with the GBWT metadata.
    pub fn set_relative_information(&mut self,
        index: &GBWT,
        paths_by_sample: &(Vec<u32>, SparseVector),
        used_paths: Option<&mut RawVector>
    ) -> Result<(), String> {
        let metadata = index.metadata().ok_or(
            String::from("The GBWT index does not contain metadata")
        )?;
        self.set_sequence_id(metadata)?;
        self.set_target_start(index, paths_by_sample, used_paths)
    }

    // TODO: back to a GAF line
}

//-----------------------------------------------------------------------------

/// Encoding / decoding the alignment in the format used in GAF-base.
impl Alignment {
    // Normalizes the interval and encodes it as (left flank, length, right flank).
    fn encode_coordinates(interval: Range<usize>, len: usize, encoder: &mut ByteCode) {
        let start = if interval.start <= interval.end { interval.start } else { interval.end };
        let end = if interval.end <= len { interval.end } else { len };
        encoder.write(start);
        encoder.write(end - start);
        encoder.write(len - end);
    }

    // Encodes a signed integer. Small absolute values are represented as small numbers.
    fn encode_signed(value: isize, encoder: &mut ByteCode) {
        let value = if value < 0 { (-2 * value - 1) as usize } else { 2 * value as usize };
        encoder.write(value);
    }

    // TODO: Include fw_offset?
    /// Encodes numerical fields as a blob.
    ///
    /// This includes the following fields:
    ///
    /// * Query sequence: `seq_len`, `seq_interval.start`, `seq_interval.end`.
    /// * Target path: `path_len`, `path_interval.start`, `path_interval.end`.
    /// * Alignment statistics: `matches`, `edits`, `mapq`, `score`.
    ///
    /// The coordinates are normalized so that `start <= end <= len`.
    pub fn encode_numbers(&self) -> Vec<u8> {
        let mut encoder = ByteCode::new();

        // Coordinates.
        Self::encode_coordinates(self.seq_interval.clone(), self.seq_len, &mut encoder);
        Self::encode_coordinates(self.path_interval.clone(), self.path_len, &mut encoder);

        // Alignment statistics.
        encoder.write(self.matches);
        encoder.write(self.edits);
        encoder.write(self.mapq.unwrap_or(Self::MISSING_MAPQ));
        if let Some(score) = self.score {
            Self::encode_signed(score, &mut encoder);
        };

        Vec::from(encoder)
    }

    /// Encodes the base quality sequence as a run-length encoded blob.
    ///
    /// Returns [`None`] if the base quality sequence is missing.
    pub fn encode_base_quality(&self) -> Option<Vec<u8>> {
        if self.base_quality.is_none() {
            return None;
        }

        // Determine the alphabet and encode it as (alphabet size, alphabet).
        let mut alphabet: Vec<u8> = Vec::new();
        for value in self.base_quality.as_ref().unwrap().iter() {
            if !alphabet.contains(value) {
                alphabet.push(*value);
            }
        }
        let mut encoder = RLE::with_sigma(alphabet.len());
        encoder.write_int(alphabet.len());
        for &value in alphabet.iter() {
            encoder.write_byte(value);
        }

        // Run-length encode the base quality sequence.
        let mut run = Run::new(0, 0);
        for &value in self.base_quality.as_ref().unwrap().iter() {
            if value as usize == run.value {
                run.len += 1;
            } else {
                if run.len > 0 {
                    run.value = alphabet.iter().position(|&x| x as usize == run.value).unwrap();
                    encoder.write(run);
                }
                run.value = value as usize;
                run.len = 1;
            }
        }
        if run.len > 0 {
            run.value = alphabet.iter().position(|&x| x as usize == run.value).unwrap();
            encoder.write(run);
        }

        Some(Vec::from(encoder))
    }

    // TODO: Use a more efficient encoding for mismatches and insertions.
    /// Encodes the difference string as a blob.
    ///
    /// Returns [`None`] if the difference string is missing.
    pub fn encode_difference(&self) -> Option<Vec<u8>> {
        if self.difference.is_none() {
            return None;
        }

        let mut encoder = RLE::with_sigma(Difference::NUM_TYPES);
        for diff in self.difference.as_ref().unwrap().iter() {
            match diff {
                Difference::Match(len) => encoder.write(Run::new(0, *len)),
                Difference::Mismatch(base) => encoder.write(Run::new(1, *base as usize)),
                Difference::Insertion(seq) => {
                    encoder.write(Run::new(2, seq.len()));
                    for byte in seq.iter() {
                        encoder.write_byte(*byte);
                    }
                },
                Difference::Deletion(len) => encoder.write(Run::new(3, *len)),
            }
        }

        Some(Vec::from(encoder))
    }

    // TODO: encode optional fields

    // Decodes an interval and a length from the iterator.
    fn decode_coordinates(iter: &mut ByteCodeIter) -> Option<(Range<usize>, usize)> {
        let left_flank = iter.next()?;
        let interval_len = iter.next()?;
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

    /// Decodes the alignment from the given fields.
    ///
    /// Returns an error if the fields cannot be decoded.
    ///
    /// # Arguments
    ///
    /// * `query`: Query sequence identifier.
    /// * `target`: Target path starting position.
    /// * `numbers`: Encoded numerical fields.
    /// * `quality`: Encoded base quality sequence.
    /// * `difference`: Encoded difference string.
    pub fn decode(query: usize, target: Pos, numbers: &[u8], quality: Option<&[u8]>, difference: Option<&[u8]>) -> Result<Self, String> {
        let name = SequenceName::Identifier(query);
        let path = TargetPath::StartPosition(target);

        // Decode the numerical fields.
        let mut number_decoder = ByteCodeIter::new(numbers);
        let (seq_interval, seq_len) = Self::decode_coordinates(&mut number_decoder).ok_or(
            String::from("Missing query sequence coordinates")
        )?;
        let (path_interval, path_len) = Self::decode_coordinates(&mut number_decoder).ok_or(
            String::from("Missing target path coordinates")
        )?;
        let matches = number_decoder.next().ok_or(String::from("Missing number of matches"))?;
        let edits = number_decoder.next().ok_or(String::from("Missing number of edits"))?;
        let mapq = match number_decoder.next() {
            Some(mapq) => if mapq == Self::MISSING_MAPQ { None } else { Some(mapq) },
            None => return Err(String::from("Missing mapping quality")),
        };
        let score = Self::decode_signed(&mut number_decoder);

        // Decode the base quality sequence.
        let base_quality = if let Some(quality) = quality {
            let mut quality_decoder = RLEIter::new(quality);
            let sigma = quality_decoder.int().ok_or(String::from("Missing base quality alphabet size"))?;
            quality_decoder.set_sigma(sigma);
            let mut alphabet: Vec<u8> = Vec::with_capacity(sigma);
            for i in 0..sigma {
                let value = quality_decoder.byte().ok_or(format!("Missing base quality symbol {}", i))?;
                alphabet.push(value);
            }
            let mut result: Vec<u8> = Vec::with_capacity(seq_len);
            while let Some(run) = quality_decoder.next() {
                let value = alphabet[run.value];
                for _ in 0..run.len {
                    result.push(value);
                }
            }
            Some(result)
        } else {
            None
        };

        // Decode the difference string.
        let difference = if let Some(difference) = difference {
            let mut difference_decoder = RLEIter::with_sigma(difference, Difference::NUM_TYPES);
            let mut result: Vec<Difference> = Vec::new();
            while let Some(run) = difference_decoder.next() {
                match run.value {
                    0 => result.push(Difference::Match(run.len)),
                    1 => result.push(Difference::Mismatch(run.len as u8)),
                    2 => {
                        let mut seq: Vec<u8> = Vec::with_capacity(run.len);
                        for _ in 0..run.len {
                            seq.push(difference_decoder.byte().ok_or(String::from("Missing insertion base"))?);
                        }
                        result.push(Difference::Insertion(seq));
                    },
                    3 => result.push(Difference::Deletion(run.len)),
                    _ => return Err(format!("Invalid difference string operation: {}", run.value)),
                }
            }
            Some(result)
        } else {
            None
        };

        // TODO: Decode optional fields
        let optional = Vec::new();

        Ok(Alignment {
            name, seq_len, seq_interval,
            path, path_len, path_interval,
            matches, edits, mapq, score,
            base_quality, difference, optional,
        })
    }
}

//-----------------------------------------------------------------------------

/// Name of a query sequence in [`Alignment`].
///
/// The name is stored either explicitly as a string or as an identifier relative to a list of sequence names.
/// When the name is for an alignment stored in [`crate::GAFBase`], the identifier is the handle in table `Sequences`.
/// During the construction of the database, this corresponds to a sample identifier in GBWT metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceName {
    /// Sequence name.
    Name(String),
    /// Offset of the name in a list of sequence names.
    Identifier(usize),
}

impl Display for SequenceName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SequenceName::Name(name) => write!(f, "(name: {})", name),
            SequenceName::Identifier(id) => write!(f, "(id: {})", id),
        }
    }
}

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
/// use gbz_base::formats::Difference;
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
}

impl Difference {
    /// Number of supported operation types.
    pub const NUM_TYPES: usize = 4;

    // TODO: This does not support `~` (intron length and splice signal) yet.
    const OPS: &'static [u8] = b"=:*+-";

    fn base_to_upper(c: u8) -> u8 {
        if c >= b'a' && c <= b'z' {
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

    /// Returns the length of the operation.
    pub fn len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(len) => *len,
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
        }
    }

    /// Returns the length of the operation in the query sequence.
    pub fn query_len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(_) => 0,
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
/// use gbz_base::formats::TypedField;
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
