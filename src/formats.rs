//! Support for writing various file formats.
//!
//! ### GFA
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
//! Weights represent the number of duplicate paths collapsed into a single line.
//! They are stored as tag `WT` of type `i`.
//! The CIGAR string is typically relative to the reference path in the same graph component.
//! It is stored as tag `CG` of type `Z`.
//!
//! ### JSON
//!
//! The JSON format is a text-based format for representing structured data.
//! The support for it is based on building a [`JSONValue`] object recursively and then writing it using the [`Display`] trait.
//! There is also a helper function [`json_path`] for building a JSON object for a path with metadata.

// TODO: tests

use crate::{GBZPath, FullPathName, GBZRecord};

use std::fmt::Display;
use std::io::{self, Write};
use std::ops::Range;

use gbwt::{Metadata, Orientation};
use gbwt::support;

//-----------------------------------------------------------------------------

/// Metadata for a walk line in a GFA file.
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
    ///
    /// # Arguments
    ///
    /// * `path`: The path.
    /// * `interval`: The interval of the path.
    /// * `weight`: Optional weight for the path.
    pub fn path_interval(path: &GBZPath, interval: Range<usize>, weight: Option<usize>) -> Self {
        let mut name = path.name.clone();
        let end = name.fragment + interval.end;
        name.fragment += interval.start;
        WalkMetadata { name, end, weight, cigar: None }
    }

    /// Creates new metadata for a haplotype using GBWT metadata.
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
    /// * `weight`: Optional weight for the path.
    pub fn anonymous(haplotype: usize, contig: &str, len: usize, weight: Option<usize>) -> Self {
        let path_name = FullPathName::haplotype("unknown", contig, haplotype, 0);
        WalkMetadata {
            name: path_name,
            end: len,
            weight,
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
pub fn write_gfa_header<T: Write>(reference_samples: Option<String>, output: &mut T) -> io::Result<()> {
    let header = if let Some(sample_names) = reference_samples {
        format!("H\tVN:Z:1.1\tRS:Z:{}\n", sample_names)
    } else {
        "H\tVN:Z:1.1\n".to_string()
    };
    output.write_all(header.as_bytes())?;
    Ok(())
}

/// Writes a GFA segment line corresponding to a node.
pub fn write_gfa_segment<T: Write>(record: &GBZRecord, output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();

    let (id, orientation) = support::decode_node(record.handle());
    buffer.extend_from_slice(b"S\t");
    buffer.extend_from_slice(id.to_string().as_bytes());
    buffer.push(b'\t');
    if orientation == Orientation::Reverse {
        buffer.extend_from_slice(record.sequence());
    } else {
        let rc = support::reverse_complement(record.sequence());
        buffer.extend_from_slice(&rc);
    }
    buffer.push(b'\n');

    output.write_all(&buffer)?;
    Ok(())
}

/// Writes a GFA link line corresponding to an edge.
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
/// * `weight`: Number of duplicate paths collapsed into a single line.
/// * `cigar`: An optional CIGAR string for the path.
/// * `path`: An array of objects with fields `id` (string) and `is_reverse` (boolean) for each node visit in the path.
pub fn json_path(path: &[usize], metadata: &WalkMetadata) -> JSONValue {
    let mut values: Vec<(String, JSONValue)> = Vec::new();
    values.push(("name".to_string(), JSONValue::String(metadata.name.path_fragment_name(metadata.end))));
    values.push(("weight".to_string(), JSONValue::Number(metadata.weight.unwrap_or(1))));
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
