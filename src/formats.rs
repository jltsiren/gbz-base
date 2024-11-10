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

use std::fmt::Display;
use std::io::{self, Write};
use std::ops::Range;
use std::str;

use gbwt::{Metadata, Orientation, FullPathName};
use gbwt::support;

#[cfg(test)]
mod tests;

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

// FIXME implement, document, examples, test
// TODO: Parse pairing information from tags.
#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    /// Name or identifier of the query sequence.
    pub name: SequenceName,
    /// Length of the query sequence.
    pub seq_len: usize,
    /// Aligned interval of the query sequence.
    pub seq_interval: Range<usize>,
    /// Target path in the orientation of the query sequence.
    pub path: Vec<(usize, Orientation)>,
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
    pub optional: Option<Vec<TypedField>>,
}

// FIXME implement, document, examples, test
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceName {
    /// Sequence name.
    Name(String),
    /// Offset of the name in a list of sequence names.
    Identifier(usize),
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
/// assert!(ops.is_some());
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
    /// Returns [`None`] if the difference string is invalid.
    pub fn parse(difference_string: &[u8]) -> Option<Vec<Self>> {
        let mut result: Vec<Self> = Vec::new();
        if difference_string.is_empty() {
            return Some(result);
        }
        if !Self::OPS.contains(&difference_string[0]) {
            return None;
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
                _ => None,
            }?;
            result.push(op);
            start = end;
        }

        Some(result)
    }

    /// Parses a difference string and returns it as a normalized vector of operations.
    ///
    /// The operations are merged and empty operations are removed.
    /// Returns [`None`] if the difference string is invalid.
    pub fn parse_normalized(difference_string: &[u8]) -> Option<Vec<Self>> {
        let ops = Self::parse(difference_string)?;
        Some(Self::normalize(ops))
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
/// assert_eq!(field, Some(TypedField::Int([b'A', b'S'], 160)));
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
    pub fn parse(field: &[u8]) -> Option<Self> {
        if field.len() < 5 || field[2] != b':' || field[4] != b':' {
            return None;
        }
        let tag = [field[0], field[1]];
        match field[3] {
            b'A' => {
                if field.len() != 6 {
                    return None;
                }
                Some(TypedField::Char(tag, field[5]))
            },
            b'Z' => Some(TypedField::String(tag, field[5..].to_vec())),
            b'i' => {
                let value = str::from_utf8(&field[5..]).ok()?;
                let value = value.parse::<isize>().ok()?;
                Some(TypedField::Int(tag, value))
            },
            b'f' => {
                let value = str::from_utf8(&field[5..]).ok()?;
                let value = value.parse::<f64>().ok()?;
                Some(TypedField::Float(tag, value))
            },
            b'b' => {
                if field.len() != 6 {
                    return None;
                }
                match field[5] {
                    b'0' => Some(TypedField::Bool(tag, false)),
                    b'1' => Some(TypedField::Bool(tag, true)),
                    _ => None,
                }
            },
            _ => None,
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
