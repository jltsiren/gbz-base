//! Support for reading and writing various file formats.
//!
//! ### GFA (writing)
//!
//! The GFA format is a text-based format for representing sequence graphs.
//! See [the specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for details.
//! The following functions support line-by-line writing of GFA version 1.1:
//!
//! * [`write_gfa_header`]: Write a GFA file header.
//! * [`write_header_lines`]: Write additional header lines.
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
//! ### GAF (reading and writing)
//!
//! The GAF format is a text-based format for representing sequence alignments to a graph.
//! See [the specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) for an overview.
//! Some details are better documented in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html#10).
//!
//! GAF header lines can be handled with the following functions:
//!
//! * [`is_gaf_header_line`]: Check if a buffer contains a GAF header line.
//! * [`peek_gaf_header_line`]: Check if the next line in a reader would be a GAF header line.
//! * [`read_gaf_header_lines`]: Read all successive GAF header lines from a reader.
//! * [`write_gaf_file_header`]: Write a GAF file header.
//! * [`write_header_lines`]: Write additional header lines.
//!
//! I/O for GAF alignment lines is currently implemented in [`crate::Alignment`].

use crate::utils;

use std::fmt::Display;
use std::io::{self, BufRead, Write};
use std::ops::Range;
use std::str;

use gbwt::{Metadata, Orientation, FullPathName};
use gbwt::support;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Appends a GFA walk / GAF path to a string represented as `Vec<u8>`.
pub fn append_walk(buffer: &mut Vec<u8>, walk: &[usize]) {
    for handle in walk {
        match support::node_orientation(*handle) {
            Orientation::Forward => buffer.push(b'>'),
            Orientation::Reverse => buffer.push(b'<'),
        }
        utils::append_usize(buffer, support::node_id(*handle));
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
/// let field = field.unwrap();
/// assert_eq!(field.to_string(), alignment_score);
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
    /// Returns an error if the field cannot be parsed or the type is unsupported.
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

    /// Appends the field to the given buffer.
    ///
    /// If `as_new_field` is `true`, a tab character is added before the field.
    pub fn append_to(&self, buffer: &mut Vec<u8>, as_new_field: bool) {
        if as_new_field {
            buffer.push(b'\t');
        }
        let _ = write!(buffer, "{}", self);
    }

    /// Appends a typed string field to the given buffer.
    ///
    /// If `as_new_field` is `true`, a tab character is added before the field.
    /// This bypasses the creation of a `TypedField` object and therefore avoids copying the value.
    pub fn append_string(buffer: &mut Vec<u8>, tag: [u8; 2], value: &[u8], as_new_field: bool) {
        if as_new_field {
            buffer.push(b'\t');
        }
        buffer.push(tag[0]);
        buffer.push(tag[1]);
        buffer.push(b':');
        buffer.push(b'Z');
        buffer.push(b':');
        buffer.extend_from_slice(value);
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

/// Returns `true` if the buffer contains a GAF header line.
///
/// A header line starts with `@`.
pub fn is_gaf_header_line(buf: &[u8]) -> bool {
    buf.first() == Some(&b'@')
}

/// Returns `true` if the next line in the reader would be a GAF header line.
///
/// The reader position is not changed.
/// Returns an I/O error if reading from the reader fails.
/// A header line starts with `@`.
pub fn peek_gaf_header_line<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    let buffer = reader.fill_buf()?;
    Ok(buffer.first() == Some(&b'@'))
}

/// Returns all successive GAF header lines from the reader.
///
/// The returned lines do not contain the trailing newline character.
/// The reader position is advanced past the header lines.
/// Returns an I/O error if reading from the reader fails.
///
/// # Examples
///
/// ```
/// use gbz_base::{formats, utils};
///
/// let filename = utils::get_test_data("good.gaf");
/// let mut reader = utils::open_file(&filename)
///     .expect("Unable to open test file");
/// let headers = formats::read_gaf_header_lines(&mut reader)
///     .expect("Unable to read GAF header lines");
/// assert_eq!(headers.len(), 2);
/// assert!(headers[0].starts_with("@HD")); // File header
/// assert!(headers[1].starts_with("@RN")); // Reference name
/// ```
pub fn read_gaf_header_lines<R: BufRead>(reader: &mut R) -> io::Result<Vec<String>> {
    let mut headers: Vec<String> = Vec::new();
    while peek_gaf_header_line(reader)? {
        let mut line: Vec<u8> = Vec::new();
        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        headers.push(String::from_utf8_lossy(&line).to_string());
    }
    Ok(headers)
}

/// Writes a GAF file header.
pub fn write_gaf_file_header<T: Write>(output: &mut T) -> io::Result<()> {
    let header = format!("@HD\tVN:Z:1.0\n");
    output.write_all(header.as_bytes())?;
    Ok(())
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

/// Writes the given header lines.
///
/// If the lines do not end with a newline character, one is added.
pub fn write_header_lines<T: Write>(header_lines: &[String], output: &mut T) -> io::Result<()> {
    for line in header_lines {
        if line.ends_with('\n') {
            output.write_all(line.as_bytes())?;
        } else {
            output.write_all(line.as_bytes())?;
            output.write_all(b"\n")?;
        }
    }
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
    utils::append_usize(&mut buffer, metadata.name.haplotype);
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.name.contig.as_bytes());
    buffer.push(b'\t');
    utils::append_usize(&mut buffer, metadata.name.fragment);
    buffer.push(b'\t');
    utils::append_usize(&mut buffer, metadata.end);
    buffer.push(b'\t');
    append_walk(&mut buffer, path);
    if let Some(weight) = metadata.weight {
        let field = TypedField::Int([b'W', b'T'], weight as isize);
        field.append_to(&mut buffer, true);
    }
    if let Some(cigar) = &metadata.cigar {
        let field = TypedField::String([b'C', b'G'], cigar.as_bytes().to_vec());
        field.append_to(&mut buffer, true);
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
