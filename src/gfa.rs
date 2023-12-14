use crate::{GBZRecord, GBZPath};

use std::collections::BTreeMap;
use std::io;
use std::io::Write;
use std::ops::Range;

use gbwt::support;
use gbwt::Orientation;

//-----------------------------------------------------------------------------

pub fn write_gfa<T: Write>(records: &BTreeMap<usize, GBZRecord>, reference_samples: Option<String>, output: &mut T) -> io::Result<()> {
    write_gfa_header(reference_samples, output)?;

    // Segments.
    for (handle, record) in records.iter() {
        if support::node_orientation(*handle) == Orientation::Forward {
            write_gfa_segment(record, output)?;
        }
    }

    // Links.
    for (handle, record) in records.iter() {
        let from = support::decode_node(*handle);
        for successor in record.successors() {
            let to = support::decode_node(successor);
            if records.contains_key(&successor) && support::edge_is_canonical(from, to) {
                write_gfa_link(
                    (from.0.to_string().as_bytes(), from.1),
                    (to.0.to_string().as_bytes(), to.1),
                    output
                )?;
            }
        }
    }

    Ok(())
}

// TODO: These GFA writing support functions should be shared with gbunzip in gbwt-rs.

pub fn write_gfa_header<T: Write>(reference_samples: Option<String>, output: &mut T) -> io::Result<()> {
    let header = if let Some(sample_names) = reference_samples {
        format!("H\tVN:Z:1.1\tRS:Z:{}\n", sample_names)
    } else {
        "H\tVN:Z:1.1\n".to_string()
    };
    output.write_all(header.as_bytes())?;
    Ok(())
}

fn write_gfa_segment<T: Write>(record: &GBZRecord, output: &mut T) -> io::Result<()> {
    let (id, orientation) = support::decode_node(record.handle());
    output.write_all(b"S\t")?;
    output.write_all(id.to_string().as_bytes())?;
    output.write_all(b"\t")?;
    if orientation == Orientation::Reverse {
        output.write_all(record.sequence())?;
    } else {
        let rc = support::reverse_complement(record.sequence());
        output.write_all(&rc)?;
    }
    output.write_all(b"\n")?;
    Ok(())
}

fn write_gfa_link<T: Write>(from: (&[u8], Orientation), to: (&[u8], Orientation), output: &mut T) -> io::Result<()> {
    output.write_all(b"L\t")?;
    output.write_all(from.0)?;
    match from.1 {
        Orientation::Forward => output.write_all(b"\t+\t")?,
        Orientation::Reverse => output.write_all(b"\t-\t")?,
    }
    output.write_all(to.0)?;
    match to.1 {
        Orientation::Forward => output.write_all(b"\t+\t0M\n")?,
        Orientation::Reverse => output.write_all(b"\t-\t0M\n")?,
    }
    Ok(())
}

pub struct WalkMetadata {
    sample: String,
    haplotype: usize,
    contig: String,
    interval: Range<usize>,
    weight: Option<usize>,
}

impl WalkMetadata {
    pub fn from_gbz_path(path: &GBZPath, interval: Range<usize>, weight: Option<usize>) -> Self {
        WalkMetadata {
            sample: path.sample.clone(),
            haplotype: path.haplotype,
            contig: path.contig.clone(),
            interval,
            weight,
        }
    }

    pub fn anonymous(haplotype: usize, contig: &str, len: usize, weight: Option<usize>) -> Self {
        WalkMetadata {
            sample: "unknown".to_string(),
            haplotype,
            contig: contig.to_owned(),
            interval: 0..len,
            weight,
        }
    }
}

pub fn write_gfa_walk<T: Write>(path: &[usize], metadata: &WalkMetadata, output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();
    buffer.push(b'W');
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.sample.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.haplotype.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.contig.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.interval.start.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(metadata.interval.end.to_string().as_bytes());
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
    buffer.push(b'\n');
    output.write_all(&buffer)?;
    Ok(())
}
