use gbz_base::{Cursor, GBZBase, GBZRecord, GBZPath};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};
use std::io::Write;
use std::{env, io, process};

use gbwt::{Orientation, REFERENCE_SAMPLES_KEY, ENDMARKER};
use getopts::Options;

use gbwt::support;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    // Parse arguments.
    let config = Config::new()?;

    // Open the database.
    let database = GBZBase::open(&config.filename)?;
    let mut cursor = Cursor::new(&database)?;

    // Find the reference path.
    let ref_path = cursor.find_path(
        &config.sample, &config.contig, config.haplotype, config.fragment
    )?;
    let ref_info = ref_path.ok_or(format!("Cannot find path {}", config.path_name()))?;

    // Find an indexed position for offset - context.
    let low = if config.offset >= config.context { config.offset - config.context } else { 0 };
    let result = cursor.indexed_position(ref_info.handle, low)?;
    let (mut initial_offset, mut pos) = result.ok_or(format!("Path {} is not indexed", ref_info.name()))?;

    // Walk the path until offset + context and find the node containing offset.
    let high = config.offset + config.context;
    let mut path_offset = initial_offset;
    let mut initial_node: Option<usize> = None;
    let mut initial_node_offset = 0;
    let mut ref_path: Vec<usize> = Vec::new();
    while path_offset < high {
        let handle = pos.node;
        let record = cursor.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset <= low {
            initial_offset = path_offset;
        }
        if path_offset <= config.offset {
            initial_node_offset = config.offset - path_offset;
        }
        path_offset += record.sequence.len();
        if path_offset > low {
            ref_path.push(handle);
            if path_offset > config.offset && initial_node.is_none() {
                initial_node = Some(handle);
            }
        }
        let gbwt_record = record.to_gbwt_record().ok_or(
            format!("The record for handle {} is invalid", handle)
        )?;
        let next = gbwt_record.lf(pos.offset);
        if next.is_none() {
            break;
        }
        pos = next.unwrap();
    }
    let initial_node = initial_node.ok_or(
        format!("Path {} does not contain offset {}", ref_info.name(), config.offset)
    )?;
    let initial_orientation = support::node_orientation(initial_node);
    let initial_node = support::node_id(initial_node);

    // Start graph traversal from the initial node.
    let mut active: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::new(); // (distance, node id)
    active.push(Reverse((0, initial_node)));

    // Traverse in both directions.
    let mut selected: BTreeMap<usize, GBZRecord> = BTreeMap::new();
    while !active.is_empty() {
        let (distance, node_id) = active.pop().unwrap().0;
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let handle = support::encode_node(node_id, orientation);
            if selected.contains_key(&handle) {
                continue;
            }
            let record = cursor.get_record(handle)?;
            let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
            let next_distance = if node_id == initial_node {
                distance_to_end(&record, initial_orientation, initial_node_offset)
            } else {
                distance + record.sequence.len()
            };
            if next_distance <= config.context {
                for edge in record.edges.iter() {
                    if edge.node != ENDMARKER && !selected.contains_key(&edge.node) {
                        active.push(Reverse((next_distance, support::node_id(edge.node))));
                    }
                }
            }
            selected.insert(handle, record);
        }
    }

    // GFA output.
    let mut output = io::stdout();
    let reference_samples = cursor.get_gbwt_tag(REFERENCE_SAMPLES_KEY)?;
    write_gfa(&selected, reference_samples, &mut output).map_err(|x| x.to_string())?;
    write_gfa_walk(&ref_path, &ref_info, initial_offset, path_offset, &mut output).map_err(|x| x.to_string())?;

    Ok(())
}

//-----------------------------------------------------------------------------

// TODO: haplotype, fragment; offset or interval
pub struct Config {
    pub filename: String,
    pub sample: String,
    pub contig: String,
    pub haplotype: usize,
    pub fragment: usize,
    pub offset: usize,
    pub context: usize,
}

impl Config {
    pub fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("s", "sample", "sample name (required)", "STR");
        opts.optopt("c", "contig", "contig name (required)", "STR");
        opts.optopt("o", "offset", "sequence offset (required)", "INT");
        opts.optopt("n", "context", "context length in bp (default 100)", "INT");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let haplotype: usize = 0;
        let fragment: usize = 0;
        let mut offset: Option<usize> = None;
        let mut context: usize = 100; // FIXME source for this
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz.db", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        let sample = matches.opt_str("s");
        let contig = matches.opt_str("c");
        if let Some(s) = matches.opt_str("o") {
            offset = Some(s.parse::<usize>().map_err(|x| format!("--offset: {}", x))?);
        }
        if let Some(s) = matches.opt_str("n") {
            context = s.parse::<usize>().map_err(|x| format!("--context: {}", x))?;
        }

        let filename = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            let header = format!("Usage: {} [options] graph.gbz.db", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };

        Ok(Config {
            filename,
            sample: sample.ok_or("Sample name must be provided with --sample".to_string())?,
            contig: contig.ok_or("Contig name must be provided with --contig".to_string())?,
            haplotype,
            fragment,
            offset: offset.ok_or("Sequence offset must be provided with --offset".to_string())?,
            context,
        })
    }

    pub fn path_name(&self) -> String {
        format!("{}#{}#{}@{}", self.sample, self.haplotype, self.contig, self.fragment)
    }
}

//-----------------------------------------------------------------------------

fn distance_to_end(record: &GBZRecord, orientation: Orientation, offset: usize) -> usize {
    if orientation == support::node_orientation(record.handle) {
        record.sequence.len() - offset
    } else {
        offset + 1
    }
}

//-----------------------------------------------------------------------------

fn write_gfa<T: Write>(records: &BTreeMap<usize, GBZRecord>, reference_samples: Option<String>, output: &mut T) -> io::Result<()> {
    write_gfa_header(reference_samples, output)?;

    // Segments.
    for (handle, record) in records.iter() {
        if support::node_orientation(*handle) == Orientation::Forward {
            write_gfa_segment(record, output)?;
        }
    }

    // Links.
    for (handle, record) in records.iter() {
        let (from_id, from_o) = support::decode_node(*handle);
        for edge in record.edges.iter() {
            let (to_id, to_o) = support::decode_node(edge.node);
            if records.contains_key(&edge.node) && is_canonical_orientation((from_id, from_o), (to_id, to_o)) {
                write_gfa_link(
                    (from_id.to_string().as_bytes(), from_o),
                    (to_id.to_string().as_bytes(), to_o),
                    output
                )?;
            }
        }
    }

    // FIXME paths, walks

    Ok(())
}

// TODO: These GFA writing support functions should be shared with gbunzip in gbwt-rs.

fn write_gfa_header<T: Write>(reference_samples: Option<String>, output: &mut T) -> io::Result<()> {
    let header = if let Some(sample_names) = reference_samples {
        format!("H\tVN:Z:1.1\tRS:Z:{}\n", sample_names)
    } else {
        "H\tVN:Z:1.1\n".to_string()
    };
    output.write_all(header.as_bytes())?;
    Ok(())
}

fn write_gfa_segment<T: Write>(record: &GBZRecord, output: &mut T) -> io::Result<()> {
    let (id, orientation) = support::decode_node(record.handle);
    output.write_all(b"S\t")?;
    output.write_all(id.to_string().as_bytes())?;
    output.write_all(b"\t")?;
    if orientation == Orientation::Reverse {
        output.write_all(&record.sequence)?;
    } else {
        let rc = support::reverse_complement(&record.sequence);
        output.write_all(&rc)?;
    }
    output.write_all(b"\n")?;
    Ok(())
}

// TODO this should be a public function
fn is_canonical_orientation(from: (usize, Orientation), to: (usize, Orientation)) -> bool {
    if from.1 == Orientation::Forward {
        to.0 >= from.0
    } else {
        (to.0 > from.0) || (to.0 == from.0 && to.1 == Orientation::Forward)
    }
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
        Orientation::Forward => output.write_all(b"\t+\t*\n")?,
        Orientation::Reverse => output.write_all(b"\t-\t*\n")?,
    }
    Ok(())
}

fn write_gfa_walk<T: Write>(path: &[usize], info: &GBZPath, from: usize, to: usize, output: &mut T) -> io::Result<()> {
    let mut buffer: Vec<u8> = Vec::new();
    buffer.push(b'W');
    buffer.push(b'\t');
    buffer.extend_from_slice(info.sample.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(info.haplotype.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(info.contig.as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(from.to_string().as_bytes());
    buffer.push(b'\t');
    buffer.extend_from_slice(to.to_string().as_bytes());
    buffer.push(b'\t');
    for handle in path.iter() {
        match support::node_orientation(*handle) {
            Orientation::Forward => buffer.push(b'>'),
            Orientation::Reverse => buffer.push(b'<'),
        }
        buffer.extend_from_slice(support::node_id(*handle).to_string().as_bytes());
    }
    buffer.push(b'\n');
    output.write_all(&buffer)?;
    Ok(())
}

//-----------------------------------------------------------------------------
