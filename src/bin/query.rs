use gbz_base::{Cursor, GBZBase, GBZRecord, GBZPath};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};
use std::io::Write;
use std::ops::Range;
use std::{env, io, process};

use gbwt::{Orientation, Pos, REFERENCE_SAMPLES_KEY, ENDMARKER};
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
    if !ref_info.is_reference {
        return Err(format!("Path {} is not a reference path", config.path_name()));
    }

    // Extract the reference path fragment.
    let (query_pos, ref_path, ref_interval, ref_pos) = extract_path(
        &mut cursor,
        &ref_info,
        config.offset,
        config.context
    )?;

    // Extract all GBWT records within the context.
    let subgraph = extract_context(&mut cursor, query_pos, config.context)?;

    // Extract all other paths.
    let other_paths = extract_other_paths(&subgraph, ref_pos)?;

    // GFA output.
    let mut output = io::stdout();
    let reference_samples = cursor.get_gbwt_tag(REFERENCE_SAMPLES_KEY)?;
    write_gfa(&subgraph, reference_samples, &mut output).map_err(|x| x.to_string())?;
    let ref_metadata = WalkMetadata::from_gbz_path(&ref_info, ref_interval);
    write_gfa_walk(&ref_path, &ref_metadata, &mut output).map_err(|x| x.to_string())?;
    for (haplotype, (path, len)) in other_paths.iter().enumerate() {
        let metadata = WalkMetadata::anonymous(haplotype + 1, &ref_info.contig, *len);
        write_gfa_walk(path, &metadata, &mut output).map_err(|x| x.to_string())?;
    }

    Ok(())
}

//-----------------------------------------------------------------------------

// TODO: haplotype, fragment; offset or interval
// TODO: ref path only, ref + distinct haplotypes, ref + all haplotypes
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

// TODO: Should this be in the library?
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct GraphPosition {
    pub node: usize,
    pub orientation: Orientation,
    pub offset: usize,
}

//-----------------------------------------------------------------------------

fn query_interval(offset: usize, context: usize) -> Range<usize> {
    let low = if offset >= context { offset - context } else { 0 };
    let high = offset + context;
    low..high
}

// FIXME: This should only return the query position and the GBWT position.
// FIXME: Then we can extract all paths later and determine which of them is the reference.
// Returns the query position, the sequence of handles, the path interval, and the GBWT position for path start.
fn extract_path(
    cursor: &mut Cursor,
    path: &GBZPath,
    offset: usize,
    context: usize
) -> Result<(GraphPosition, Vec<usize>, Range<usize>, Pos), String> {
    // Find an indexed position for the start of the query interval.
    let query_interval = query_interval(offset, context);
    let result = cursor.indexed_position(path.handle, query_interval.start)?;
    let (mut path_offset, mut pos) = result.ok_or(format!("Path {} is not indexed", path.name()))?;

    // Walk the path until offset + context and find the node containing offset.
    let mut path_handles: Vec<usize> = Vec::new();
    let mut query_pos: Option<GraphPosition> = None;
    let mut path_interval = 0..0;
    let mut start_pos = pos;
    while path_offset < query_interval.end {
        let handle = pos.node;
        let record = cursor.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset <= query_interval.start {
            path_interval.start = path_offset;
            start_pos = pos;
        }
        path_offset += record.sequence.len();
        path_interval.end = path_offset;
        if path_offset > query_interval.start {
            path_handles.push(handle);
            if path_offset > offset && query_pos.is_none() {
                query_pos = Some(GraphPosition {
                    node: support::node_id(handle),
                    orientation: support::node_orientation(handle),
                    offset: offset - (path_offset - record.sequence.len())
                });
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

    let query_pos = query_pos.ok_or(
        format!("Path {} does not contain offset {}", path.name(), offset)
    )?;
    Ok((query_pos, path_handles, path_interval, start_pos))
}

//-----------------------------------------------------------------------------

fn distance_to_end(record: &GBZRecord, orientation: Orientation, offset: usize) -> usize {
    if orientation == support::node_orientation(record.handle) {
        record.sequence.len() - offset
    } else {
        offset + 1
    }
}

fn extract_context(
    cursor: &mut Cursor,
    from: GraphPosition,
    context: usize
) -> Result<BTreeMap<usize, GBZRecord>, String> {
    // Start graph traversal from the initial node.
    let mut active: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::new(); // (distance, node id)
    active.push(Reverse((0, from.node)));

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
            let next_distance = if node_id == from.node {
                distance_to_end(&record, from.orientation, from.offset)
            } else {
                distance + record.sequence.len()
            };
            if next_distance <= context {
                for edge in record.edges.iter() {
                    if edge.node != ENDMARKER && !selected.contains_key(&edge.node) {
                        active.push(Reverse((next_distance, support::node_id(edge.node))));
                    }
                }
            }
            selected.insert(handle, record);
        }
    }

    Ok(selected)
}

//-----------------------------------------------------------------------------

fn next_pos(pos: Pos, successors: &BTreeMap<usize, Vec<(Pos, bool)>>) -> Option<Pos> {
    if let Some(v) = successors.get(&pos.node) {
        let (next, _) = v[pos.offset];
        if next.node == ENDMARKER || !successors.contains_key(&next.node) {
            None
        } else {
            Some(next)
        }
    } else {
        None
    }
}

// TODO this should be a public function
fn edge_is_canonical(from: (usize, Orientation), to: (usize, Orientation)) -> bool {
    if from.1 == Orientation::Forward {
        to.0 >= from.0
    } else {
        (to.0 > from.0) || (to.0 == from.0 && to.1 == Orientation::Forward)
    }
}

// TODO this should be a public function
fn path_is_canonical(path: &[usize]) -> bool {
    if path.is_empty() {
        return true;
    }
    if path.len() == 1 && support::node_orientation(path[0]) == Orientation::Forward {
        return true;
    }

    // In the general case, we consider the path a single edge.
    let first = path[0];
    let last = path[path.len() - 1];
    edge_is_canonical(
        (support::node_id(first), support::node_orientation(first)),
        (support::node_id(last), support::node_orientation(last))
    )
}

// Extract the handle sequences and lengths for all paths in the subgraph except
// the one passing through the given GBWT position.
// FIXME this should only return distinct paths
fn extract_other_paths(
    subgraph: &BTreeMap<usize, GBZRecord>,
    ref_pos: Pos
) -> Result<Vec<(Vec<usize>, usize)>, String> {
    // Decompress the GBWT node records for the subgraph.
    let mut keys: Vec<usize> = Vec::new();
    let mut successors: BTreeMap<usize, Vec<(Pos, bool)>> = BTreeMap::new();
    for (handle, gbz_record) in subgraph.iter() {
        let gbwt_record = gbz_record.to_gbwt_record().unwrap();
        let decompressed: Vec<(Pos, bool)> = gbwt_record.decompress().into_iter().map(|x| (x, false)).collect();
        keys.push(*handle);
        successors.insert(*handle, decompressed);
    }

    // Mark the positions that have predecessors in the subgraph.
    for handle in keys.iter() {
        let decompressed = successors.get(handle).unwrap().clone();
        for (pos, _) in decompressed.iter() {
            if let Some(v) = successors.get_mut(&pos.node) {
                v[pos.offset].1 = true;
            }
        }
    }

    // FIXME: Check for infinite loops.
    // Extract all paths that do not pass through `ref_pos`.
    let mut result: Vec<(Vec<usize>, usize)> = Vec::new();
    for (handle, positions) in successors.iter() {
        for (offset, (_, has_predecessor)) in positions.iter().enumerate() {
            if *has_predecessor {
                continue;
            }
            let mut curr = Some(Pos::new(*handle, offset));
            let mut ok = true;
            let mut path: Vec<usize> = Vec::new();
            let mut len = 0;
            while let Some(pos) = curr {
                if pos == ref_pos {
                    ok = false;
                    break;
                }
                path.push(pos.node);
                len += subgraph.get(&pos.node).unwrap().sequence.len();
                curr = next_pos(pos, &successors);
            }
            if ok && path_is_canonical(&path) {
                result.push((path, len));
            }
        }
    }

    Ok(result)
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
            if records.contains_key(&edge.node) && edge_is_canonical((from_id, from_o), (to_id, to_o)) {
                write_gfa_link(
                    (from_id.to_string().as_bytes(), from_o),
                    (to_id.to_string().as_bytes(), to_o),
                    output
                )?;
            }
        }
    }

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

struct WalkMetadata {
    sample: String,
    haplotype: usize,
    contig: String,
    interval: Range<usize>,
}

impl WalkMetadata {
    fn from_gbz_path(path: &GBZPath, interval: Range<usize>) -> Self {
        WalkMetadata {
            sample: path.sample.clone(),
            haplotype: path.haplotype,
            contig: path.contig.clone(),
            interval,
        }
    }

    fn anonymous(haplotype: usize, contig: &str, len: usize) -> Self {
        WalkMetadata { sample: "unknown".to_string(), haplotype, contig: contig.to_owned(), interval: 0..len }
    }
}

fn write_gfa_walk<T: Write>(path: &[usize], metadata: &WalkMetadata, output: &mut T) -> io::Result<()> {
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
    buffer.push(b'\n');
    output.write_all(&buffer)?;
    Ok(())
}

//-----------------------------------------------------------------------------
