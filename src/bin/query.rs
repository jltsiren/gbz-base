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

    // Find the query position on the reference path.
    let (query_pos, gbwt_pos) = query_position(
        &mut cursor,
        &ref_info,
        config.offset
    )?;

    // Extract all GBWT records within the context.
    let subgraph = extract_context(&mut cursor, query_pos, config.context)?;

    // Extract paths.
    let (paths, (ref_id, path_offset)) = extract_paths(&subgraph, gbwt_pos)?;
    let ref_start = config.offset - distance_to(&subgraph, &paths[ref_id].0, path_offset, query_pos.offset);
    let ref_interval = ref_start..ref_start + paths[ref_id].1;

    // GFA output.
    let mut output = io::stdout();
    let reference_samples = cursor.get_gbwt_tag(REFERENCE_SAMPLES_KEY)?;
    write_gfa(&subgraph, reference_samples, &mut output).map_err(|x| x.to_string())?;
    let ref_metadata = WalkMetadata::from_gbz_path(&ref_info, ref_interval);
    write_gfa_walk(&paths[ref_id].0, &ref_metadata, &mut output).map_err(|x| x.to_string())?;
    let mut haplotype = 1;
    for (id, (path, len)) in paths.iter().enumerate() {
        if id == ref_id {
            continue;
        }
        let metadata = WalkMetadata::anonymous(haplotype, &ref_info.contig, *len);
        write_gfa_walk(path, &metadata, &mut output).map_err(|x| x.to_string())?;
        haplotype += 1;
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

// Returns the graph position and the GBWT position for the given offset.
fn query_position(cursor: &mut Cursor, path: &GBZPath, query_offset: usize) -> Result<(GraphPosition, Pos), String> {
    let result = cursor.indexed_position(path.handle, query_offset)?;
    let (mut path_offset, mut pos) = result.ok_or(format!("Path {} is not indexed", path.name()))?;

    let mut graph_pos: Option<GraphPosition> = None;
    let mut gbwt_pos: Option<Pos> = None;
    loop {
        let handle = pos.node;
        let record = cursor.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset + record.sequence.len() > query_offset {
            graph_pos = Some(GraphPosition {
                node: support::node_id(handle),
                orientation: support::node_orientation(handle),
                offset: query_offset - path_offset,
            });
            gbwt_pos = Some(pos);
            break;
        }
        path_offset += record.sequence.len();
        let gbwt_record = record.to_gbwt_record().ok_or(
            format!("The record for handle {} is invalid", handle)
        )?;
        let next = gbwt_record.lf(pos.offset);
        if next.is_none() {
            break;
        }
        pos = next.unwrap();
    }

    let graph_pos = graph_pos.ok_or(
        format!("Path {} does not contain offset {}", path.name(), query_offset)
    )?;
    let gbwt_pos = gbwt_pos.unwrap();
    Ok((graph_pos, gbwt_pos))
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

    // If both ends are in the same orientation, the forward orientation is canonical.
    let first = path[0];
    let last = path[path.len() - 1];
    let first_o = support::node_orientation(first);
    let last_o = support::node_orientation(last);
    if first_o == last_o {
        return first_o == Orientation::Forward;
    }

    // Otherwise we consider the path a single edge.
    edge_is_canonical(
        (support::node_id(first), first_o),
        (support::node_id(last), last_o)
    )
}

// Extract the handle sequences and lengths for all paths. The second return
// value is (offset in result, offset on that path) for the handle corresponding
// to `ref_pos`.
fn extract_paths(
    subgraph: &BTreeMap<usize, GBZRecord>,
    ref_pos: Pos
) -> Result<(Vec<(Vec<usize>, usize)>, (usize, usize)), String> {
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
    // Extract all paths and note if one of them passes through `ref_pos`.
    let mut result: Vec<(Vec<usize>, usize)> = Vec::new();
    let mut ref_id_offset: Option<(usize, usize)> = None;
    for (handle, positions) in successors.iter() {
        for (offset, (_, has_predecessor)) in positions.iter().enumerate() {
            if *has_predecessor {
                continue;
            }
            let mut curr = Some(Pos::new(*handle, offset));
            let mut is_ref = false;
            let mut path: Vec<usize> = Vec::new();
            let mut len = 0;
            while let Some(pos) = curr {
                if pos == ref_pos {
                    ref_id_offset = Some((result.len(), path.len()));
                    is_ref = true;
                }
                path.push(pos.node);
                len += subgraph.get(&pos.node).unwrap().sequence.len();
                curr = next_pos(pos, &successors);
            }
            if is_ref {
                if !path_is_canonical(&path) {
                    eprintln!("Warning: the reference path is not in canonical orientation");
                }
                result.push((path, len));
            } else if path_is_canonical(&path) {
                result.push((path, len));
            }
        }
    }

    let ref_id_offset = ref_id_offset.ok_or(format!("Could not find the reference path"))?;
    Ok((result, ref_id_offset))
}

fn distance_to(subgraph: &BTreeMap<usize, GBZRecord>, path: &[usize], path_offset: usize, node_offset: usize) -> usize {
    let mut result = node_offset;
    for handle in path.iter().take(path_offset) {
        result += subgraph.get(handle).unwrap().sequence.len();
    }
    result
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
