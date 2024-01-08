use gbz_base::{GBZBase, GBZRecord, GBZPath, GraphInterface, GBZPathName};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};
use std::fmt::Display;
use std::io::Write;
use std::ops::Range;
use std::time::Instant;
use std::{env, io, process};

use gbwt::{Orientation, Pos, ENDMARKER};
use getopts::Options;

use gbwt::support;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = Config::new()?;

    // Open the database.
    let database = GBZBase::open(&config.filename)?;
    let mut graph = GraphInterface::new(&database)?;

    // Extract the subgraph.
    let subgraph = Subgraph::new(
        &mut graph,
        &config.path_name, config.offset, config.context,
        config.output
    )?;

    // Write the output.
    let mut output = io::stdout();
    match config.format {
        OutputFormat::GFA => subgraph.write_gfa(&mut output).map_err(|x| x.to_string())?,
        OutputFormat::JSON => subgraph.write_json(&mut output).map_err(|x| x.to_string())?,
    }

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
}

//-----------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum HaplotypeOutput {
    All,
    Distinct,
    ReferenceOnly,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum OutputFormat {
    GFA,
    JSON,
}

struct Config {
    filename: String,
    path_name: GBZPathName,
    offset: usize,
    context: usize,
    output: HaplotypeOutput,
    format: OutputFormat,
}

impl Config {
    // Default context length in bp.
    const DEFAULT_CONTEXT: usize = 100;

    fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("s", "sample", "sample name (required)", "STR");
        opts.optopt("c", "contig", "contig name (required)", "STR");
        opts.optopt("o", "offset", "sequence offset (-o or -i is required)", "INT");
        opts.optopt("i", "interval", "sequence interval (-o or -i is required)", "INT..INT");
        let context_desc = format!("context length in bp (default: {})", Self::DEFAULT_CONTEXT);
        opts.optopt("n", "context", &context_desc, "INT");
        opts.optflag("d", "distinct", "output distinct haplotypes with weights");
        opts.optflag("r", "reference-only", "output the reference but no other haplotypes");
        opts.optopt("f", "format", "output format (gfa or json, default: gfa)", "STR");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let filename = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            let header = format!("Usage: {} [options] graph.gbz.db", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };

        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz.db", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }

        let sample = matches.opt_str("s").ok_or("Sample name must be provided with --sample".to_string())?;
        let contig = matches.opt_str("c").ok_or("Contig name must be provided with --contig".to_string())?;
        let path_name = GBZPathName::reference(&sample, &contig);
        let (offset, context) = Self::parse_offset_and_context(&matches)?;

        let mut output = HaplotypeOutput::All;
        if matches.opt_present("d") {
            output = HaplotypeOutput::Distinct;
        }
        if matches.opt_present("r") {
            output = HaplotypeOutput::ReferenceOnly;
        }

        let mut format: OutputFormat = OutputFormat::GFA;
        if let Some(s) = matches.opt_str("f") {
            match s.to_lowercase().as_str() {
                "gfa" => format = OutputFormat::GFA,
                "json" => format = OutputFormat::JSON,
                _ => return Err(format!("Invalid output format: {}", s)),
            }
        }

        Ok(Config { filename, path_name, offset, context, output, format, })
    }

    fn parse_interval(s: &str) -> Result<Range<usize>, String> {
        let mut parts = s.split("..");
        let start = parts.next().ok_or(format!("Invalid interval: {}", s))?;
        let start = start.parse::<usize>().map_err(|x| format!("{}", x))?;
        let end = parts.next().ok_or(format!("Invalid interval: {}", s))?;
        let end = end.parse::<usize>().map_err(|x| format!("{}", x))?;
        if parts.next().is_some() {
            return Err(format!("Invalid interval: {}", s));
        }
        Ok(start..end)
    }

    fn parse_offset_and_context(matches: &getopts::Matches) -> Result<(usize, usize), String> {
        let mut context = if let Some(s) = matches.opt_str("n") {
            s.parse::<usize>().map_err(|x| format!("--context: {}", x))?
        } else {
            Self::DEFAULT_CONTEXT
        };

        let mut offset: Option<usize> = None;
        let mut interval: Option<Range<usize>> = None;
        if let Some(s) = matches.opt_str("o") {
            offset = Some(s.parse::<usize>().map_err(|x| format!("--offset: {}", x))?);
        }
        if let Some(s) = matches.opt_str("i") {
            let parsed = Self::parse_interval(&s)?;
            if parsed.is_empty() {
                return Err("--interval: Sequence interval cannot be empty".to_string());
            }
            interval = Some(parsed);
        }
        if offset.is_none() && interval.is_none() {
            return Err("Either --offset or --interval must be provided".to_string());
        }
        if offset.is_some() && interval.is_some() {
            return Err("Only one of --offset and --interval can be provided".to_string());
        }
        if let Some(interval) = interval {
            context += interval.len() / 2;
            offset = Some(interval.start + interval.len() / 2);
        }
        let offset = offset.unwrap();

        Ok((offset, context))
    }
}

//-----------------------------------------------------------------------------

// TODO: Should this be in the library?
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct GraphPosition {
    node: usize,
    orientation: Orientation,
    offset: usize,
}

//-----------------------------------------------------------------------------

// TODO: Maybe this should be in the library.
struct Subgraph {
    records: BTreeMap<usize, GBZRecord>,
    paths: Vec<PathInfo>,
    ref_id: usize,
    ref_path: GBZPath,
    ref_interval: Range<usize>,
}

impl Subgraph {
    fn new(graph: &mut GraphInterface, path_name: &GBZPathName, offset: usize, context: usize, haplotype_output: HaplotypeOutput) -> Result<Self, String> {
        // Find the reference path.
        let ref_path = graph.find_path(path_name)?;
        let ref_path = ref_path.ok_or(format!("Cannot find path {}", path_name))?;
        if !ref_path.is_indexed {
            return Err(format!("Path {} has not been indexed for random access", path_name));
        }

        let (query_pos, gbwt_pos) = query_position(graph, &ref_path, offset)?;
        let records = extract_context(graph, query_pos, context)?;
        let (mut paths, (mut ref_id, path_offset)) = extract_paths(&records, gbwt_pos)?;
        let ref_start = offset - distance_to(&records, &paths[ref_id].path, path_offset, query_pos.offset);
        let ref_interval = ref_start..ref_start + paths[ref_id].len;

        if haplotype_output == HaplotypeOutput::Distinct {
            // TODO: We should use bidirectional search in GBWT to find the distinct paths directly.
            (paths, ref_id) = make_distinct(paths, ref_id);
        } else if haplotype_output == HaplotypeOutput::ReferenceOnly {
            paths = vec![paths[ref_id].clone()];
            ref_id = 0;
        }

        Ok(Subgraph { records, paths, ref_id, ref_path, ref_interval, })
    }

    fn write_gfa<T: Write>(&self, output: &mut T) -> io::Result<()> {
        // Header.
        let reference_samples = Some(self.ref_path.name.sample.clone());
        write_gfa_header(reference_samples, output)?;

        // Segments.
        for (handle, record) in self.records.iter() {
            if support::node_orientation(*handle) == Orientation::Forward {
                write_gfa_segment(record, output)?;
            }
        }

        // Links.
        for (handle, record) in self.records.iter() {
            let from = support::decode_node(*handle);
            for successor in record.successors() {
                let to = support::decode_node(successor);
                if self.records.contains_key(&successor) && support::edge_is_canonical(from, to) {
                    write_gfa_link(
                        (from.0.to_string().as_bytes(), from.1),
                        (to.0.to_string().as_bytes(), to.1),
                        output
                    )?;
                }
            }
        }

        // Paths.
        let ref_metadata = WalkMetadata::path_interval(
            &self.ref_path,
            self.ref_interval.clone(),
            self.paths[self.ref_id].weight
        );
        write_gfa_walk(&self.paths[self.ref_id].path, &ref_metadata, output)?;
        if self.paths.len() > 1 {
            let mut haplotype = 1;
            for (id, path_info) in self.paths.iter().enumerate() {
                if id == self.ref_id {
                    continue;
                }
                let metadata = WalkMetadata::anonymous(
                    haplotype,
                    &self.ref_path.name.contig,
                    path_info.len,
                    path_info.weight
                );
                write_gfa_walk(&path_info.path, &metadata, output)?;
                haplotype += 1;
            }
        }

        Ok(())
    }

    fn write_json<T: Write>(&self, output: &mut T) -> io::Result<()> {
        // Nodes.
        let mut nodes: Vec<JSONValue> = Vec::new();
        for (_, record) in self.records.iter() {
            let (id, orientation) = support::decode_node(record.handle());
            if orientation == Orientation::Reverse {
                continue;
            }
            let node = JSONValue::Object(vec![
                ("id".to_string(), JSONValue::String(id.to_string())),
                ("sequence".to_string(), JSONValue::String(String::from_utf8_lossy(record.sequence()).to_string())),
            ]);
            nodes.push(node);
        }

        // Edges.
        let mut edges: Vec<JSONValue> = Vec::new();
        for (handle, record) in self.records.iter() {
            let from = support::decode_node(*handle);
            for successor in record.successors() {
                let to = support::decode_node(successor);
                if self.records.contains_key(&successor) && support::edge_is_canonical(from, to) {
                    let edge = JSONValue::Object(vec![
                        ("from".to_string(), JSONValue::String(from.0.to_string())),
                        ("from_is_reverse".to_string(), JSONValue::Boolean(from.1 == Orientation::Reverse)),
                        ("to".to_string(), JSONValue::String(to.0.to_string())),
                        ("to_is_reverse".to_string(), JSONValue::Boolean(to.1 == Orientation::Reverse)),
                    ]);
                    edges.push(edge);
                }
            }
        }

        // Paths.
        let mut paths: Vec<JSONValue> = Vec::new();
        let ref_metadata = WalkMetadata::path_interval(
            &self.ref_path,
            self.ref_interval.clone(),
            self.paths[self.ref_id].weight
        );
        let ref_path = json_path(&self.paths[self.ref_id].path, &ref_metadata);
        paths.push(ref_path);
        let mut haplotype = 1;
        for (id, path_info) in self.paths.iter().enumerate() {
            if id == self.ref_id {
                continue;
            }
            let metadata = WalkMetadata::anonymous(
                haplotype,
                &self.ref_path.name.contig,
                path_info.len,
                path_info.weight
            );
            let path = json_path(&path_info.path, &metadata);
            paths.push(path);
            haplotype += 1;
        }

        let result = JSONValue::Object(vec![
            ("nodes".to_string(), JSONValue::Array(nodes)),
            ("edges".to_string(), JSONValue::Array(edges)),
            ("paths".to_string(), JSONValue::Array(paths)),
        ]);
        output.write_all(result.to_string().as_bytes())?;

        Ok(())
    }
}

//-----------------------------------------------------------------------------

// Returns the graph position and the GBWT position for the given offset.
fn query_position(graph: &mut GraphInterface, path: &GBZPath, query_offset: usize) -> Result<(GraphPosition, Pos), String> {
    let result = graph.indexed_position(path.handle, query_offset)?;
    let (mut path_offset, mut pos) = result.ok_or(format!("Path {} is not indexed", path.name()))?;

    let mut graph_pos: Option<GraphPosition> = None;
    let mut gbwt_pos: Option<Pos> = None;
    loop {
        let handle = pos.node;
        let record = graph.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset + record.sequence_len() > query_offset {
            graph_pos = Some(GraphPosition {
                node: support::node_id(handle),
                orientation: support::node_orientation(handle),
                offset: query_offset - path_offset,
            });
            gbwt_pos = Some(pos);
            break;
        }
        path_offset += record.sequence_len();
        let gbwt_record = record.to_gbwt_record();
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
    if orientation == support::node_orientation(record.handle()) {
        record.sequence_len() - offset
    } else {
        offset + 1
    }
}

fn extract_context(
    graph: &mut GraphInterface,
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
            let record = graph.get_record(handle)?;
            let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
            let next_distance = if node_id == from.node {
                distance_to_end(&record, from.orientation, from.offset)
            } else {
                distance + record.sequence_len()
            };
            if next_distance <= context {
                for successor in record.successors() {
                    if !selected.contains_key(&successor) {
                        active.push(Reverse((next_distance, support::node_id(successor))));
                    }
                }
            }
            selected.insert(handle, record);
        }
    }

    Ok(selected)
}

//-----------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PathInfo {
    path: Vec<usize>,
    len: usize,
    weight: Option<usize>,
}

impl PathInfo {
    fn new(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: None }
    }

    fn weighted(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: Some(1) }
    }

    fn increment_weight(&mut self) {
        if let Some(weight) = self.weight {
            self.weight = Some(weight + 1);
        }
    }
}

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

// Extract all paths in the subgraph. The second return value is
// (offset in result, offset on that path) for the handle corresponding to `ref_pos`.
fn extract_paths(
    subgraph: &BTreeMap<usize, GBZRecord>,
    ref_pos: Pos
) -> Result<(Vec<PathInfo>, (usize, usize)), String> {
    // Decompress the GBWT node records for the subgraph.
    let mut keys: Vec<usize> = Vec::new();
    let mut successors: BTreeMap<usize, Vec<(Pos, bool)>> = BTreeMap::new();
    for (handle, gbz_record) in subgraph.iter() {
        let gbwt_record = gbz_record.to_gbwt_record();
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
    let mut result: Vec<PathInfo> = Vec::new();
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
                len += subgraph.get(&pos.node).unwrap().sequence_len();
                curr = next_pos(pos, &successors);
            }
            if is_ref {
                if !support::encoded_path_is_canonical(&path) {
                    eprintln!("Warning: the reference path is not in canonical orientation");
                }
                result.push(PathInfo::new(path, len));
            } else if support::encoded_path_is_canonical(&path) {
                result.push(PathInfo::new(path, len));
            }
        }
    }

    let ref_id_offset = ref_id_offset.ok_or("Could not find the reference path".to_string())?;
    Ok((result, ref_id_offset))
}

fn distance_to(subgraph: &BTreeMap<usize, GBZRecord>, path: &[usize], path_offset: usize, node_offset: usize) -> usize {
    let mut result = node_offset;
    for handle in path.iter().take(path_offset) {
        result += subgraph.get(handle).unwrap().sequence_len();
    }
    result
}

// Returns all distinct paths and uses the weight field for storing their counts.
// Also updates `ref_id`.
fn make_distinct(
    paths: Vec<PathInfo>,
    ref_id: usize
) -> (Vec<PathInfo>, usize) {
    let ref_path = paths[ref_id].path.clone();
    let mut paths = paths;
    paths.sort_unstable();

    let mut new_paths: Vec<PathInfo> = Vec::new();
    let mut ref_id = 0;
    for info in paths.iter() {
        if new_paths.is_empty() || new_paths.last().unwrap().path != info.path {
            if info.path == ref_path {
                ref_id = new_paths.len();
            }
            new_paths.push(PathInfo::weighted(info.path.clone(), info.len));
        } else {
            new_paths.last_mut().unwrap().increment_weight();
        }
    }

    (new_paths, ref_id)
}

//-----------------------------------------------------------------------------

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

fn write_gfa_link<T: Write>(from: (&[u8], Orientation), to: (&[u8], Orientation), output: &mut T) -> io::Result<()> {
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

struct WalkMetadata {
    name: GBZPathName,
    end: usize,
    weight: Option<usize>,
}

impl WalkMetadata {
    fn path_interval(path: &GBZPath, interval: Range<usize>, weight: Option<usize>) -> Self {
        let mut name = path.name.clone();
        let end = name.fragment + interval.end;
        name.fragment += interval.start;
        WalkMetadata { name, end, weight }
    }

    fn anonymous(haplotype: usize, contig: &str, len: usize, weight: Option<usize>) -> Self {
        let path_name = GBZPathName::haplotype("unknown", contig, haplotype, 0);
        WalkMetadata {
            name: path_name,
            end: len,
            weight,
        }
    }
}

fn write_gfa_walk<T: Write>(path: &[usize], metadata: &WalkMetadata, output: &mut T) -> io::Result<()> {
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
    buffer.push(b'\n');

    output.write_all(&buffer)?;
    Ok(())
}

//-----------------------------------------------------------------------------

// TODO: This JSON output should also be in the library.

#[derive(Clone, Debug, PartialEq, Eq)]
enum JSONValue {
    Boolean(bool),
    String(String),
    Number(usize),
    Array(Vec<JSONValue>),
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

fn json_path(path: &[usize], metadata: &WalkMetadata) -> JSONValue {
    JSONValue::Object(vec![
        ("name".to_string(), JSONValue::String(metadata.name.path_fragment_name(metadata.end))),
        ("weight".to_string(), JSONValue::Number(metadata.weight.unwrap_or(1))),
        ("path".to_string(), JSONValue::Array(path.iter().map(
            |x| JSONValue::Object(vec![
                ("id".to_string(), JSONValue::String(support::node_id(*x).to_string())),
                ("is_reverse".to_string(), JSONValue::Boolean(support::node_orientation(*x) == Orientation::Reverse)),
            ])
        ).collect())),
    ])
}

//-----------------------------------------------------------------------------
