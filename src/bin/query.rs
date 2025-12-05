use gbz_base::{GBZBase, GraphInterface, GraphReference, PathIndex, Chains};
use gbz_base::{Subgraph, SubgraphQuery, HaplotypeOutput};
use gbz_base::{GAFBase, ReadSet, AlignmentOutput};

use gbwt::{FullPathName, Orientation, GBZ, REF_SAMPLE};
use gbwt::support;

use simple_sds::serialize;

use std::fs::OpenOptions;
use std::ops::Range;
use std::time::Instant;
use std::{env, io, process};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = Config::new()?;

    // Determine the type of the input file and extract the subgraph accordingly.
    let use_gbz = GBZ::is_gbz(&config.filename);
    let mut subgraph = Subgraph::new();
    if use_gbz {
        let graph: GBZ = serialize::load_from(&config.filename).map_err(|x| x.to_string())?;
        let path_index = PathIndex::new(&graph, GBZBase::INDEX_INTERVAL, false)?;
        let chains = match &config.chains {
            Some(file) => Some(Chains::load_from(file.as_ref())?),
            None => None,
        };
        subgraph.from_gbz(&graph, Some(&path_index), chains.as_ref(), &config.query)?;
        subgraph_statistics(&subgraph);
        write_subgraph(&subgraph, &config)?;
        extract_gaf(GraphReference::Gbz(&graph), &subgraph, &config)?;
    } else {
        let database = GBZBase::open(&config.filename)?;
        let mut graph = GraphInterface::new(&database)?;
        subgraph.from_db(&mut graph, &config.query)?;
        subgraph_statistics(&subgraph);
        write_subgraph(&subgraph, &config)?;
        extract_gaf(GraphReference::Db(&mut graph), &subgraph, &config)?;
    }

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
}

//-----------------------------------------------------------------------------

fn subgraph_statistics(subgraph: &Subgraph) {
    eprintln!("Subgraph contains {} nodes and {} paths", subgraph.nodes(), subgraph.paths());
}

fn write_subgraph(subgraph: &Subgraph, config: &Config) -> Result<(), String> {
    let mut output = io::stdout().lock();
    match config.format {
        OutputFormat::Gfa => subgraph.write_gfa(&mut output, config.cigar).map_err(|x| x.to_string()),
        OutputFormat::Json => subgraph.write_json(&mut output, config.cigar).map_err(|x| x.to_string()),
    }
}

fn extract_gaf(graph: GraphReference<'_, '_>, subgraph: &Subgraph, config: &Config) -> Result<(), String> {
    if !config.write_gaf() {
        return Ok(());
    }

    let gaf_base_file = config.gaf_base.as_ref().unwrap();
    let gaf_base = GAFBase::open(gaf_base_file)?;
    let read_set = ReadSet::new(graph, subgraph, &gaf_base, config.alignment_output)?;
    if config.alignment_output == AlignmentOutput::Clipped {
        eprintln!(
            "Extracted {} fragments for {} reads in {} alignment blocks with {} node records in {} clusters",
            read_set.len(), read_set.unclipped(), read_set.blocks(), read_set.node_records(), read_set.clusters()
        );
    } else {
        eprintln!(
            "Extracted {} reads in {} alignment blocks with {} node records in {} clusters",
            read_set.len(), read_set.blocks(), read_set.node_records(), read_set.clusters()
        );
    }

    let gaf_output_file = config.gaf_output.as_ref().unwrap();
    let mut options = OpenOptions::new();
    options.write(true).create(true).truncate(true);
    let mut gaf_output = options.open(gaf_output_file).map_err(|x| x.to_string())?;
    read_set.to_gaf(&mut gaf_output).map_err(
        |x| format!("Failed to write GAF to {}: {}", gaf_output_file, x)
    )
}

//-----------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum OutputFormat {
    Gfa,
    Json,
}

struct Config {
    filename: String,
    query: SubgraphQuery,
    chains: Option<String>,
    cigar: bool,
    format: OutputFormat,
    gaf_base: Option<String>,
    gaf_output: Option<String>,
    alignment_output: AlignmentOutput,
}

impl Config {
    // Default context length in bp.
    const DEFAULT_CONTEXT: usize = 100;

    fn new() -> Result<Config, String> {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("", "sample", "sample name (default: no sample name)", "STR");
        opts.optopt("", "contig", "contig name (required for -o and -i)", "STR");
        opts.optopt("o", "offset", "sequence offset", "INT");
        opts.optopt("i", "interval", "half-open sequence interval", "INT..INT");
        opts.optmulti("n", "node", "node identifier (may repeat)", "INT");
        opts.optopt("b", "between", "subgraph between boundary nodes", "INT[+-]:INT[+-]");
        opts.optopt("", "limit", "safety limit for the number of nodes in -b", "INT");
        let context_desc = format!("context length in bp (not for -b; default: {})", Self::DEFAULT_CONTEXT);
        opts.optopt("", "context", &context_desc, "INT");
        opts.optflag("", "snarls", "include nodes in covered top-level snarls");
        opts.optopt("", "chains", "top-level chains file (for --snarls with a GBZ graph)", "FILE");
        opts.optflag("", "distinct", "output distinct haplotypes with weights");
        opts.optflag("", "reference-only", "output the reference but no other haplotypes");
        opts.optflag("", "cigar", "output CIGAR strings for the haplotypes");
        opts.optopt("", "format", "output format (gfa or json; default: gfa)", "STR");
        opts.optopt("", "gaf-base", "GAF-base file (for GAF output)", "FILE");
        opts.optopt("", "gaf-output", "GAF output file (for GAF output)", "FILE");
        opts.optopt("", "alignments", "alignment selection (overlapping, clipped, or contained; default: clipped)", "STR");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let header = format!("Usage: {} [options] graph.gbz[.db]\n\nQuery type must be speficied using one of -o, -i, -n, and -b.", program);
        if matches.opt_present("help") {
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }

        let filename = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };

        let query = Self::parse_query(&matches)?;
        let chains = matches.opt_str("chains");
        let cigar = matches.opt_present("cigar");
        let mut format: OutputFormat = OutputFormat::Gfa;
        if let Some(s) = matches.opt_str("format") {
            match s.to_lowercase().as_str() {
                "gfa" => format = OutputFormat::Gfa,
                "json" => format = OutputFormat::Json,
                _ => return Err(format!("Invalid output format: {}", s)),
            }
        }

        let gaf_base = matches.opt_str("gaf-base");
        let gaf_output = matches.opt_str("gaf-output");

        let alignment_output = if let Some(s) = matches.opt_str("alignments") {
            match s.to_lowercase().as_str() {
                "overlapping" => AlignmentOutput::Overlapping,
                "clipped" => AlignmentOutput::Clipped,
                "contained" => AlignmentOutput::Contained,
                _ => return Err(format!("Invalid alignment selection: {}", s)),
            }
        } else {
            AlignmentOutput::Clipped
        };

        Ok(Config { filename, query, chains, cigar, format, gaf_base, gaf_output, alignment_output })
    }

    fn write_gaf(&self) -> bool {
        self.gaf_base.is_some() && self.gaf_output.is_some()
    }

    fn parse_interval(s: &str) -> Result<Range<usize>, String> {
        let mut parts = s.split("..");
        let start = parts.next().ok_or(format!("Invalid interval: {}", s))?;
        let start = start.parse::<usize>().map_err(|x| format!("Failed to parse interval start: {}", x))?;
        let end = parts.next().ok_or(format!("Invalid interval: {}", s))?;
        let end = end.parse::<usize>().map_err(|x| format!("Failed to parse interval end: {}", x))?;
        if parts.next().is_some() {
            return Err(format!("Invalid interval: {}", s));
        }
        Ok(start..end)
    }

    // Parses a node id that may be followed by a + or a -.
    fn parse_handle(s: &str) -> Result<usize, String> {
        let mut len = s.len();
        let orientation = if s.ends_with('+') {
            len -= 1;
            Orientation::Forward
        } else if s.ends_with('-') {
            len -= 1;
            Orientation::Reverse
        } else {
            Orientation::Forward
        };
        let id = s[..len].parse::<usize>().map_err(|x| format!("Failed to parse (oriented) node: {}", x))?;
        Ok(support::encode_node(id, orientation))
    }

    fn parse_between(s: &str) -> Result<(usize, usize), String> {
        let mut parts = s.split(':');
        let start = parts.next().ok_or(format!("Invalid pair of (oriented) nodes: {}", s))?;
        let start = Self::parse_handle(start)?;
        let end = parts.next().ok_or(format!("Invalid pair of (oriented) nodes: {}", s))?;
        let end = Self::parse_handle(end)?;
        if parts.next().is_some() {
            return Err(format!("Invalid pair of (oriented) nodes: {}", s));
        }
        Ok((start, end))
    }

    fn parse_integer(s: &str, option: &str) -> Result<usize, String> {
        s.parse::<usize>().map_err(|x| format!("Failed to parse --{}: {}", option, x))
    }

    fn parse_query(matches: &getopts::Matches) -> Result<SubgraphQuery, String> {
        let mut count = 0;
        let mut needs_path_name = false;
        if matches.opt_present("offset") { count += 1; needs_path_name = true; }
        if matches.opt_present("interval") { count += 1; needs_path_name = true; }
        if matches.opt_present("node") { count += 1; }
        if matches.opt_present("between") { count += 1; }
        if count != 1 {
            return Err("Exactly one of --offset, --interval, --node, and --between must be provided".to_string());
        }

        let path_name = if needs_path_name {
            let sample = matches.opt_str("sample").unwrap_or(String::from(REF_SAMPLE));
            let contig = matches.opt_str("contig").ok_or(String::from("Contig name must be provided with --contig"))?;
            Some(FullPathName::reference(&sample, &contig))
        } else {
            None
        };
        let context = if let Some(s) = matches.opt_str("context") {
            Self::parse_integer(&s, "context")?
        } else {
            Self::DEFAULT_CONTEXT
        };
        let limit = if let Some(s) = matches.opt_str("limit") {
            Some(Self::parse_integer(&s, "limit")?)
        } else {
            None
        };
        let snarls = matches.opt_present("snarls");
        let mut output = HaplotypeOutput::All;
        if matches.opt_present("distinct") {
            output = HaplotypeOutput::Distinct;
        }
        if matches.opt_present("reference-only") {
            output = HaplotypeOutput::ReferenceOnly;
        }

        let query = if let Some(s) = matches.opt_str("offset") {
            let offset = Self::parse_integer(&s, "offset")?;
            SubgraphQuery::path_offset(&path_name.unwrap(), offset)
        }
        else if let Some(s) = matches.opt_str("interval") {
            let interval = Self::parse_interval(&s)?;
            SubgraphQuery::path_interval(&path_name.unwrap(), interval)
        } else if let Some(s) = matches.opt_str("between") {
            let (start, end) = Self::parse_between(&s)?;
            SubgraphQuery::between(start, end, limit)
        } else {
            let node_strings = matches.opt_strs("node");
            let mut nodes = Vec::with_capacity(node_strings.len());
            for node in node_strings {
                let id = Self::parse_integer(&node, "node")?;
                nodes.push(id);
            }
            SubgraphQuery::nodes(nodes)
        };

        Ok(query.with_context(context).with_snarls(snarls).with_output(output))
    }
}

//-----------------------------------------------------------------------------
