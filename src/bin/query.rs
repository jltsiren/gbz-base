use gbz_base::{GBZBase, GraphInterface, PathIndex, Subgraph, SubgraphQuery, HaplotypeOutput};

use gbwt::{FullPathName, GBZ, REF_SAMPLE};

use simple_sds::serialize;

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
        subgraph.from_gbz(&graph, Some(&path_index), &config.query)?;
    } else {
        let database = GBZBase::open(&config.filename)?;
        let mut graph = GraphInterface::new(&database)?;
        subgraph.from_db(&mut graph, &config.query)?;
    }

    // Write the output.
    let mut output = io::stdout();
    match config.format {
        OutputFormat::Gfa => subgraph.write_gfa(&mut output, config.cigar).map_err(|x| x.to_string())?,
        OutputFormat::Json => subgraph.write_json(&mut output, config.cigar).map_err(|x| x.to_string())?,
    }

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
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
    cigar: bool,
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
        opts.optopt("", "sample", "sample name (default: no sample name)", "STR");
        opts.optopt("", "contig", "contig name (required for -o and -i)", "STR");
        opts.optopt("o", "offset", "sequence offset", "INT");
        opts.optopt("i", "interval", "sequence interval", "INT..INT");
        opts.optmulti("n", "node", "node identifier (may repeat)", "INT");
        let context_desc = format!("context length in bp (default: {})", Self::DEFAULT_CONTEXT);
        opts.optopt("", "context", &context_desc, "INT");
        opts.optflag("", "distinct", "output distinct haplotypes with weights");
        opts.optflag("", "reference-only", "output the reference but no other haplotypes");
        opts.optflag("", "cigar", "output CIGAR strings for the haplotypes");
        opts.optopt("", "format", "output format (gfa or json, default: gfa)", "STR");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let header = format!("Usage: {} [options] graph.gbz[.db]\n\nQuery type must be speficied using one of -o, -i, and -n.", program);
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
        let cigar = matches.opt_present("cigar");
        let mut format: OutputFormat = OutputFormat::Gfa;
        if let Some(s) = matches.opt_str("format") {
            match s.to_lowercase().as_str() {
                "gfa" => format = OutputFormat::Gfa,
                "json" => format = OutputFormat::Json,
                _ => return Err(format!("Invalid output format: {}", s)),
            }
        }

        Ok(Config { filename, query, cigar, format, })
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

    fn parse_integer(s: &str, option: &str) -> Result<usize, String> {
        s.parse::<usize>().map_err(|x| format!("--{}: {}", option, x))
    }

    fn parse_query(matches: &getopts::Matches) -> Result<SubgraphQuery, String> {
        let mut count = 0;
        if matches.opt_present("offset") { count += 1; }
        if matches.opt_present("interval") { count += 1; }
        if matches.opt_present("node") { count += 1; }
        if count != 1 {
            return Err("Exactly one of --offset, --interval, and --node must be provided".to_string());
        }

        let path_name = if matches.opt_present("node") {
            None
        } else {
            let sample = matches.opt_str("sample").unwrap_or(String::from(REF_SAMPLE));
            let contig = matches.opt_str("contig").ok_or(String::from("Contig name must be provided with --contig"))?;
            Some(FullPathName::reference(&sample, &contig))
        };
        let context = if let Some(s) = matches.opt_str("context") {
            Self::parse_integer(&s, "context")?
        } else {
            Self::DEFAULT_CONTEXT
        };
        let mut output = HaplotypeOutput::All;
        if matches.opt_present("distinct") {
            output = HaplotypeOutput::Distinct;
        }
        if matches.opt_present("reference-only") {
            output = HaplotypeOutput::ReferenceOnly;
        }

        if let Some(s) = matches.opt_str("offset") {
            let offset = Self::parse_integer(&s, "offset")?;
            let query = SubgraphQuery::path_offset(&path_name.unwrap(), offset, context, output);
            Ok(query)
        }
        else if let Some(s) = matches.opt_str("interval") {
            let interval = Self::parse_interval(&s)?;
            let query = SubgraphQuery::path_interval(&path_name.unwrap(), interval, context, output);
            Ok(query)
        } else {
            let node_strings = matches.opt_strs("node");
            let mut nodes = Vec::with_capacity(node_strings.len());
            for node in node_strings {
                let id = Self::parse_integer(&node, "node")?;
                nodes.push(id);
            }
            let query = SubgraphQuery::nodes(nodes, context, output);
            Ok(query)
        }
    }
}

//-----------------------------------------------------------------------------
