use gbz_base::{GBZBase, GraphInterface, GraphReference, PathIndex, Subgraph, SubgraphQuery, HaplotypeOutput, Chains};
use gbz_base::{GAFBase, ReadSet};

use gbwt::{FullPathName, GBZ, REF_SAMPLE};

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
            Some(file) => Some(Chains::from_file(file.as_ref())?),
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
    let mut output = io::stdout();
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
    let read_set = ReadSet::new(graph, subgraph, &gaf_base, config.contained)?;
    eprintln!("Extracted {} reads in {} alignment blocks with {} node records", read_set.len(), read_set.blocks(), read_set.node_records());

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
    contained: bool,
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
        opts.optflag("", "snarls", "include nodes in covered top-level snarls");
        opts.optopt("", "chains", "top-level chains file (for --snarls with a GBZ graph)", "FILE");
        opts.optflag("", "distinct", "output distinct haplotypes with weights");
        opts.optflag("", "reference-only", "output the reference but no other haplotypes");
        opts.optflag("", "cigar", "output CIGAR strings for the haplotypes");
        opts.optopt("", "format", "output format (gfa or json, default: gfa)", "STR");
        opts.optopt("", "gaf-base", "GAF-base file (for GAF output)", "FILE");
        opts.optopt("", "gaf-output", "GAF output file (for GAF output)", "FILE");
        opts.optflag("", "contained", "output only reads that are fully within the subgraph");
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
        let contained = matches.opt_present("contained");

        Ok(Config { filename, query, chains, cigar, format, gaf_base, gaf_output, contained })
    }

    fn write_gaf(&self) -> bool {
        self.gaf_base.is_some() && self.gaf_output.is_some()
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
