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
    let subgraph = if use_gbz {
        let graph: GBZ = serialize::load_from(&config.filename).map_err(|x| x.to_string())?;
        let path_index = PathIndex::new(&graph, GBZBase::INDEX_INTERVAL, false)?;
        Subgraph::from_gbz(&graph, Some(&path_index), &config.query)?
    } else {
        let database = GBZBase::open(&config.filename)?;
        let mut graph = GraphInterface::new(&database)?;
        Subgraph::from_db(&mut graph, &config.query)?
    };

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

        // FIXME node-based queries
        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("s", "sample", "sample name (default: no sample name)", "STR");
        opts.optopt("c", "contig", "contig name (required)", "STR");
        opts.optopt("o", "offset", "sequence offset (-o or -i is required)", "INT");
        opts.optopt("i", "interval", "sequence interval (-o or -i is required)", "INT..INT");
        let context_desc = format!("context length in bp (default: {})", Self::DEFAULT_CONTEXT);
        opts.optopt("n", "context", &context_desc, "INT");
        opts.optflag("d", "distinct", "output distinct haplotypes with weights");
        opts.optflag("C", "cigar", "output CIGAR strings for the haplotypes");
        opts.optflag("r", "reference-only", "output the reference but no other haplotypes");
        opts.optopt("f", "format", "output format (gfa or json, default: gfa)", "STR");
        let matches = opts.parse(&args[1..]).map_err(|x| x.to_string())?;

        let filename = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            let header = format!("Usage: {} [options] graph.gbz[.db]", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };

        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz[.db]", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }

        let sample = matches.opt_str("s").unwrap_or(String::from(REF_SAMPLE));
        let contig = matches.opt_str("c").ok_or(String::from("Contig name must be provided with --contig"))?;
        let path_name = FullPathName::reference(&sample, &contig);
        let (offset, context) = Self::parse_offset_and_context(&matches)?;

        let mut output = HaplotypeOutput::All;
        if matches.opt_present("d") {
            output = HaplotypeOutput::Distinct;
        }
        if matches.opt_present("r") {
            output = HaplotypeOutput::ReferenceOnly;
        }
        // FIXME: Use path_offset(), path_interval(), or node() depending on arguments.
        let query = SubgraphQuery::path_offset(&path_name, offset, context, output);
        let cigar = matches.opt_present("C");

        let mut format: OutputFormat = OutputFormat::Gfa;
        if let Some(s) = matches.opt_str("f") {
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
