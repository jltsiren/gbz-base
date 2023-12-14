use gbz_base::{GBZBase, GraphInterface};
use gbz_base::gfa::{WalkMetadata, write_gfa, write_gfa_walk};
use gbz_base::query::{query_position, extract_context, extract_paths, distinct_paths, distance_to};

use std::ops::Range;
use std::time::Instant;
use std::{env, io, process};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = Config::new()?;

    // Open the database.
    let database = GBZBase::open(&config.filename)?;
    let mut graph = GraphInterface::new(&database)?;

    // Find the reference path.
    let ref_path = graph.find_path(
        &config.sample, &config.contig, config.haplotype, config.fragment
    )?;
    let ref_info = ref_path.ok_or(format!("Cannot find path {}", config.path_name()))?;
    if !ref_info.is_indexed {
        return Err(format!("Path {} has not been indexed for random access", config.path_name()));
    }

    // Find the query position on the reference path.
    let (query_pos, gbwt_pos) = query_position(
        &mut graph,
        &ref_info,
        config.offset
    )?;

    // Extract all GBWT records within the context.
    let subgraph = extract_context(&mut graph, query_pos, config.context)?;

    // Extract paths.
    let (mut paths, (mut ref_id, path_offset)) = extract_paths(&subgraph, gbwt_pos)?;
    let ref_start = config.offset - distance_to(&subgraph, &paths[ref_id].path, path_offset, query_pos.offset);
    let ref_interval = ref_start..ref_start + paths[ref_id].len;
    if config.output == HaplotypeOutput::Distinct {
        // TODO: We should use bidirectional search in GBWT to find the distinct paths directly.
        (paths, ref_id) = distinct_paths(paths, ref_id);
    }

    // GFA output: segments and links.
    let mut output = io::stdout();
    let reference_samples = Some(ref_info.sample.clone());
    write_gfa(&subgraph, reference_samples, &mut output).map_err(|x| x.to_string())?;

    // GFA output: walks. 
    let ref_metadata = WalkMetadata::from_gbz_path(&ref_info, ref_interval, paths[ref_id].weight);
    write_gfa_walk(&paths[ref_id].path, &ref_metadata, &mut output).map_err(|x| x.to_string())?;
    if config.output != HaplotypeOutput::ReferenceOnly {
        let mut haplotype = 1;
        for (id, path_info) in paths.iter().enumerate() {
            if id == ref_id {
                continue;
            }
            let metadata = WalkMetadata::anonymous(haplotype, &ref_info.contig, path_info.len, path_info.weight);
            write_gfa_walk(&path_info.path, &metadata, &mut output).map_err(|x| x.to_string())?;
            haplotype += 1;
        }
    }

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
}

//-----------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HaplotypeOutput {
    All,
    Distinct,
    ReferenceOnly,
}

// TODO: haplotype, fragment
pub struct Config {
    pub filename: String,
    pub sample: String,
    pub contig: String,
    pub haplotype: usize,
    pub fragment: usize,
    pub offset: usize,
    pub context: usize,
    pub output: HaplotypeOutput,
}

impl Config {
    // Default context length in bp.
    pub const DEFAULT_CONTEXT: usize = 100;

    pub fn new() -> Result<Config, String> {
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
        let haplotype: usize = 0;
        let fragment: usize = 0;
        let (offset, context) = Self::parse_offset_and_context(&matches)?;

        let mut output = HaplotypeOutput::All;
        if matches.opt_present("d") {
            output = HaplotypeOutput::Distinct;
        }
        if matches.opt_present("r") {
            output = HaplotypeOutput::ReferenceOnly;
        }

        Ok(Config { filename, sample, contig, haplotype, fragment, offset, context, output, })
    }

    pub fn path_name(&self) -> String {
        format!("{}#{}#{}@{}", self.sample, self.haplotype, self.contig, self.fragment)
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

