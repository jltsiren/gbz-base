use gbz_base::{GBZBase, GraphInterface, GraphReference, PathIndex};
use gbz_base::{Subgraph, SubgraphQuery, HaplotypeOutput, SnarlOutput};
use gbz_base::{GAFBase, ReadSet, AlignmentOutput};
use gbz_base::{formats, utils};
use gbz_base::{Error, Result};

use gbz::{FullPathName, Orientation, GBZ, GENERIC_SAMPLE};
use gbz::support;

use simple_sds::{binaries, serialize};

use std::fs::{self, OpenOptions};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::io;

use clap::{Args, Parser, Subcommand, ValueEnum};

//-----------------------------------------------------------------------------

#[derive(Parser)]
#[command(
    name = "gbz-base",
    version,
    about = "Pangenome graphs in SQLite",
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
enum Commands {
    /// Build a GBZ-base from a GBZ graph
    Construct(ConstructArgs),
    /// Extract a subgraph from a GBZ-base or GBZ graph
    Query(QueryArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Construct(args) => construct(args),
        Commands::Query(args) => query(args),
    }
}

// Parses an unsigned integer that may use suffixes such as `k` or `M`.
fn parse_quantity(s: &str) -> Result<usize> {
    binaries::parse_unsigned(s).map_err(Error::invalid_query)
}

//-----------------------------------------------------------------------------

#[derive(Args)]
struct ConstructArgs {
    /// GBZ graph file
    #[arg(value_name = "graph.gbz")]
    graph: PathBuf,

    /// Top-level chain file (optional)
    #[arg(short, long, value_name = "FILE")]
    chains: Option<PathBuf>,

    /// Output file name (default: <input>.db)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Overwrite the database file if it exists
    #[arg(long)]
    overwrite: bool,
}

fn construct(args: ConstructArgs) -> Result<()> {
    let start_time = Instant::now();

    let db_file = args.output.unwrap_or_else(|| {
        let mut name = args.graph.clone();
        name.add_extension("db");
        name
    });

    // Check if the database already exists.
    if binaries::file_exists(&db_file) {
        if args.overwrite {
            eprintln!("Overwriting database {}", db_file.display());
            fs::remove_file(&db_file)?;
        } else {
            return Err(Error::invalid_query(format!("Database {} already exists", db_file.display())));
        }
    }

    // Create the database.
    let chains_file: Option<&Path> = args.chains.as_deref();
    GBZBase::create_from_files(args.graph.as_ref(), chains_file, &db_file)?;

    // Statistics.
    let database = GBZBase::open(&db_file)?;
    eprintln!("The graph contains {} nodes in {} chains with {} links",
        database.nodes(), database.chains(), database.chain_links()
    );
    eprintln!("There are {} paths representing {} samples, {} haplotypes, and {} contigs",
        database.paths(), database.samples(), database.haplotypes(), database.contigs()
    );
    let size = database.file_size().unwrap_or_else(|| String::from("unknown"));
    eprintln!("Final database size: {}", size);

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);
    utils::report_peak_memory_usage();

    Ok(())
}

//-----------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum HaplotypeSelection {
    All,
    Distinct,
    ReferenceOnly,
    None,
}

impl From<HaplotypeSelection> for HaplotypeOutput {
    fn from(value: HaplotypeSelection) -> Self {
        match value {
            HaplotypeSelection::All => HaplotypeOutput::All,
            HaplotypeSelection::Distinct => HaplotypeOutput::Distinct,
            HaplotypeSelection::ReferenceOnly => HaplotypeOutput::ReferenceOnly,
            HaplotypeSelection::None => HaplotypeOutput::None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum OutputFormat {
    Gfa,
    Json,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum AlignmentSelection {
    Overlapping,
    Clipped,
    Contained,
}

impl From<AlignmentSelection> for AlignmentOutput {
    fn from(value: AlignmentSelection) -> Self {
        match value {
            AlignmentSelection::Overlapping => AlignmentOutput::Overlapping,
            AlignmentSelection::Clipped => AlignmentOutput::Clipped,
            AlignmentSelection::Contained => AlignmentOutput::Contained,
        }
    }
}

// Default context length in bp.
const DEFAULT_CONTEXT: usize = 100;

#[derive(Args)]
struct QueryArgs {
    /// GBZ graph or GBZ-base file
    #[arg(value_name = "graph.gbz[.db]")]
    filename: String,

    /// Sample name (default: no sample name)
    #[arg(long, value_name = "STR")]
    sample: Option<String>,

    /// Contig name (required for --offset and --interval)
    #[arg(long, value_name = "STR")]
    contig: Option<String>,

    /// Sequence offset
    #[arg(short, long, value_name = "INT")]
    offset: Option<usize>,

    /// Half-open sequence interval
    #[arg(short, long, value_name = "INT..INT")]
    interval: Option<String>,

    /// Node identifier (may repeat)
    #[arg(short, long, value_name = "INT")]
    node: Vec<usize>,

    /// Node corresponding to a handle (may repeat)
    #[arg(long, value_name = "INT")]
    handle: Vec<usize>,

    /// Subgraph between boundary nodes
    #[arg(short, long, value_name = "INT[+-]:INT[+-]")]
    between: Option<String>,

    /// Safety limit for the number of nodes in the subgraph
    #[arg(long, value_name = "INT", value_parser = parse_quantity)]
    limit: Option<usize>,

    /// Context length in bp (not for --between)
    #[arg(long, value_name = "INT", value_parser = parse_quantity, default_value_t = DEFAULT_CONTEXT)]
    context: usize,

    /// Extend subgraph to include contained top-level snarls
    #[arg(long)]
    snarls: bool,

    /// Extend subgraph to include overlapping top-level snarls
    #[arg(long)]
    extend_snarls: bool,

    /// Top-level chains file (for --snarls/--extend-snarls with a GBZ graph)
    #[arg(long, value_name = "FILE")]
    chains: Option<String>,

    /// Haplotype output selection
    #[arg(long, value_name = "SELECTION", value_enum, default_value_t = HaplotypeSelection::All)]
    haplotypes: HaplotypeSelection,

    /// Output CIGAR strings for the haplotypes
    #[arg(long)]
    cigar: bool,

    /// Output format
    #[arg(long, value_name = "FORMAT", value_enum, default_value_t = OutputFormat::Gfa)]
    format: OutputFormat,

    /// GAF-base file (for GAF output)
    #[arg(long, value_name = "FILE")]
    gaf_base: Option<String>,

    /// GAF output file (for GAF output)
    #[arg(long, value_name = "FILE")]
    gaf_output: Option<String>,

    /// Write GAF output to stdout (no subgraph output)
    #[arg(long)]
    gaf_only: bool,

    /// Alignment selection (for GAF output)
    #[arg(long, value_name = "SELECTION", value_enum, default_value_t = AlignmentSelection::Clipped)]
    alignments: AlignmentSelection,
}

struct QueryConfig {
    filename: String,
    query: SubgraphQuery,
    chains: Option<String>,
    cigar: bool,
    format: OutputFormat,
    gaf_base: Option<String>,
    gaf_output: Option<String>,
    gaf_only: bool,
    alignment_output: AlignmentOutput,
}

impl QueryConfig {
    fn write_gaf(&self) -> bool {
        self.gaf_base.is_some()
    }
}

fn query(args: QueryArgs) -> Result<()> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = build_query_config(args)?;

    // Determine the type of the input file and extract the subgraph accordingly.
    let use_gbz = GBZ::is_gbz(&config.filename);
    let mut subgraph = Subgraph::new();
    if use_gbz {
        let graph: GBZ = serialize::load_from(&config.filename)?;
        let path_index = PathIndex::new(&graph, GBZBase::INDEX_INTERVAL, false)?;
        let chains = match &config.chains {
            Some(file) => Some(serialize::load_from(file)?),
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
    utils::report_peak_memory_usage();

    Ok(())
}

fn build_query_config(args: QueryArgs) -> Result<QueryConfig> {
    validate_gaf_output(&args)?;
    let query = build_subgraph_query(&args)?;
    Ok(QueryConfig {
        filename: args.filename,
        query,
        chains: args.chains,
        cigar: args.cigar,
        format: args.format,
        gaf_base: args.gaf_base,
        gaf_output: args.gaf_output,
        gaf_only: args.gaf_only,
        alignment_output: args.alignments.into(),
    })
}

fn validate_gaf_output(args: &QueryArgs) -> Result<()> {
    if args.gaf_base.is_none() && (args.gaf_output.is_some() || args.gaf_only) {
        return Err(Error::invalid_query("GAF output requires --gaf-base"));
    }
    if args.gaf_base.is_some() && !(args.gaf_output.is_some() || args.gaf_only) {
        return Err(Error::invalid_query("Option --gaf-base requires --gaf-output or --gaf-only"));
    }
    if args.gaf_only && args.gaf_output.is_some() {
        return Err(Error::invalid_query("Options --gaf-only and --gaf-output are mutually exclusive"));
    }
    Ok(())
}

fn build_subgraph_query(args: &QueryArgs) -> Result<SubgraphQuery> {
    let mut count = 0;
    let mut needs_path_name = false;
    if args.offset.is_some() { count += 1; needs_path_name = true; }
    if args.interval.is_some() { count += 1; needs_path_name = true; }
    if !args.node.is_empty() || !args.handle.is_empty() { count += 1; }
    if args.between.is_some() { count += 1; }
    if count != 1 {
        return Err(Error::invalid_query("Exactly one of --offset, --interval, --node (or --handle), and --between must be provided"));
    }

    let path_name = if needs_path_name {
        let sample = args.sample.clone().unwrap_or_else(|| String::from(GENERIC_SAMPLE));
        let contig = args.contig.clone().ok_or_else(|| Error::invalid_query("Contig name must be provided with --contig"))?;
        Some(FullPathName::reference(&sample, &contig))
    } else {
        None
    };
    let snarls = match (args.snarls, args.extend_snarls) {
        (true, true) => SnarlOutput::Overlapping,
        (true, false) => SnarlOutput::Contained,
        (false, true) => SnarlOutput::Overlapping,
        (false, false) => SnarlOutput::None,
    };
    let output = args.haplotypes.into();

    let query = if let Some(offset) = args.offset {
        SubgraphQuery::path_offset(&path_name.unwrap(), offset)
    } else if let Some(s) = &args.interval {
        let interval = parse_interval(s)?;
        SubgraphQuery::path_interval(&path_name.unwrap(), interval)
    } else if let Some(s) = &args.between {
        let (start, end) = parse_between(s)?;
        SubgraphQuery::between(start, end)
    } else {
        let mut nodes = Vec::with_capacity(args.node.len() + args.handle.len());
        for id in &args.node {
            nodes.push(*id);
        }
        for handle in &args.handle {
            nodes.push(support::node_id(*handle));
        }
        SubgraphQuery::nodes(nodes)
    };

    Ok(query.with_limit(args.limit).with_context(args.context).with_snarls(snarls).with_haplotypes(output))
}

fn parse_interval(s: &str) -> Result<Range<usize>> {
    let mut parts = s.split("..");
    let start = parts.next().ok_or_else(|| Error::invalid_query(format!("Invalid interval: {}", s)))?;
    let start = start.parse::<usize>().map_err(|x| Error::invalid_query(format!("Failed to parse interval start: {}", x)))?;
    let end = parts.next().ok_or_else(|| Error::invalid_query(format!("Invalid interval: {}", s)))?;
    let end = end.parse::<usize>().map_err(|x| Error::invalid_query(format!("Failed to parse interval end: {}", x)))?;
    if parts.next().is_some() {
        return Err(Error::invalid_query(format!("Invalid interval: {}", s)));
    }
    Ok(start..end)
}

// Parses a node id that may be followed by a + or a -.
fn parse_handle(s: &str) -> Result<usize> {
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
    let id = s[..len].parse::<usize>().map_err(|x| Error::invalid_query(format!("Failed to parse (oriented) node: {}", x)))?;
    Ok(support::encode_node(id, orientation))
}

fn parse_between(s: &str) -> Result<(usize, usize)> {
    let mut parts = s.split(':');
    let start = parts.next().ok_or_else(|| Error::invalid_query(format!("Invalid pair of (oriented) nodes: {}", s)))?;
    let start = parse_handle(start)?;
    let end = parts.next().ok_or_else(|| Error::invalid_query(format!("Invalid pair of (oriented) nodes: {}", s)))?;
    let end = parse_handle(end)?;
    if parts.next().is_some() {
        return Err(Error::invalid_query(format!("Invalid pair of (oriented) nodes: {}", s)));
    }
    Ok((start, end))
}

//-----------------------------------------------------------------------------

fn subgraph_statistics(subgraph: &Subgraph) {
    eprintln!("Subgraph contains {} nodes and {} paths", subgraph.nodes(), subgraph.paths());
}

fn write_subgraph(subgraph: &Subgraph, config: &QueryConfig) -> Result<()> {
    if config.gaf_only {
        return Ok(());
    }
    let mut output = io::stdout().lock();
    match config.format {
        OutputFormat::Gfa => subgraph.write_gfa(&mut output, config.cigar).map_err(Error::io),
        OutputFormat::Json => subgraph.write_json(&mut output, config.cigar).map_err(Error::io),
    }
}

fn extract_gaf(graph: GraphReference<'_, '_>, subgraph: &Subgraph, config: &QueryConfig) -> Result<()> {
    if !config.write_gaf() {
        return Ok(());
    }

    // Open the database and check that it is compatible with the graph.
    let gaf_base_file = config.gaf_base.as_ref().unwrap();
    let gaf_base = GAFBase::open(gaf_base_file)?;
    let mut graph = graph;
    let reference = graph.graph_name()?;
    let alignments = gaf_base.graph_name()?;
    utils::require_valid_reference(&alignments, &reference)?;

    // Extract the reads.
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

    let (gaf_output_file, mut gaf_output) = if config.gaf_only {
        (String::from("stdout"), Box::new(io::stdout().lock()) as Box<dyn io::Write>)
    } else {
        let gaf_output_file = config.gaf_output.as_ref().unwrap();
        let mut options = OpenOptions::new();
        options.write(true).create(true).truncate(true);
        let gaf_output = options.open(gaf_output_file)?;
        (gaf_output_file.clone(), Box::new(gaf_output) as Box<dyn io::Write>)
    };
    formats::write_gaf_file_header(&mut gaf_output, None).map_err(
        |x| Error::io(format!("Failed to write GAF header to {}: {}", gaf_output_file, x))
    )?;
    let graph_name = subgraph.graph_name();
    if let Some(name) = graph_name {
        let header_lines = name.to_gaf_header_lines();
        formats::write_header_lines(&header_lines, &mut gaf_output).map_err(
            |x| Error::io(format!("Failed to write GAF header lines to {}: {}", gaf_output_file, x))
        )?;
    }
    read_set.to_gaf(&mut gaf_output).map_err(
        |x| Error::io(format!("Failed to write GAF to {}: {}", gaf_output_file, x))
    )
}

//-----------------------------------------------------------------------------
