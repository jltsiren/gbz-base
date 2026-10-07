use std::fs::{self, File};
use std::io::{self, Write, BufWriter};
use std::path::PathBuf;
use std::sync::mpsc;
use std::time::Instant;
use std::thread;

use gbz_base::{GBZBase, GraphInterface};
use gbz_base::{GAFBase, GAFBaseParams, GraphReference};
use gbz_base::db::FileType;
use gbz_base::gaf_sort::{sort_gaf, KeyType, SortParameters};
use gbz_base::{db, formats, utils};
use gbz_base::ReadSet;
use gbz_base::{Error, Result};

use gbz::{GBZ, GraphName};

use simple_sds::{binaries, serialize};

use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum};
use clap::parser::ValueSource;

//-----------------------------------------------------------------------------

#[derive(Parser)]
#[command(
    name = "gaf-base",
    version,
    about = "Sequence alignments in SQLite",
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a GAF-base from sorted GAF alignments
    Construct(ConstructArgs),
    /// Convert a GAF-base back to GAF
    Decompress(DecompressArgs),
    /// Sort a GAF file for GAF-base construction
    Sort(SortArgs),
}

fn main() -> Result<()> {
    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).map_err(Error::invalid_query)?;
    match cli.command {
        Commands::Construct(args) => {
            // `subcommand_matches` is guaranteed to succeed, as we just parsed the `construct` subcommand.
            let sub_matches = matches.subcommand_matches("construct").unwrap();
            construct(args, sub_matches)
        },
        Commands::Decompress(args) => decompress(args),
        Commands::Sort(args) => {
            // `subcommand_matches` is guaranteed to succeed, as we just parsed the `sort` subcommand.
            let sub_matches = matches.subcommand_matches("sort").unwrap();
            sort(args, sub_matches)
        },
    }
}

// Parses an unsigned integer that may use suffixes such as `k` or `M`.
fn parse_quantity(s: &str) -> Result<usize> {
    binaries::parse_unsigned(s).map_err(Error::invalid_query)
}

// Parameter preset shared by the `construct` and `sort` subcommands.
// The variants match the presets in `GAFBaseParams` and `SortParameters`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Preset {
    Default,
    Short,
    Long,
}

impl Preset {
    // Returns the preset name understood by `GAFBaseParams::with_preset` and
    // `SortParameters::with_preset`.
    fn as_str(self) -> &'static str {
        match self {
            Preset::Default => "default",
            Preset::Short => "short",
            Preset::Long => "long",
        }
    }
}

//-----------------------------------------------------------------------------

#[derive(Args)]
struct ConstructArgs {
    /// GAF alignment file (may be gzip-compressed)
    #[arg(value_name = "alignments.gaf[.gz]")]
    gaf: PathBuf,

    /// GBWT file name
    #[arg(short, long, value_name = "FILE")]
    gbwt: Option<PathBuf>,

    /// Build a reference-free GAF-base using this graph
    #[arg(short = 'r', long = "ref-free", value_name = "FILE")]
    ref_free: Option<PathBuf>,

    /// Parameter preset
    #[arg(long, value_name = "PRESET", value_enum, default_value_t = Preset::Default)]
    preset: Preset,

    /// Number of alignments per block
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = GAFBaseParams::BLOCK_SIZE)]
    block_size: usize,

    /// Do not store quality strings
    #[arg(long)]
    no_quality: bool,

    /// Do not store unsupported optional fields
    #[arg(long)]
    no_optional: bool,

    /// Output file name (default: <input>.db)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Overwrite the database file if it exists
    #[arg(long)]
    overwrite: bool,
}

fn construct(args: ConstructArgs, matches: &clap::ArgMatches) -> Result<()> {
    let start_time = Instant::now();

    // Start from the preset and let each option override it only if the user passed
    // it explicitly. Otherwise clap's default values would silently override the preset.
    let mut params = GAFBaseParams::with_preset(args.preset.as_str())?;
    let explicit = |id: &str| matches.value_source(id) == Some(ValueSource::CommandLine);
    if explicit("block_size") {
        params.block_size = args.block_size;
    }
    if args.no_quality && explicit("no_quality") {
        params.store_quality_strings = false;
    }
    if args.no_optional && explicit("no_optional") {
        params.store_optional_fields = false;
    }
    if args.ref_free.is_some() && explicit("ref_free") {
        params.reference_free = true;
    }

    let db_file = args.output.unwrap_or_else(|| {
        let mut name = args.gaf.clone();
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
    if let Some(graph_file) = &args.ref_free {
        match db::identify_file(graph_file) {
            FileType::Gbz => {
                eprintln!("Loading GBZ graph {}", graph_file.display());
                let graph: GBZ = serialize::load_from(graph_file)?;
                GAFBase::create_from_files(
                    &args.gaf, args.gbwt.as_deref(), &db_file,
                    GraphReference::Gbz(&graph), &params
                )?;
            },
            FileType::Version(v) => {
                if v != GBZBase::VERSION {
                    return Err(Error::unsupported(format!("File {} is {}; expected {}", graph_file.display(), v, GBZBase::VERSION)));
                }
                eprintln!("Opening GBZ-base {}", graph_file.display());
                let database = GBZBase::open(graph_file)?;
                let mut graph = GraphInterface::new(&database)?;
                GAFBase::create_from_files(
                    &args.gaf, args.gbwt.as_deref(), &db_file,
                    GraphReference::Db(&mut graph), &params
                )?;
            },
            _ => {
                return Err(Error::unsupported(format!("File {} is not a valid graph", graph_file.display())));
            }
        };
    } else {
        GAFBase::create_from_files(&args.gaf, args.gbwt.as_deref(), &db_file, GraphReference::None, &params)?;
    }

    // Statistics.
    let database = GAFBase::open(&db_file)?;
    eprintln!(
        "The database contains {} nodes and {} alignments in {} blocks",
        database.nodes(), database.alignments(), database.blocks()
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

// Default number of blocks in a single ReadSet.
const DEFAULT_CHUNK_SIZE: usize = 100;

#[derive(Args)]
struct DecompressArgs {
    /// GAF-base file
    #[arg(value_name = "gaf_base.db")]
    gaf_base: PathBuf,

    /// Use this GBZ graph as the reference
    #[arg(short, long, value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Chunk size in blocks
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = DEFAULT_CHUNK_SIZE)]
    chunk_size: usize,

    /// Output file name (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
}

fn decompress(args: DecompressArgs) -> Result<()> {
    let start_time = Instant::now();

    if args.chunk_size < 1 {
        return Err(Error::invalid_query("--chunk-size must be positive"));
    }

    // Inputs. Increasing cache size helps a bit.
    let database = GAFBase::open(&args.gaf_base)?;
    let cache_size: isize = -1048576; // 1 GiB
    database.set_pragma("cache_size", cache_size)?;
    let graph = if let Some(gbz_file) = &args.reference {
        Some(serialize::load_from(gbz_file)?)
    } else {
        None
    };

    // Check that the inputs are compatible.
    let alignments = database.graph_name()?;
    if let Some(graph) = &graph {
        let reference = GraphName::from_gbz(graph);
        utils::require_valid_reference(&alignments, &reference)?;
    }

    write_gaf(&database, &alignments, graph.as_ref(), args.chunk_size, args.output.as_deref())?;

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);
    utils::report_peak_memory_usage();

    Ok(())
}

fn write_gaf(
    database: &GAFBase,
    alignments: &GraphName,
    graph: Option<&GBZ>,
    chunk_size: usize,
    output: Option<&std::path::Path>,
) -> Result<()> {
    // Open the output as either a file or stdout.
    let writer: Box<dyn Write + Send> = match output {
        Some(path) if path != std::path::Path::new("-") => {
            let file = File::create(path).map_err(|x| Error::io(format!("Failed to create {}: {}", path.display(), x)))?;
            Box::new(file)
        },
        _ => Box::new(io::stdout()),
    };

    // Decoded ReadSets, with an empty ReadSet signaling the end of input.
    let (to_output, from_decoder) = mpsc::sync_channel(4);

    // Status of the output thread as Result<()>.
    let (to_decoder, from_output) = mpsc::sync_channel(1);

    // Determine header lines first and pass them to the output thread.
    let header_lines = alignments.to_gaf_header_lines();

    // Output thread.
    let output_thread = thread::spawn(move || {
        let mut output = BufWriter::new(writer);
        let mut status = formats::write_gaf_file_header(&mut output)
            .map_err(Error::io);
        if status.is_ok() {
            status = formats::write_header_lines(&header_lines, &mut output)
                .map_err(Error::io);
        }
        while status.is_ok() {
            let read_set: ReadSet = from_decoder.recv().unwrap_or_else(|_| ReadSet::default());
            if read_set.is_empty() {
                break;
            }
            status = read_set.to_gaf(&mut output);
        }
        if status.is_ok() {
            status = output.flush().map_err(Error::io);
        }
        let _ = to_decoder.send(status);
    });

    let mut found_alns = 0;
    let mut rowid = 1; // SQLite row ids start from 1.
    let mut status = Ok(());
    while found_alns < database.alignments() {
        let range = rowid..(rowid + chunk_size);
        let read_set = ReadSet::from_rows(database, range.clone(), graph);
        if let Err(msg) = &read_set {
            status = Err(msg.clone());
            let _ = to_output.send(ReadSet::default()); // Signal end of input.
            break;
        }
        let read_set = read_set.unwrap();
        if read_set.is_empty() {
            status = Err(Error::invalid_data(format!("No reads found in rows {}..{}", range.start, range.end)));
            let _ = to_output.send(ReadSet::default()); // Signal end of input.
            break;
        }
        found_alns += read_set.len();
        let _ = to_output.send(read_set);
        rowid += chunk_size;
    }
    if status.is_ok() {
        let _ = to_output.send(ReadSet::default()); // Signal end of input.
        if found_alns != database.alignments() {
            status = Err(Error::invalid_data(format!("Expected {} alignments, but found {}", database.alignments(), found_alns)));
        }
    }

    // Wait for the output thread to finish.
    let output_result = from_output.recv().unwrap_or(Ok(()));
    let _ = output_thread.join();
    output_result?;

    status
}

//-----------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum SortKey {
    Interval,
    Hash,
}

impl From<SortKey> for KeyType {
    fn from(value: SortKey) -> Self {
        match value {
            SortKey::Interval => KeyType::NodeInterval,
            SortKey::Hash => KeyType::Hash,
        }
    }
}

#[derive(Args)]
struct SortArgs {
    /// GAF alignment file (may be gzip-compressed)
    #[arg(value_name = "input.gaf[.gz]")]
    input: PathBuf,

    /// Output file name (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Parameter preset
    #[arg(long, value_name = "PRESET", value_enum, default_value_t = Preset::Default)]
    preset: Preset,

    /// Sorting key type
    #[arg(short, long, value_name = "TYPE", value_enum, default_value_t = SortKey::Interval)]
    key_type: SortKey,

    /// Number of records per file in initial sort
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = SortParameters::DEFAULT_RECORDS_PER_FILE)]
    records_per_file: usize,

    /// Number of files to merge at once
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = SortParameters::DEFAULT_FILES_PER_MERGE)]
    files_per_merge: usize,

    /// Buffer size for reading/writing records
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = SortParameters::DEFAULT_BUFFER_SIZE)]
    buffer_size: usize,

    /// Number of worker threads
    #[arg(short, long, value_name = "INT", value_parser = parse_quantity, default_value_t = 1)]
    threads: usize,

    /// Use stable sorting (slower but preserves order of equal keys)
    #[arg(short, long)]
    stable: bool,

    /// Print progress information to stderr
    #[arg(short, long)]
    progress: bool,
}

fn sort(args: SortArgs, matches: &clap::ArgMatches) -> Result<()> {
    let start_time = Instant::now();

    let output_file = args.output.unwrap_or_else(|| PathBuf::from("-"));

    let mut params = SortParameters::with_preset(args.preset.as_str())?;

    // Start from the preset and let each option override it only if the user passed
    // it explicitly. Otherwise clap's default values would silently override the preset.
    let explicit = |id: &str| matches.value_source(id) == Some(ValueSource::CommandLine);
    if explicit("key_type") {
        params.key_type = args.key_type.into();
    }
    if explicit("records_per_file") {
        params.records_per_file = args.records_per_file;
    }
    if explicit("files_per_merge") {
        params.files_per_merge = args.files_per_merge;
    }
    if explicit("buffer_size") {
        params.buffer_size = args.buffer_size;
    }
    if explicit("threads") {
        params.threads = args.threads;
    }
    if explicit("stable") {
        params.stable = args.stable;
    }
    if explicit("progress") {
        params.progress = args.progress;
    }

    // Validate the effective parameters.
    if params.records_per_file < 1000 {
        return Err(Error::invalid_query("--records-per-file must be at least 1000"));
    }
    if params.files_per_merge < 2 {
        return Err(Error::invalid_query("--files-per-merge must be at least 2"));
    }
    if params.buffer_size < 1 {
        return Err(Error::invalid_query("--buffer-size must be positive"));
    }

    // Sort the GAF file.
    sort_gaf(&args.input, &output_file, &params)?;

    if params.progress {
        let end_time = Instant::now();
        let seconds = end_time.duration_since(start_time).as_secs_f64();
        eprintln!("Total time: {:.3} seconds", seconds);
        utils::report_peak_memory_usage();
    }

    Ok(())
}

//-----------------------------------------------------------------------------
