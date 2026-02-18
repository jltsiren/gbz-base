use std::path::PathBuf;
use std::time::Instant;
use std::{env, fs, process};

use gbz_base::{GBZBase, GraphInterface};
use gbz_base::{GAFBase, GAFBaseParams, GraphReference};
use gbz_base::db::FileType;
use gbz_base::{db, utils};

use gbz::GBZ;

use simple_sds::serialize;

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = Config::new();

    // Check if the database already exists.
    if utils::file_exists(&config.db_file) {
        if config.overwrite {
            eprintln!("Overwriting database {}", config.db_file.display());
            fs::remove_file(&config.db_file).map_err(|x| x.to_string())?;
        } else {
            return Err(format!("Database {} already exists", config.db_file.display()));
        }
    }

    // Create the database.
    if let Some(graph_file) = &config.graph_file {
        match db::identify_file(graph_file) {
            FileType::Gbz => {
                eprintln!("Loading GBZ graph {}", graph_file.display());
                let graph: GBZ = serialize::load_from(graph_file).map_err(|x| x.to_string())?;
                GAFBase::create_from_files(
                    &config.gaf_file, &config.gbwt_file, &config.db_file,
                    GraphReference::Gbz(&graph), &config.params
                )?;
            },
            FileType::Version(v) => {
                if v != GBZBase::VERSION {
                    let msg = format!("File {} is {}; expected {}", graph_file.display(), v, GBZBase::VERSION);
                    return Err(msg);
                }
                eprintln!("Opening GBZ-base {}", graph_file.display());
                let database = GBZBase::open(graph_file)?;
                let mut graph = GraphInterface::new(&database)?;
                GAFBase::create_from_files(
                    &config.gaf_file, &config.gbwt_file, &config.db_file,
                    GraphReference::Db(&mut graph), &config.params
                )?;
            },
            _ => {
                return Err(format!("File {} is not a valid graph", graph_file.display()));
            }
        };
    } else {
        GAFBase::create_from_files(&config.gaf_file, &config.gbwt_file, &config.db_file, GraphReference::None, &config.params)?;
    }

    // Statistics.
    let database = GAFBase::open(&config.db_file)?;
    eprintln!(
        "The database contains {} nodes and {} alignments in {} blocks",
        database.nodes(), database.alignments(), database.blocks()
    );
    let size = database.file_size().unwrap_or(String::from("unknown"));
    eprintln!("Final database size: {}", size);

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
}

//-----------------------------------------------------------------------------

struct Config {
    pub gaf_file: PathBuf,
    pub gbwt_file: PathBuf,
    pub graph_file: Option<PathBuf>,
    pub db_file: PathBuf,
    pub overwrite: bool,
    pub params: GAFBaseParams,
}

impl Config {
    pub fn new() -> Config {
        let mut params = GAFBaseParams::default();

        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();
        let header = format!("Usage: {} [options] alignments.gaf[.gz]", program);

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("g", "gbwt", "GBWT file name (required)", "FILE");
        opts.optopt("r", "ref-free", "build a reference-free GAF-base using this graph", "FILE");
        let block_desc = format!("number of alignments per block (default: {})", params.block_size);
        opts.optopt("b", "block-size", &block_desc, "INT");
        opts.optflag("", "no-quality", "do not store quality strings");
        opts.optflag("", "no-optional", "do not store unsupported optional fields");
        opts.optopt("o", "output", "output file name (default: <input>.db)", "FILE");
        opts.optflag("", "overwrite", "overwrite the database file if it exists");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let mut db_file: Option<PathBuf> = None;
        if matches.opt_present("h") {
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        let gbwt_file = if let Some(s) = matches.opt_str("g") {
            PathBuf::from(s)
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        let graph_file = if let Some(s) = matches.opt_str("r") {
            params.reference_free = true;
            Some(PathBuf::from(s))
        } else {
            None
        };
        if let Some(s) = matches.opt_str("o") {
            db_file = Some(PathBuf::from(s));
        }

        let gaf_file = if let Some(s) = matches.free.first() {
            PathBuf::from(s)
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        if db_file.is_none() {
            let mut name = gaf_file.clone();
            name.add_extension("db");
            db_file = Some(name);
        }

        // Parameters.
        if let Some(s) = matches.opt_str("b") {
            match s.parse::<usize>() {
                Ok(size) => params.block_size = size,
                Err(_) => {
                    eprintln!("Invalid block size: {}", s);
                    process::exit(1);
                }
            }
        }
        if matches.opt_present("no-quality") {
            params.store_quality_strings = false;
        }
        if matches.opt_present("no-optional") {
            params.store_optional_fields = false;
        }

        let overwrite = matches.opt_present("overwrite");

        Config {
            gaf_file, gbwt_file, graph_file,
            db_file: db_file.unwrap(),
            overwrite,
            params,
        }
    }
}

//-----------------------------------------------------------------------------
