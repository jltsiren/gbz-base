use std::{env, process};

use gbwt::GBZ;
use gbz_base::GBZBase;
use getopts::Options;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    // Parse arguments.
    let config = Config::new();

    // Check if the database already exists.
    if GBZBase::exists(&config.db_file) {
        return Err(format!("Database {} already exists", config.db_file));
    }

    // Load input graph.
    eprintln!("Loading GBZ graph {}", config.gbz_file);
    let graph: GBZ = serialize::load_from(&config.gbz_file).map_err(|x| x.to_string())?;
    sanity_checks(&graph)?;

    // Create database.
    eprintln!("Creating database {}", config.db_file);
    GBZBase::create(&graph, &config.db_file)?;

    // Statistics.
    let database = GBZBase::open(&config.db_file)?;
    eprintln!(
        "The graph contains {} nodes, {} samples, {} haplotypes, {} contigs, and {} paths",
        database.nodes(), database.samples(), database.haplotypes(), database.contigs(), database.paths()
    );

    Ok(())
}

//-----------------------------------------------------------------------------

pub struct Config {
    pub gbz_file: String,
    pub db_file: String,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("o", "output", "output file name (default: <input>.db)", "FILE");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f.to_string());
                process::exit(1);
            }
        };

        let mut db_file:  Option<String> = None;
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(s) = matches.opt_str("o") {
            db_file = Some(s);
        }

        let gbz_file = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            let header = format!("Usage: {} [options] graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        if db_file.is_none() {
            db_file = Some(format!("{}.db", gbz_file));
        }

        Config {
            gbz_file,
            db_file: db_file.unwrap(),
        }
    }
}

//-----------------------------------------------------------------------------

fn sanity_checks(graph: &GBZ) -> Result<(), String> {
    if !graph.has_metadata() {
        return Err("The graph does not contain metadata".to_string());
    }
    let metadata = graph.metadata().unwrap();

    if !metadata.has_path_names() {
        return Err("The metadata does not contain path names".to_string());
    }
    if !metadata.has_sample_names() {
        return Err("The metadata does not contain sample names".to_string());
    }
    if !metadata.has_contig_names() {
        return Err("The metadata does not contain contig names".to_string());
    }

    Ok(())
}

//-----------------------------------------------------------------------------
