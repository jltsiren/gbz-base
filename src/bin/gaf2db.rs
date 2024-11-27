use std::time::Instant;
use std::{env, fs, process};

use gbz_base::GAFBase;
use gbz_base::utils;

use getopts::Options;

//-----------------------------------------------------------------------------

// TODO: Handle gzip-compressed GAF.

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    // Parse arguments.
    let config = Config::new();

    // Check if the database already exists.
    if utils::file_exists(&config.db_file) {
        if config.overwrite {
            eprintln!("Overwriting database {}", config.db_file);
            fs::remove_file(&config.db_file).map_err(|x| x.to_string())?;
        } else {
            return Err(format!("Database {} already exists", config.db_file));
        }
    }

    // Create the database.
    GAFBase::create_from_files(&config.gaf_file, &config.gbwt_file, &config.db_file)?;

    // Statistics.
    let database = GAFBase::open(&config.db_file)?;
    eprintln!(
        "The database contains {} nodes, {} alignments, and {} sequences",
        database.nodes(), database.alignments(), database.sequences()
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
    pub gaf_file: String,
    pub gbwt_file: String,
    pub db_file: String,
    pub overwrite: bool,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();
        let header = format!("Usage: {} [options] alignments.gaf[.gz]", program);

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("g", "gbwt", "GBWT file name (required)", "FILE");
        opts.optopt("o", "output", "output file name (default: <input>.db)", "FILE");
        opts.optflag("", "overwrite", "overwrite the database file if it exists");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let mut db_file: Option<String> = None;
        if matches.opt_present("h") {
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        let gbwt_file = if let Some(s) = matches.opt_str("g") {
            s.clone()
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        if let Some(s) = matches.opt_str("o") {
            db_file = Some(s);
        }

        let gaf_file = if let Some(s) = matches.free.first() {
            s.clone()
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        if db_file.is_none() {
            db_file = Some(format!("{}.db", gaf_file));
        }

        let overwrite = matches.opt_present("overwrite");

        Config {
            gaf_file, gbwt_file,
            db_file: db_file.unwrap(),
            overwrite,
        }
    }
}

//-----------------------------------------------------------------------------
