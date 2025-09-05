use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{env, fs, process};

use gbz_base::GBZBase;
use gbz_base::utils;

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
    GBZBase::create_from_files(config.gbz_file.as_ref(), config.chains_file(), &config.db_file)?;

    // Statistics.
    let database = GBZBase::open(&config.db_file)?;
    eprintln!("The graph contains {} nodes in {} chains with {} links",
        database.nodes(), database.chains(), database.chain_links()
    );
    eprintln!("There are {} paths representing {} samples, {} haplotypes, and {} contigs",
        database.paths(), database.samples(), database.haplotypes(), database.contigs()
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
    pub gbz_file: PathBuf,
    pub chains_file: Option<PathBuf>,
    pub db_file: PathBuf,
    pub overwrite: bool,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("c", "chains", "top-level chain file (optional)", "FILE");
        opts.optopt("o", "output", "output file name (default: <input>.db)", "FILE");
        opts.optflag("", "overwrite", "overwrite the database file if it exists");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let chains_file = matches.opt_str("c").map(PathBuf::from);

        let mut db_file: Option<PathBuf> = None;
        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(s) = matches.opt_str("o") {
            db_file = Some(PathBuf::from(s));
        }

        let gbz_file = if let Some(s) = matches.free.first() {
            PathBuf::from(s)
        } else {
            let header = format!("Usage: {} [options] graph.gbz", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        if db_file.is_none() {
            // TODO: add_extension is still experimental
            db_file = Some(PathBuf::from(format!("{}.db", gbz_file.display())));
        }

        let overwrite = matches.opt_present("overwrite");

        Config {
            gbz_file,
            chains_file,
            db_file: db_file.unwrap(),
            overwrite,
        }
    }

    fn chains_file(&self) -> Option<&Path> {
        self.chains_file.as_ref().map(|x| x.as_ref())
    }
}

//-----------------------------------------------------------------------------
