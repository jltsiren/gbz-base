use std::path::PathBuf;
use std::time::Instant;
use std::{env, process};

use gbz_base::gaf_sort::{sort_gaf, KeyType, SortParameters};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    let config = Config::new();

    // Sort the GAF file
    sort_gaf(&config.input_file, &config.output_file, &config.params)?;

    if config.params.progress {
        let end_time = Instant::now();
        let seconds = end_time.duration_since(start_time).as_secs_f64();
        eprintln!("Total time: {:.3} seconds", seconds);
    }

    Ok(())
}

//-----------------------------------------------------------------------------

struct Config {
    input_file: PathBuf,
    output_file: PathBuf,
    params: SortParameters,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();
        let header = format!("Usage: {} [options] input.gaf[.gz] > output.gaf", program);

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("o", "output", "output file name (default: stdout)", "FILE");
        opts.optopt(
            "k",
            "key-type",
            "sorting key type: 'interval' (default) or 'hash'",
            "TYPE",
        );
        opts.optopt(
            "r",
            "records-per-file",
            &format!(
                "number of records per file in initial sort (default: {})",
                SortParameters::DEFAULT_RECORDS_PER_FILE
            ),
            "INT",
        );
        opts.optopt(
            "f",
            "files-per-merge",
            &format!(
                "number of files to merge at once (default: {})",
                SortParameters::DEFAULT_FILES_PER_MERGE
            ),
            "INT",
        );
        opts.optopt(
            "b",
            "buffer-size",
            &format!(
                "buffer size for reading/writing records (default: {})",
                SortParameters::DEFAULT_BUFFER_SIZE
            ),
            "INT",
        );
        opts.optopt(
            "t",
            "threads",
            "number of worker threads (default: 1)",
            "INT",
        );
        opts.optflag("s", "stable", "use stable sorting (slower but preserves order of equal keys)");
        opts.optflag("p", "progress", "print progress information to stderr");

        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        // Parse options.
        if matches.opt_present("h") {
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }

        // Parse positional arguments
        if matches.free.len() != 1 {
            eprintln!("Error: Expected 1 positional argument (input file)\n");
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }

        let input_file = PathBuf::from(&matches.free[0]);
        let output_file = if let Some(o) = matches.opt_str("o") {
            PathBuf::from(o)
        } else {
            PathBuf::from("-") // Use stdout if no output file is specified
        };

        let mut params = SortParameters::default();

        // Parse key type
        params.key_type = if let Some(s) = matches.opt_str("k") {
            match s.to_lowercase().as_str() {
                "interval" => KeyType::NodeInterval,
                "hash" => KeyType::Hash,
                _ => {
                    eprintln!("Error: Invalid key type '{}'. Must be 'interval' or 'hash'", s);
                    process::exit(1);
                }
            }
        } else {
            KeyType::NodeInterval
        };

        // Parse records per file
        params.records_per_file = if let Some(s) = matches.opt_str("r") {
            match s.parse::<usize>() {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("Error: Failed to parse --records-per-file: {}", e);
                    process::exit(1);
                }
            }
        } else {
            SortParameters::DEFAULT_RECORDS_PER_FILE
        };

        // Parse files per merge
        params.files_per_merge = if let Some(s) = matches.opt_str("f") {
            match s.parse::<usize>() {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("Error: Failed to parse --files-per-merge: {}", e);
                    process::exit(1);
                }
            }
        } else {
            SortParameters::DEFAULT_FILES_PER_MERGE
        };

        // Parse buffer size
        params.buffer_size = if let Some(s) = matches.opt_str("b") {
            match s.parse::<usize>() {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("Error: Failed to parse --buffer-size: {}", e);
                    process::exit(1);
                }
            }
        } else {
            SortParameters::DEFAULT_BUFFER_SIZE
        };

        // Parse threads
        params.threads = if let Some(s) = matches.opt_str("t") {
            match s.parse::<usize>() {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("Error: Failed to parse --threads: {}", e);
                    process::exit(1);
                }
            }
        } else {
            1
        };

        params.stable = matches.opt_present("s");
        params.progress = matches.opt_present("p");

        // Validate options.
        if params.records_per_file < 1 {
            eprintln!("Error: --records-per-file must be positive");
            process::exit(1);
        }
        if params.files_per_merge < 2 {
            eprintln!("Error: --files-per-merge must be at least 2");
            process::exit(1);
        }
        if params.buffer_size < 1 {
            eprintln!("Error: --buffer-size must be positive");
            process::exit(1);
        }

        Config {
            input_file,
            output_file,
            params,
        }
    }
}

//-----------------------------------------------------------------------------
