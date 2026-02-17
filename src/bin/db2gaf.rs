use std::io::{Write, BufWriter};
use std::path::PathBuf;
use std::sync::mpsc;
use std::time::Instant;
use std::{env, io, process, thread};

use gbz::GBZ;

use gbz_base::{GAFBase, ReadSet};
use gbz_base::{formats, utils};

use pggname::GraphName;

use simple_sds::serialize;

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start_time = Instant::now();

    let config = Config::new();

    // Inputs.
    let database = GAFBase::open(&config.gaf_base_file)?;
    let graph = if let Some(gbz_file) = &config.gbz_file {
        Some(serialize::load_from(gbz_file).map_err(|x| x.to_string())?)
    } else {
        None
    };

    // Check that the inputs are compatible.
    let alignments = database.graph_name()?;
    if let Some(graph) = &graph {
        let reference = GraphName::from_gbz(graph);
        let result = utils::require_valid_reference(&alignments, &reference);
        if let Err(e) = result {
            // Print the error manually, as it contains multiple lines.
            eprint!("Error: {}", e);
            process::exit(1);
        }
    }

    write_gaf(&database, &alignments, graph.as_ref(), &config)?;

    let end_time = Instant::now();
    let seconds = end_time.duration_since(start_time).as_secs_f64();
    eprintln!("Used {:.3} seconds", seconds);

    Ok(())
}

//-----------------------------------------------------------------------------

fn write_gaf(database: &GAFBase, alignments: &GraphName, graph: Option<&GBZ>, config: &Config) -> Result<(), String> {
    // Decoded ReadSets, with an empty ReadSet signaling the end of input.
    let (to_output, from_decoder) = mpsc::sync_channel(4);

    // Status of the output thread as Result<(), String>.
    let (to_decoder, from_output) = mpsc::sync_channel(1);

    // Determine header lines first and pass them to the output thread.
    let header_lines = alignments.to_gaf_header_lines();

    // Output thread.
    let output_thread = thread::spawn(move || {
        let mut output = BufWriter::new(io::stdout().lock());
        let mut status = formats::write_gaf_file_header(&mut output)
            .map_err(|e| e.to_string());
        if status.is_ok() {
            status = formats::write_header_lines(&header_lines, &mut output)
                .map_err(|e| e.to_string());
        }
        while status.is_ok() {
            let read_set: ReadSet = from_decoder.recv().unwrap_or(ReadSet::default());
            if read_set.is_empty() {
                break;
            }
            status = read_set.to_gaf(&mut output);
        }
        if status.is_ok() {
            status = output.flush().map_err(|e| e.to_string());
        }
        let _ = to_decoder.send(status);
    });

    let mut found_alns = 0;
    let mut rowid = 1; // SQLite row ids start from 1.
    let mut status = Ok(());
    while found_alns < database.alignments() {
        let range = rowid..(rowid + config.chunk_size);
        let read_set = ReadSet::from_rows(database, range.clone(), graph);
        if let Err(msg) = &read_set {
            status = Err(msg.clone());
            let _ = to_output.send(ReadSet::default()); // Signal end of input.
            break;
        }
        let read_set = read_set.unwrap();
        if read_set.is_empty() {
            status = Err(format!("No reads found in rows {}..{}", range.start, range.end));
            let _ = to_output.send(ReadSet::default()); // Signal end of input.
            break;
        }
        found_alns += read_set.len();
        let _ = to_output.send(read_set);
        rowid += config.chunk_size;
    }
    if status.is_ok() {
        let _ = to_output.send(ReadSet::default()); // Signal end of input.
        if found_alns != database.alignments() {
            status = Err(format!("Expected {} alignments, but found {}", database.alignments(), found_alns));
        }
    }

    // Wait for the output thread to finish.
    let output_result = from_output.recv().unwrap_or(Ok(()));
    let _ = output_thread.join();
    output_result?;

    status
}

//-----------------------------------------------------------------------------

struct Config {
    gaf_base_file: PathBuf,
    gbz_file: Option<PathBuf>,
    chunk_size: usize,
}

impl Config {
    // Default number of blocks in a single ReadSet.
    const CHUNK_SIZE: usize = 100;

    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();
        let header = format!("Usage: {} [options] gaf_base.db > output.gaf", program);

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("r", "reference", "use this GBZ graph as the reference", "FILE");
        opts.optopt("c", "chunk-size", "chunk size in blocks (default: 10)", "INT");
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
        let gaf_base_file = if let Some(s) = matches.free.first() {
            PathBuf::from(s)
        } else {
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        };
        let gbz_file = if let Some(s) = matches.opt_str("r") {
            Some(PathBuf::from(s))
        } else {
            None
        };
        let chunk_size = if let Some(s) = matches.opt_str("c") {
            match s.parse::<usize>() {
                Ok(x) => x,
                Err(e) => {
                    eprintln!("Failed to parse --chunk-size: {}", e);
                    process::exit(1);
                }
            }
        } else {
            Config::CHUNK_SIZE
        };

        // Validate options.
        if chunk_size < 1 {
            eprintln!("--chunk-size must be positive");
            process::exit(1);
        }

        Config {
            gaf_base_file,
            gbz_file,
            chunk_size,
        }
    }
}

//-----------------------------------------------------------------------------
