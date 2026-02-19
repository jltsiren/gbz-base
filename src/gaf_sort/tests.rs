use std::io::BufRead;
use super::*;

//-----------------------------------------------------------------------------
// Helpers
//-----------------------------------------------------------------------------

/// Reads all non-header data lines from a GAF file.
/// Lines are returned as-is (including the trailing newline), matching what
/// GAFRecord stores internally, so that key_of() computes correct hash keys.
fn read_data_lines(path: &Path) -> Vec<Vec<u8>> {
    let mut reader = utils::open_file(path).unwrap();
    let mut lines = Vec::new();
    let mut buf = Vec::new();
    loop {
        buf.clear();
        match reader.read_until(b'\n', &mut buf) {
            Ok(0) => break,
            Ok(_) => {
                // Determine content length, excluding the trailing newline.
                let content_end = if buf.last() == Some(&b'\n') { buf.len() - 1 } else { buf.len() };
                if content_end > 0 && buf[0] != b'@' {
                    lines.push(buf.clone());
                }
            }
            Err(e) => panic!("read error: {}", e),
        }
    }
    lines
}

/// Reads header lines from a GAF file.
fn read_headers(path: &Path) -> Vec<String> {
    let mut reader = utils::open_file(path).unwrap();
    formats::read_gaf_header_lines(&mut reader).unwrap()
}

/// Runs sort_gaf on a test input file and returns the path to the output file.
fn run_sort(input: &'static str, params: &SortParameters) -> PathBuf {
    let input_path = utils::get_test_data(input);
    let output_path = serialize::temp_file_name("gaf-sort-test");
    sort_gaf(&input_path, &output_path, params).expect("sort_gaf failed");
    output_path
}

/// Returns the sort key for a GAF line (same logic as GAFRecord).
fn key_of(line: &[u8], key_type: KeyType) -> u64 {
    GAFRecord::new(line.to_vec(), key_type).key
}

/// Asserts that lines are in non-decreasing key order.
fn assert_sorted(lines: &[Vec<u8>], key_type: KeyType) {
    let mut prev = 0u64;
    for (i, line) in lines.iter().enumerate() {
        let k = key_of(line, key_type);
        assert!(
            k >= prev,
            "record {} is out of order: prev_key={:#x}, key={:#x}",
            i, prev, k
        );
        prev = k;
    }
}

/// Asserts that the lines match.
fn assert_lines_equal(lines1: &[Vec<u8>], lines2: &[Vec<u8>], first: &str, second: &str) {
    assert_eq!(lines1.len(), lines2.len(), "line count differs between {} and {}", first, second);
    for (i, (line1, line2)) in lines1.iter().zip(lines2.iter()).enumerate() {
        assert_eq!(
            line1, line2,
            "line {} differs between {} and {}: {:?} != {:?}",
            i, first, second, String::from_utf8_lossy(line1), String::from_utf8_lossy(line2)
        );
    }
}

/// Number of data records in shuffled.gaf (and shuffled.gaf.gz).
const RECORD_COUNT: usize = 12439;

//-----------------------------------------------------------------------------
// Unit tests
//-----------------------------------------------------------------------------

#[test]
fn test_node_interval_key() {
    let line = b"query\t100\t0\t100\t+\t>1>2>3\t300\t0\t300\t100\t100\t60".to_vec();
    let record = GAFRecord::new(line, KeyType::NodeInterval);
    // Handles: >1=2, >2=4, >3=6
    // min_handle=2, max_handle=6 -> (2 << 32) | 6
    let min_handle = support::encode_node(1, Orientation::Forward) as u64;
    let max_handle = support::encode_node(3, Orientation::Forward) as u64;
    assert_eq!(record.key, (min_handle << 32) | max_handle);

    let line = b"query\t100\t0\t100\t+\t>5<10>15\t300\t0\t300\t100\t100\t60".to_vec();
    let record = GAFRecord::new(line, KeyType::NodeInterval);
    // Handles: >5=10, <10=21, >15=30
    // min_handle=10, max_handle=30 -> (10 << 32) | 30
    let min_handle = support::encode_node(5, Orientation::Forward) as u64;
    let max_handle = support::encode_node(15, Orientation::Forward) as u64;
    assert_eq!(record.key, (min_handle << 32) | max_handle);
}

#[test]
fn test_hash_key() {
    let line = b"query\t100\t0\t100\t+\t>1>2>3\t300\t0\t300\t100\t100\t60".to_vec();
    let record = GAFRecord::new(line, KeyType::Hash);
    assert_ne!(record.key, GAFRecord::MISSING_KEY);
}

#[test]
fn test_serialization() {
    let line = b"query\t100\t0\t100\t+\t>1>2>3\t300\t0\t300\t100\t100\t60".to_vec();
    let record = GAFRecord::new(line.clone(), KeyType::NodeInterval);

    let mut buffer = Vec::new();
    record.serialize(&mut buffer).unwrap();

    let mut cursor = std::io::Cursor::new(buffer);
    let deserialized = GAFRecord::deserialize(&mut cursor).unwrap();

    assert_eq!(record.key, deserialized.key);
    assert_eq!(record.value, deserialized.value);
}

//-----------------------------------------------------------------------------
// sort_gaf correctness tests
//-----------------------------------------------------------------------------

// Single batch: all records fit in one batch and are sorted directly to the
// output file, exercising the sort_to_output path.
#[test]
fn sort_single_batch() {
    let params = SortParameters {
        records_per_file: RECORD_COUNT + 1,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(lines.len(), RECORD_COUNT);
    assert_sorted(&lines, params.key_type);
}

// Multiple batches merged in a single round: ceil(12439 / 5000) = 3 initial
// temp files, all merged at once because files_per_merge = 32.
// Exercises sort_to_temp + merge_to_output.
#[test]
fn sort_multi_batch_single_merge() {
    let params = SortParameters {
        records_per_file: 5000,
        files_per_merge: 32,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(lines.len(), RECORD_COUNT);
    assert_sorted(&lines, params.key_type);
}

// Multiple batches and multiple merge rounds: ceil(12439 / 1000) = 13 initial
// temp files reduced to 2 via several rounds of merge_files before the final
// merge_to_output.
// Exercises sort_to_temp + merge_files + merge_to_output.
#[test]
fn sort_multi_batch_multi_round() {
    let params = SortParameters {
        records_per_file: 1000,
        files_per_merge: 2,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(lines.len(), RECORD_COUNT);
    assert_sorted(&lines, params.key_type);
}

// Repeat the above with two threads.
#[test]
fn sort_multithreaded() {
    let params = SortParameters {
        records_per_file: 1000,
        files_per_merge: 2,
        threads: 2,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(lines.len(), RECORD_COUNT);
    assert_sorted(&lines, params.key_type);
}

// Gzipped input: shuffled.gaf.gz should produce the same set of records as
// the plain shuffled.gaf.
#[test]
fn sort_gzipped_input() {
    let params = SortParameters {
        records_per_file: 5000,
        ..SortParameters::default()
    };
    let out_plain = run_sort("shuffled.gaf", &params);
    let out_gz = run_sort("shuffled.gaf.gz", &params);

    let mut lines_plain = read_data_lines(&out_plain);
    let mut lines_gz = read_data_lines(&out_gz);
    let _ = fs::remove_file(&out_plain);
    let _ = fs::remove_file(&out_gz);

    lines_plain.sort_unstable();
    lines_gz.sort_unstable();
    assert_lines_equal(&lines_plain, &lines_gz, "plain", "gzipped");
}

// Hash key: output should be sorted by hash of the record value.
#[test]
fn sort_hash_key() {
    let params = SortParameters {
        key_type: KeyType::Hash,
        records_per_file: 5000,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(lines.len(), RECORD_COUNT);
    assert_sorted(&lines, KeyType::Hash);
}

// Header lines from the input should appear unchanged at the top of the output.
#[test]
fn sort_preserves_headers() {
    let input_path = utils::get_test_data("shuffled.gaf");
    let input_headers = read_headers(&input_path);

    let params = SortParameters {
        records_per_file: 5000,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let output_headers = read_headers(&output);
    let _ = fs::remove_file(&output);

    assert_eq!(input_headers, output_headers);
}

// The sorted output should contain exactly the same set of records as the
// input (no records lost, duplicated, or corrupted).
#[test]
fn sort_preserves_all_records() {
    let input_path = utils::get_test_data("shuffled.gaf");
    let mut input_lines = read_data_lines(&input_path);

    let params = SortParameters {
        records_per_file: 5000,
        ..SortParameters::default()
    };
    let output = run_sort("shuffled.gaf", &params);
    let mut output_lines = read_data_lines(&output);
    let _ = fs::remove_file(&output);

    input_lines.sort_unstable();
    output_lines.sort_unstable();
    assert_lines_equal(&input_lines, &output_lines, "input", "output");
}

// Different parameter configurations should produce identical sets of records.
#[test]
fn sort_consistent_across_configs() {
    // Single batch, multi-batch single merge, multi-batch multi-round.
    let params = [
        SortParameters { records_per_file: RECORD_COUNT + 1, ..SortParameters::default() },
        SortParameters { records_per_file: 5000, files_per_merge: 32, ..SortParameters::default() },
        SortParameters { records_per_file: 1000, files_per_merge: 2, ..SortParameters::default() },
    ];

    let outputs: Vec<PathBuf> = params.iter().map(|p| run_sort("shuffled.gaf", p)).collect();
    let sorted_lines: Vec<Vec<Vec<u8>>> = outputs.iter().map(|path| {
        let mut lines = read_data_lines(path);
        lines.sort_unstable();
        lines
    }).collect();
    for path in &outputs {
        let _ = fs::remove_file(path);
    }

    assert_lines_equal(&sorted_lines[0], &sorted_lines[1], "config 0", "config 1");
    assert_lines_equal(&sorted_lines[0], &sorted_lines[2], "config 0", "config 2");
}

//-----------------------------------------------------------------------------
