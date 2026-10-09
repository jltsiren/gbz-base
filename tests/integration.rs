//! Integration tests for binaries.

use std::collections::HashSet;
use std::fs;
use std::ops::Range;
use std::path::PathBuf;
use std::process::{Command, Output, ExitStatus};

use gbz_base::{GBZBase, GraphInterface, SnarlOutput, PathIndex};
use gbz_base::{GAFBase, GAFBaseParams, GraphReference};
use gbz_base::{Subgraph, SubgraphQuery, HaplotypeOutput, ReadSet, AlignmentOutput};
use gbz_base::gaf_sort::{sort_gaf, SortParameters, KeyType};
use gbz_base::{formats, utils};
use gbz_base::Error;

use gbz::{GBZ, FullPathName, Orientation};
use gbz::support::{self, Chains};
use simple_sds::serialize;

use rand::prelude::IndexedRandom;
use rand::Rng;

//-----------------------------------------------------------------------------

// TODO: Should this be in a library somewhere?
struct TempFileHandler {
    files: Vec<PathBuf>,
}

impl TempFileHandler {
    fn new() -> Self {
        TempFileHandler { files: Vec::new() }
    }

    fn add_file(&mut self, file: &PathBuf) {
        self.files.push(file.clone());
    }

    fn new_file(&mut self, name_part: &str) -> PathBuf {
        let file = serialize::temp_file_name(name_part);
        assert!(!file.exists(), "Temporary file {} already exists", file.display());
        self.files.push(file.clone());
        file
    }

    fn create_copy(&mut self, source: &PathBuf) -> PathBuf {
        let copy = serialize::temp_file_name("copy");
        assert!(!copy.exists(), "Temporary file {} already exists", copy.display());
        fs::copy(source, &copy).expect("Failed to copy file");
        self.files.push(copy.clone());
        copy
    }
}

impl Drop for TempFileHandler {
    fn drop(&mut self) {
        for file in &self.files {
            let _ = fs::remove_file(file);
        }
    }
}

fn compare_files(file1: &PathBuf, file2: &PathBuf) -> bool {
    let content1 = fs::read(file1).expect(&format!("Failed to read file 1: {}", file1.display()));
    let content2 = fs::read(file2).expect(&format!("Failed to read file 2: {}", file2.display()));
    content1 == content2
}

fn count_lines_in_file(file: &PathBuf) -> usize {
    let content = fs::read_to_string(file).expect(&format!("Failed to read file: {}", file.display()));
    content.lines().count()
}

fn count_lines(data: &[u8]) -> usize {
    let content = String::from_utf8_lossy(data);
    content.lines().count()
}

//-----------------------------------------------------------------------------

fn get_binary_path(binary_name: &str) -> std::path::PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push(binary_name);
    path
}

fn build_gbz_base_truth(
    temp_files: &mut TempFileHandler,
    graph_file: &PathBuf, chains_file: Option<&PathBuf>
) -> PathBuf {
    let truth_file = temp_files.new_file("gbz-base-truth");
    let result = GBZBase::create_from_files(
        graph_file, chains_file.map(|p| p.as_path()), &truth_file
    );
    assert!(result.is_ok(), "Failed to create truth database: {}", result.unwrap_err());
    truth_file
}

fn sort_gaf_to_temp(temp_files: &mut TempFileHandler, input_file: &PathBuf) -> PathBuf {
    let sorted_file = temp_files.new_file("sorted");
    sort_gaf(input_file, &sorted_file, &SortParameters::default())
        .expect("Failed to sort GAF file");
    sorted_file
}

fn build_gaf_base_truth(
    temp_files: &mut TempFileHandler,
    sorted_input: &PathBuf, gbwt_file: Option<&PathBuf>, graph_file: Option<&PathBuf>,
    params: &GAFBaseParams
) -> PathBuf {
    let truth_file = temp_files.new_file("gaf-base-truth");

    let result = if let Some(graph_file) = graph_file {
        let graph: GBZ = serialize::load_from(graph_file).expect("Failed to read graph file");
        GAFBase::create_from_files(
            sorted_input, gbwt_file.map(|p| p.as_path()), &truth_file,
            GraphReference::Gbz(&graph), params
        )
    } else {
        GAFBase::create_from_files(
            sorted_input, gbwt_file.map(|p| p.as_path()), &truth_file,
            GraphReference::None, params
        )
    };
    assert!(result.is_ok(), "Failed to create truth database: {}", result.unwrap_err());

    truth_file
}

//-----------------------------------------------------------------------------

// Tests for `gbz-base construct`.

// We assume that the graph file is a temporary file if no output file is provided.
// Otherwise there may be conflicts with multiple tests running in parallel.
fn run_gbz_base_construct(
    temp_files: &mut TempFileHandler,
    graph_file: &PathBuf, chains_file: Option<&PathBuf>,
    overwrite: bool, output_file: Option<&PathBuf>
) -> (PathBuf, ExitStatus) {
    let mut args = vec![String::from("construct")];
    let output_file = match output_file {
        Some(path) => {
            args.push(String::from("--output"));
            args.push(path.to_str().unwrap().to_string());
            path.clone()
        },
        None => {
            let path = graph_file.with_added_extension("db");
            temp_files.add_file(&path);
            path
        },
    };
    if overwrite {
        args.push(String::from("--overwrite"));
    }
    if let Some(chains_file) = chains_file {
        args.push(String::from("--chains"));
        args.push(chains_file.to_str().unwrap().to_string());
    }
    args.push(graph_file.to_str().unwrap().to_string());

    let binary = get_binary_path("gbz-base");
    let result = Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gbz-base construct");
    (output_file, result.status)
}

#[test]
fn gbz_base_construct() {
    let mut temp_files = TempFileHandler::new();
    let original_graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let graph_file = temp_files.create_copy(&original_graph_file);
    let chains_file = None;

    // We expect that the GBZ-base built with `gbz-base construct` is identical to
    // one built with a direct library call.
    let truth_file = build_gbz_base_truth(&mut temp_files, &graph_file, chains_file);
    let (output_file, status) = run_gbz_base_construct(
        &mut temp_files, &graph_file, chains_file, false, None
    );
    assert!(status.success(), "gbz-base construct failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
    let output_metadata = fs::metadata(&output_file).expect("Failed to read output file metadata");
    let old_timestamp = output_metadata.modified().expect("Failed to get output file timestamp");

    // Try rebuilding the database without overwrite, expect failure.
    let (_, status) = run_gbz_base_construct(
        &mut temp_files, &graph_file, chains_file, false, None
    );
    assert!(!status.success(), "gbz-base construct should have failed without overwrite");

    // Rebuild with overwrite, expect success.
    let (_, status) = run_gbz_base_construct(
        &mut temp_files, &graph_file, chains_file, true, None
    );
    assert!(status.success(), "gbz-base construct failed with overwrite with status: {}", status);

    // Check that we did actually overwrite the output file.
    let output_metadata = fs::metadata(&output_file).expect("Failed to read output file metadata");
    let new_timestamp = output_metadata.modified().expect("Failed to get output file timestamp");
    assert!(new_timestamp > old_timestamp, "Output file timestamp did not change after overwrite");
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file after overwrite");

    // Try specifying an explicit output file name.
    let explicit_output_file = temp_files.new_file("gbz-base-output");
    let (_, status) = run_gbz_base_construct(
        &mut temp_files, &graph_file, chains_file, false, Some(&explicit_output_file)
    );
    assert!(status.success(), "gbz-base construct failed with explicit output file with status: {}", status);
    assert!(compare_files(&explicit_output_file, &truth_file), "Explicit output file does not match truth file");
}

#[test]
fn gbz_base_construct_with_chains() {
    let mut temp_files = TempFileHandler::new();
    let original_graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let graph_file = temp_files.create_copy(&original_graph_file);
    let chains_file = Some(utils::get_test_data("micb-kir3dl1.chains"));

    let truth_file = build_gbz_base_truth(&mut temp_files, &graph_file, chains_file.as_ref());
    let (output_file, status) = run_gbz_base_construct(
        &mut temp_files, &graph_file, chains_file.as_ref(), false, None
    );
    assert!(status.success(), "gbz-base construct failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
}

//-----------------------------------------------------------------------------

// Tests for `gaf-base sort`.

fn run_gaf_base_sort(input_file: &PathBuf, output_file: Option<&PathBuf>, params: &SortParameters) -> Output {
    let mut args = vec![String::from("sort")];
    if let Some(output_file) = output_file {
        args.push(String::from("--output"));
        args.push(output_file.to_str().unwrap().to_string());
    }
    args.push(String::from("--key-type"));
    match params.key_type {
        KeyType::NodeInterval => args.push(String::from("interval")),
        KeyType::Hash => args.push(String::from("hash")),
    }
    args.push(String::from("--records-per-file"));
    args.push(params.records_per_file.to_string());
    args.push(String::from("--files-per-merge"));
    args.push(params.files_per_merge.to_string());
    args.push(String::from("--buffer-size"));
    args.push(params.buffer_size.to_string());
    args.push(String::from("--threads"));
    args.push(params.threads.to_string());
    if params.stable {
        args.push(String::from("--stable"));
    }
    if params.progress {
        args.push(String::from("--progress"));
    }
    args.push(input_file.to_str().unwrap().to_string());

    let binary = get_binary_path("gaf-base");
    Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gaf-base sort")
}

fn run_gaf_base_sort_with_preset(
    input_file: &PathBuf, output_file: &PathBuf,
    preset: &str, other_args: &[String]
) {
    let mut args = vec![String::from("sort")];
    args.push(String::from("--output"));
    args.push(output_file.to_str().unwrap().to_string());
    args.push(String::from("--preset"));
    args.push(preset.to_string());
    args.extend_from_slice(other_args);
    args.push(input_file.to_str().unwrap().to_string());

    let binary = get_binary_path("gaf-base");
    let result = Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gaf-base sort with preset");
    assert!(result.status.success(), "gaf-base sort with preset {} and args {:?} failed with status: {}", preset, other_args, result.status);
}

fn count_merge_rounds(stderr: &[u8]) -> usize {
    let stderr_str = String::from_utf8_lossy(stderr);
    let lines = stderr_str.lines();
    let mut rounds = 0;
    for line in lines {
        if line.starts_with("Round") {
            if line.ends_with("files per batch") {
                rounds += 1;
            }
        } else if line.starts_with("Starting the final merge") {
            rounds += 1;
        }
    }
    rounds
}

#[test]
fn gaf_base_sort() {
    let mut temp_files = TempFileHandler::new();
    let input_file = utils::get_test_data("shuffled.gaf");
    let params = SortParameters::default();

    let truth_file = temp_files.new_file("sorted");
    sort_gaf(&input_file, &truth_file, &params).expect("Failed to sort GAF file for truth");
    let truth = fs::read(&truth_file).expect("Failed to read truth file");

    // Check that the default arguments correspond to the default parameters.
    let default_args = vec![String::from("sort"), input_file.to_str().unwrap().to_string()];
    let binary = get_binary_path("gaf-base");
    let default_output = Command::new(&binary)
        .args(&default_args)
        .output()
        .expect("Failed to execute gaf-base sort with default args");
    assert!(default_output.status.success(), "gaf-base sort failed with default args");
    let correct_output = default_output.stdout == truth; // These are large, so we don't want to print them out.
    assert!(correct_output, "gaf-base sort produced incorrect output with default args");
    assert!(default_output.stderr.is_empty(), "gaf-base sort produced stderr output with default args");

    let explicit_output = run_gaf_base_sort(&input_file, None, &params);
    assert!(explicit_output.status.success(), "gaf-base sort failed with explicit args");
    let correct_output = explicit_output.stdout == truth;
    assert!(correct_output, "gaf-base sort produced incorrect output with explicit args");
    assert!(explicit_output.stderr.is_empty(), "gaf-base sort produced stderr output with explicit args");

    let output_file = temp_files.new_file("sorted");
    let file_output = run_gaf_base_sort(&input_file, Some(&output_file), &params);
    assert!(file_output.status.success(), "gaf-base sort failed with file output");
    assert!(file_output.stdout.is_empty(), "gaf-base sort produced stdout output with file output");
    assert!(file_output.stderr.is_empty(), "gaf-base sort produced stderr output with file output");
    assert!(compare_files(&output_file, &truth_file), "gaf-base sort produced incorrect output with file output");

    let params = SortParameters { progress: true, ..params };
    let progress_output = run_gaf_base_sort(&input_file, None, &params);
    assert!(progress_output.status.success(), "gaf-base sort failed with --progress");
    let correct_output = progress_output.stdout == truth;
    assert!(correct_output, "gaf-base sort produced incorrect output with --progress");
    assert!(!progress_output.stderr.is_empty(), "gaf-base sort did not produce stderr output with --progress");
}

#[test]
fn gaf_base_sort_stable() {
    let mut temp_files = TempFileHandler::new();
    let input_file = utils::get_test_data("shuffled.gaf");
    let params = SortParameters { stable: true, ..SortParameters::default() };

    let truth_file = temp_files.new_file("stable");
    sort_gaf(&input_file, &truth_file, &params).expect("Failed to sort GAF file for truth");
    let truth = fs::read(&truth_file).expect("Failed to read truth file");

    let stable_output = run_gaf_base_sort(&input_file, None, &params);
    assert!(stable_output.status.success(), "gaf-base sort failed with --stable");
    let correct_output = stable_output.stdout == truth;
    assert!(correct_output, "gaf-base sort produced incorrect output with --stable");
    assert!(stable_output.stderr.is_empty(), "gaf-base sort produced stderr output with --stable");
}

#[test]
fn gaf_base_sort_shuffle() {
    let mut temp_files = TempFileHandler::new();
    let input_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let params = SortParameters { key_type: KeyType::Hash, ..SortParameters::default() };

    let truth_file = temp_files.new_file("shuffle");
    sort_gaf(&input_file, &truth_file, &params).expect("Failed to shuffle GAF file for truth");
    let truth = fs::read(&truth_file).expect("Failed to read truth file");

    let shuffled_output = run_gaf_base_sort(&input_file, None, &params);
    assert!(shuffled_output.status.success(), "gaf-base sort failed with --key-type hash");
    let correct_output = shuffled_output.stdout == truth;
    assert!(correct_output, "gaf-base sort produced incorrect output with --key-type hash");
    assert!(shuffled_output.stderr.is_empty(), "gaf-base sort produced stderr output with --key-type hash");
}

fn expect_merge_rounds(input_file: &PathBuf, truth: &[u8], params: &SortParameters, expected_rounds: usize) {
    let output = run_gaf_base_sort(input_file, None, params);
    assert!(output.status.success(), "gaf-base sort failed with {} merge round(s)", expected_rounds);
    let correct_output = output.stdout == truth;
    assert!(correct_output, "gaf-base sort produced incorrect output with {} merge round(s)", expected_rounds);
    let rounds = count_merge_rounds(&output.stderr);
    assert_eq!(rounds, expected_rounds, "wrong number of merge rounds");
}

#[test]
fn gaf_base_sort_params() {
    let mut temp_files = TempFileHandler::new();
    // The input file has 12439 GAF records.
    let input_file = utils::get_test_data("shuffled.gaf");
    let params = SortParameters { stable: true, ..Default::default() };

    // We can reuse the truth, as we use stable sorting.
    let truth_file = temp_files.new_file("sort-params");
    sort_gaf(&input_file, &truth_file, &params).expect("Failed to sort GAF file for truth");
    let truth = fs::read(&truth_file).expect("Failed to read truth file");

    // We didn't want to print progress when sorting the truth.
    let params = SortParameters { progress: true, ..params };

    // Single batch.
    {
        expect_merge_rounds(&input_file, &truth, &params, 0);
    }

    // One round of merges.
    {
        let params = SortParameters { records_per_file: 2000, ..params };
        expect_merge_rounds(&input_file, &truth, &params, 1);
    }

    // Two rounds of merges.
    {
        let params = SortParameters { records_per_file: 2000, files_per_merge: 3, ..params };
        expect_merge_rounds(&input_file, &truth, &params, 2);
    }

    // Multithreaded sort with two rounds of merges.
    {
        let params = SortParameters {
            records_per_file: 2000, files_per_merge: 3, threads: 2, progress: false, ..params
        };
        let output = run_gaf_base_sort(&input_file, None, &params);
        assert!(output.status.success(), "gaf-base sort failed with 2 threads");
        let correct_output = output.stdout == truth;
        assert!(correct_output, "gaf-base sort produced incorrect output with 2 threads");
        assert!(output.stderr.is_empty(), "gaf-base sort produced stderr output with 2 threads");
    }
}

#[test]
fn gaf_base_sort_presets() {
    let mut temp_files = TempFileHandler::new();
    let input_file = utils::get_test_data("shuffled.gaf");

    for preset in SortParameters::PRESETS {
        // First we use just the preset itself.
        let mut params = SortParameters::with_preset(preset).unwrap();
        let truth_file = temp_files.new_file("sort-preset");
        sort_gaf(&input_file, &truth_file, &params).expect("Failed to sort GAF file for truth");
        let truth = fs::read(&truth_file).expect("Failed to read truth file");
        let output_file = temp_files.new_file("sort-preset");
        run_gaf_base_sort_with_preset(&input_file, &output_file, preset, &[]);
        let output = fs::read(&output_file).expect("Failed to read output file");
        let correct_output = output == truth;
        assert!(correct_output, "gaf-base sort produced incorrect output with preset {}", preset);

        // Then we try overriding something in the preset.
        params.records_per_file /= 2;
        let truth_file = temp_files.new_file("sort-preset-override");
        sort_gaf(&input_file, &truth_file, &params).expect("Failed to sort GAF file for truth");
        let truth = fs::read(&truth_file).expect("Failed to read truth file");
        let output_file = temp_files.new_file("sort-preset-override");
        let other_args = vec![String::from("--records-per-file"), params.records_per_file.to_string()];
        run_gaf_base_sort_with_preset(&input_file, &output_file, preset, &other_args);
        let output = fs::read(&output_file).expect("Failed to read output file");
        let correct_output = output == truth;
        assert!(correct_output, "gaf-base sort produced incorrect output with preset {} and override", preset);
    }
}

//-----------------------------------------------------------------------------

// Tests for `gaf-base construct`.
// NOTE: micb-kir3dl1_HG003.gaf is already sorted with `vg gamsort`, which uses a slightly
// different sort key than `gaf-base sort`. micb-kir3dl1_HG003.gbwt uses the same order.

// We assume that the sorted input file is a temporary file if no output file is provided.
// Otherwise there may be conflicts with multiple tests running in parallel.
fn run_gaf_base_construct(
    temp_files: &mut TempFileHandler,
    sorted_input: &PathBuf, gbwt_file: Option<&PathBuf>, graph_file: Option<&PathBuf>, output_file: Option<&PathBuf>,
    params: &GAFBaseParams, overwrite: bool
) -> (PathBuf, ExitStatus) {
    let mut args = vec![String::from("construct")];
    let output_file = match output_file {
        Some(path) => {
            args.push(String::from("--output"));
            args.push(path.to_str().unwrap().to_string());
            path.clone()
        },
        None => {
            let path = sorted_input.with_added_extension("db");
            temp_files.add_file(&path);
            path
        },
    };
    if let Some(gbwt_file) = gbwt_file {
        args.push(String::from("--gbwt"));
        args.push(gbwt_file.to_str().unwrap().to_string());
    }
    if let Some(graph_file) = graph_file {
        args.push(String::from("--ref-free"));
        args.push(graph_file.to_str().unwrap().to_string());
    }
    args.push(String::from("--block-size"));
    args.push(params.block_size.to_string());
    if !params.store_quality_strings {
        args.push(String::from("--no-quality"));
    }
    if !params.store_optional_fields {
        args.push(String::from("--no-optional"));
    }
    if overwrite {
        args.push(String::from("--overwrite"));
    }
    args.push(sorted_input.to_str().unwrap().to_string());

    let binary = get_binary_path("gaf-base");
    let result = Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gaf-base construct");
    (output_file, result.status)
}

fn run_gaf_base_construct_with_preset(
    sorted_input: &PathBuf, output_file: &PathBuf,
    preset: &str, other_args: &[String]
) {
    let mut args = vec![String::from("construct")];
    args.push(String::from("--preset"));
    args.push(preset.to_string());
    args.push(String::from("--output"));
    args.push(output_file.to_str().unwrap().to_string());
    args.extend_from_slice(other_args);
    args.push(sorted_input.to_str().unwrap().to_string());

    let binary = get_binary_path("gaf-base");
    let result = Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gaf-base construct");
    assert!(result.status.success(), "gaf-base construct with preset {} and args {:?} failed with status: {}", preset, other_args, result.status);
}

#[test]
fn gaf_base_construct() {
    // In this test, we also check that we can build GAF-base from the output of
    // `gaf-base sort`. Other tests start from a pre-sorted GAF file.
    let mut temp_files = TempFileHandler::new();
    let original_input = utils::get_test_data("shuffled.gaf");
    let sorted_input = sort_gaf_to_temp(&mut temp_files, &original_input);
    let gbwt_file = None;
    let graph_file = None;
    let params = GAFBaseParams::default();

    // We expect that the GAF-base built with `gaf-base construct` is identical to
    // one built with a direct library call.
    let truth_file = build_gaf_base_truth(
        &mut temp_files, &sorted_input, None, None, &GAFBaseParams::default()
    );
    let (output_file, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, None, &params, false
    );
    assert!(status.success(), "gaf-base construct failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
    let output_metadata = fs::metadata(&output_file).expect("Failed to read output file metadata");
    let old_timestamp = output_metadata.modified().expect("Failed to get output file timestamp");

    // Try rebuilding the database without overwrite, expect failure.
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, None, &params, false
    );
    assert!(!status.success(), "gaf-base construct should have failed without overwrite");

    // Rebuild with overwrite, expect success.
    let (output_file, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, None, &params, true
    );
    assert!(status.success(), "gaf-base construct failed with overwrite with status: {}", status);

    // Check that we did actually overwrite the output file.
    let output_metadata = fs::metadata(&output_file).expect("Failed to read output file metadata");
    let new_timestamp = output_metadata.modified().expect("Failed to get output file timestamp");
    assert!(new_timestamp > old_timestamp, "Output file timestamp did not change after overwrite");
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file after overwrite");

    // Check that the default parameters correspond to the default arguments.
    let input_copy = temp_files.create_copy(&sorted_input);
    let default_args = vec![String::from("construct"), input_copy.to_str().unwrap().to_string()];
    let binary = get_binary_path("gaf-base");
    let default_output = Command::new(&binary)
        .args(&default_args)
        .output()
        .expect("Failed to execute gaf-base construct with default args");
    assert!(default_output.status.success(), "gaf-base construct failed with default args");
    let correct_output = compare_files(&output_file, &truth_file);
    assert!(correct_output, "gaf-base construct produced incorrect output with default args");

    // Try specifying an explicit output file name.
    let explicit_output_file = temp_files.new_file("gaf-base-output");
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, Some(&explicit_output_file), &params, false
    );
    assert!(status.success(), "gaf-base construct failed with explicit output file with status: {}", status);
    assert!(compare_files(&explicit_output_file, &truth_file), "Explicit output file does not match truth file");
}

#[test]
fn gaf_base_construct_ref_free() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let params = GAFBaseParams { reference_free: true, ..GAFBaseParams::default() };

    let truth_file = build_gaf_base_truth(
        &mut temp_files, &sorted_input, gbwt_file, Some(&graph_file), &params
    );
    let output_file = temp_files.new_file("gaf-base");
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, Some(&graph_file), Some(&output_file), &params, false
    );
    assert!(status.success(), "gaf-base construct ref-free failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
}

#[test]
fn gaf_base_construct_with_gbwt() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
    let graph_file = None;
    let params = GAFBaseParams::default();

    let truth_file = build_gaf_base_truth(
        &mut temp_files, &sorted_input, Some(&gbwt_file), graph_file, &params
    );
    let output_file = temp_files.new_file("gaf-base");
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, Some(&gbwt_file), graph_file, Some(&output_file), &params, false
    );
    assert!(status.success(), "gaf-base construct with GBWT failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
}

#[test]
fn gaf_base_construct_no_quality() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = None;
    let params = GAFBaseParams { store_quality_strings: false, ..GAFBaseParams::default() };

    let truth_file = build_gaf_base_truth(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, &params
    );
    let output_file = temp_files.new_file("gaf-base");
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, Some(&output_file), &params, false
    );
    assert!(status.success(), "gaf-base construct no quality failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
}

#[test]
fn gaf_base_construct_no_optional() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = None;
    let params = GAFBaseParams { store_optional_fields: false, ..GAFBaseParams::default() };

    let truth_file = build_gaf_base_truth(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, &params
    );
    let output_file = temp_files.new_file("gaf-base");
    let (_, status) = run_gaf_base_construct(
        &mut temp_files, &sorted_input, gbwt_file, graph_file, Some(&output_file), &params, false
    );
    assert!(status.success(), "gaf-base construct no optional failed with status: {}", status);
    assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file");
}

#[test]
fn gaf_base_construct_block_sizes() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = None;
    let block_sizes = vec![10, 100, 1000];

    for block_size in block_sizes {
        let params = GAFBaseParams { block_size, ..GAFBaseParams::default() };
        let truth_file = build_gaf_base_truth(
            &mut temp_files, &sorted_input, gbwt_file, graph_file, &params
        );
        let output_file = temp_files.new_file("gaf-base");
        let (_, status) = run_gaf_base_construct(
            &mut temp_files, &sorted_input, gbwt_file, graph_file, Some(&output_file), &params, false
        );
        assert!(status.success(), "gaf-base construct with block size {} failed with status: {}", block_size, status);
        assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file for block size {}", block_size);
    }
}

#[test]
fn gaf_base_construct_presets() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = None;

    for preset in GAFBaseParams::PRESETS {
        // First we use just the preset itself.
        let mut params = GAFBaseParams::with_preset(preset).unwrap();
        let truth_file = build_gaf_base_truth(
            &mut temp_files, &sorted_input, gbwt_file, graph_file, &params
        );
        let output_file = temp_files.new_file("gaf-base");
        run_gaf_base_construct_with_preset(&sorted_input, &output_file, preset, &[]);
        assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file for preset {}", preset);

        // Then we try overriding something in the preset.
        params.block_size *= 2;
        let truth_file = build_gaf_base_truth(
            &mut temp_files, &sorted_input, gbwt_file, graph_file, &params
        );
        let output_file = temp_files.new_file("gaf-base");
        let other_args = vec![String::from("--block-size"), params.block_size.to_string()];
        run_gaf_base_construct_with_preset(&sorted_input, &output_file, preset, &other_args);
        assert!(compare_files(&output_file, &truth_file), "Output file does not match truth file for preset {} with override", preset);
    }
}

//-----------------------------------------------------------------------------

// Tests for `gaf-base decompress`.
// There is no corresponding library function, and the output is not guaranteed to be
// identical to the input. (For example, the representation is normalized and the order
// of optional fields may change).

fn run_gaf_base_decompress(
    input_file: &PathBuf, output_file: Option<&PathBuf>, graph_file: Option<&PathBuf>,
    chunk_size: Option<usize>
) -> Output {
    let mut args = vec![String::from("decompress")];
    if let Some(output_file) = output_file {
        args.push(String::from("--output"));
        args.push(output_file.to_str().unwrap().to_string());
    }
    if let Some(graph_file) = graph_file {
        args.push(String::from("--reference"));
        args.push(graph_file.to_str().unwrap().to_string());
    }
    if let Some(chunk_size) = chunk_size {
        args.push(String::from("--chunk-size"));
        args.push(chunk_size.to_string());
    }
    args.push(input_file.to_str().unwrap().to_string());

    let binary = get_binary_path("gaf-base");
    Command::new(binary)
        .args(&args)
        .output()
        .expect("Failed to execute gaf-base decompress")
}

#[test]
fn gaf_base_decompress() {
    let mut temp_files = TempFileHandler::new();
    let sorted_input = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gbwt_file = None;
    let graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let params = GAFBaseParams::default();
    let expected_lines = count_lines_in_file(&sorted_input);
    let input_file = build_gaf_base_truth(&mut temp_files, &sorted_input, gbwt_file, None, &params);
    let output_file = None;
    let chunk_size = None;

    // Default arguments, store the output as the expected baseline.
    let baseline;
    {
        let output = run_gaf_base_decompress(&input_file, output_file, Some(&graph_file), chunk_size);
        assert!(output.status.success(), "gaf-base decompress failed with status: {}", output.status);
        let output_lines = count_lines(&output.stdout);
        assert_eq!(output_lines, expected_lines, "Decompressed output has incorrect number of lines");
        baseline = output.stdout;
    }

    // Try different chunk sizes.
    for chunk_size in [1, 10, 100] {
        let output = run_gaf_base_decompress(&input_file, output_file, Some(&graph_file), Some(chunk_size));
        assert!(output.status.success(), "gaf-base decompress with chunk size {} failed with status: {}", chunk_size, output.status);
        let correct_output = output.stdout == baseline;
        assert!(correct_output, "gaf-base decompress produced unexpected output with chunk size {}", chunk_size);
    }

    // Specify an explicit output file.
    {
        let output_file = temp_files.new_file("decompressed");
        let output = run_gaf_base_decompress(&input_file, Some(&output_file), Some(&graph_file), chunk_size);
        assert!(output.status.success(), "gaf-base decompress with file output failed with status: {}", output.status);
        let output_data = fs::read(&output_file).expect("Failed to read output file");
        let correct_output = output_data == baseline;
        assert!(correct_output, "gaf-base decompress produced unexpected output in file");
    }

    // Build and decompress a reference-free GAF-base.
    {
        let params = GAFBaseParams { reference_free: true, ..GAFBaseParams::default() };
        let input_file = build_gaf_base_truth(&mut temp_files, &sorted_input, gbwt_file, Some(&graph_file), &params);
        let output = run_gaf_base_decompress(&input_file, output_file, None, chunk_size);
        assert!(output.status.success(), "reference-free gaf-base decompress failed with status: {}", output.status);
        let correct_output = output.stdout == baseline;
        assert!(correct_output, "gaf-base decompress produced unexpected output for reference-free GAF-base");
    }
}

//-----------------------------------------------------------------------------

// Tests for `gbz-base query`.

fn parse_gfa(gfa_file: &PathBuf) -> (Vec<usize>, Vec<(FullPathName, Range<usize>)>) {
    let gfa_data = fs::read(gfa_file).expect("Failed to read GFA file");
    let gfa_str = String::from_utf8_lossy(&gfa_data);

    let mut ref_samples = HashSet::new();
    let mut node_ids = Vec::new();
    let mut reference_paths = Vec::new();
    for gfa_line in gfa_str.lines() {
        if gfa_line.starts_with("H") {
            let fields = gfa_line.split('\t').collect::<Vec<_>>();
            for field in &fields[1..] {
                let tag_parts = field.split(':').collect::<Vec<_>>();
                if tag_parts.len() == 3 && tag_parts[0] == "RS" && tag_parts[1] == "Z" {
                    for sample in tag_parts[2].split(' ') {
                        ref_samples.insert(sample.to_string());
                    }
                }
            }
        } else if gfa_line.starts_with('S') {
            let fields = gfa_line.split('\t').collect::<Vec<_>>();
            if fields.len() >= 2 {
                if let Ok(node_id) = fields[1].parse::<usize>() {
                    node_ids.push(node_id);
                }
            }
        } else if gfa_line.starts_with('W') {
            let fields = gfa_line.split('\t').collect::<Vec<_>>();
            if fields.len() >= 6 && ref_samples.contains(fields[1]) {
                let path_name = FullPathName::reference(fields[1], fields[3]);
                let start = fields[4].parse::<usize>().expect("Failed to parse start position");
                let end = fields[5].parse::<usize>().expect("Failed to parse end position");
                reference_paths.push((path_name, start..end));
            }
        }
    }

    (node_ids, reference_paths)
}

fn random_reference_offset(reference_paths: &[(FullPathName, Range<usize>)], rng: &mut impl Rng) -> (SubgraphQuery, Vec<String>) {
    let (path_name, range) = reference_paths.choose(rng).expect("No reference paths available");
    let offset = rng.random_range(range.clone());
    let args = vec![
        String::from("--sample"),
        path_name.sample.clone(),
        String::from("--contig"),
        path_name.contig.clone(),
        String::from("--offset"),
        offset.to_string(),
    ];
    (SubgraphQuery::path_offset(path_name, offset), args)
}

fn random_reference_interval(reference_paths: &[(FullPathName, Range<usize>)], rng: &mut impl Rng) -> (SubgraphQuery, Vec<String>) {
    let (path_name, range) = reference_paths.choose(rng).expect("No reference paths available");
    let len = rng.random_range(10..300);
    let offset = rng.random_range(range.start..(range.end - len));
    let args = vec![
        String::from("--sample"),
        path_name.sample.clone(),
        String::from("--contig"),
        path_name.contig.clone(),
        String::from("--interval"),
        format!("{}..{}", offset, offset + len),
    ];
    (SubgraphQuery::path_interval(path_name, offset..(offset + len)), args)
}

fn random_node(node_ids: &[usize], rng: &mut impl Rng) -> (SubgraphQuery, Vec<String>) {
    let node_id = *node_ids.choose(rng).expect("No node IDs available");
    (SubgraphQuery::nodes([node_id]), vec![String::from("--node"), node_id.to_string()])
}

fn random_node_with_neighbors(node_ids: &[usize], rng: &mut impl Rng) -> (SubgraphQuery, Vec<String>) {
    let initial_offset = rng.random_range(0..node_ids.len());
    let start = if initial_offset > 0 { initial_offset - 1 } else { 0 };
    let end = if initial_offset + 1 <= node_ids.len() { initial_offset + 1 } else { node_ids.len() };
    let mut args = Vec::new();
    for node_id in &node_ids[start..end] {
        args.push(String::from("--node"));
        args.push(node_id.to_string());
    }
    (SubgraphQuery::nodes(node_ids[start..end].to_vec()), args)
}

fn append_query_variants(
    queries: &mut Vec<SubgraphQuery>, arg_lists: &mut Vec<Vec<String>>,
    query: SubgraphQuery, args: Vec<String>
) {
    for context in [0, 30, 100] {
        for snarls in [SnarlOutput::None, SnarlOutput::Contained, SnarlOutput::Overlapping] {
            for output in [HaplotypeOutput::All, HaplotypeOutput::Distinct, HaplotypeOutput::ReferenceOnly, HaplotypeOutput::None] {
                if query.is_node_based() && output == HaplotypeOutput::ReferenceOnly {
                    // No reference path.
                    continue;
                }
                if query.is_multi_node() && snarls == SnarlOutput::Overlapping {
                    // Not implemented.
                    continue;
                }
                queries.push(query.clone().with_context(context).with_snarls(snarls).with_haplotypes(output));
                let mut args = args.to_vec();
                args.push(String::from("--context"));
                args.push(context.to_string());
                match snarls {
                    SnarlOutput::None => {},
                    SnarlOutput::Contained => {
                        args.push(String::from("--snarls"));
                    },
                    SnarlOutput::Overlapping => {
                        args.push(String::from("--extend-snarls"));
                    },
                }
                args.push(String::from("--haplotypes"));
                match output {
                    HaplotypeOutput::All => {
                        args.push(String::from("all"));
                    },
                    HaplotypeOutput::Distinct => {
                        args.push(String::from("distinct"));
                    },
                    HaplotypeOutput::ReferenceOnly => {
                        args.push(String::from("reference-only"));
                    },
                    HaplotypeOutput::None => {
                        args.push(String::from("none"));
                    },
                }
                arg_lists.push(args);
            }
        }
    }
}

// This does not generate between queries, because they only make sense with specific nodes.
fn generate_queries(gfa_file: &PathBuf, include_variants: bool) -> (Vec<SubgraphQuery>, Vec<Vec<String>>) {
    const NUM_QUERIES: usize = 3;

    let (node_ids, reference_paths) = parse_gfa(gfa_file);
    let mut queries = Vec::new();
    let mut arg_lists = Vec::new();
    let mut rng = rand::rng();
    for _ in 0..NUM_QUERIES {
        let (query, args) = random_reference_offset(&reference_paths, &mut rng);
        // We also want a version of the query with default params and args.
        queries.push(query.clone());
        arg_lists.push(args.clone());
        if include_variants {
            append_query_variants(&mut queries, &mut arg_lists, query, args);
        }
    }
    for _ in 0..NUM_QUERIES {
        let (query, args) = random_reference_interval(&reference_paths, &mut rng);
        queries.push(query.clone());
        arg_lists.push(args.clone());
        if include_variants {
            append_query_variants(&mut queries, &mut arg_lists, query, args);
        }
    }
    for _ in 0..NUM_QUERIES {
        let (query, args) = random_node(&node_ids, &mut rng);
        queries.push(query.clone());
        arg_lists.push(args.clone());
        if include_variants {
            append_query_variants(&mut queries, &mut arg_lists, query, args);
        }
    }
    for _ in 0..NUM_QUERIES {
        let (query, args) = random_node_with_neighbors(&node_ids, &mut rng);
        queries.push(query.clone());
        arg_lists.push(args.clone());
        if include_variants {
            append_query_variants(&mut queries, &mut arg_lists, query, args);
        }
    }

    (queries, arg_lists)
}

fn run_gbz_base_query(
    graph_file: &PathBuf, chains_file: Option<&PathBuf>,
    query_args: &[String], output_format: Option<&str>, with_cigar: bool
) -> Output {
    let mut args = vec![String::from("query")];
    args.extend_from_slice(query_args);
    if let Some(format) = output_format {
        args.push(String::from("--format"));
        args.push(format.to_string());
    }
    if with_cigar {
        args.push(String::from("--cigar"));
    }
    if let Some(chains_file) = chains_file {
        args.push(String::from("--chains"));
        args.push(chains_file.to_str().unwrap().to_string());
    }
    args.push(graph_file.to_str().unwrap().to_string());

    let binary = get_binary_path("gbz-base");
    Command::new(&binary)
        .args(&args)
        .output()
        .expect("Failed to execute gbz-base query")
}

// Returns an error if the query itself fails.
fn run_subgraph_query_with_gbz(
    graph: &GBZ, path_index: Option<&PathIndex>, chains: Option<&Chains>,
    query: &SubgraphQuery, output_format: &str, with_cigar: bool
) -> Result<Vec<u8>, Error> {
    let mut subgraph = Subgraph::new();
    subgraph.from_gbz(graph, path_index, chains, query)?;

    let mut output = Vec::new();
    match output_format {
        "gfa" => {
            let result = subgraph.write_gfa(&mut output, with_cigar);
            assert!(result.is_ok(), "Failed to write GFA output for query {}: {}", query, result.unwrap_err());
        },
        "json" => {
            let result = subgraph.write_json(&mut output, with_cigar);
            assert!(result.is_ok(), "Failed to write JSON output for query {}: {}", query, result.unwrap_err());
        },
        _ => panic!("Unsupported output format: {}", output_format),
    }

    Ok(output)
}

// Returns an error if the query itself fails.
fn run_subgraph_query_with_db(
    interface: &mut GraphInterface, query: &SubgraphQuery, output_format: &str, with_cigar: bool
) -> Result<Vec<u8>, Error> {
    let mut subgraph = Subgraph::new();
    subgraph.from_db(interface, query)?;

    let mut output = Vec::new();
    match output_format {
        "gfa" => {
            let result = subgraph.write_gfa(&mut output, with_cigar);
            assert!(result.is_ok(), "Failed to write GFA output for query {}: {}", query, result.unwrap_err());
        },
        "json" => {
            let result = subgraph.write_json(&mut output, with_cigar);
            assert!(result.is_ok(), "Failed to write JSON output for query {}: {}", query, result.unwrap_err());
        },
        _ => panic!("Unsupported output format: {}", output_format),
    }

    Ok(output)
}

#[test]
fn gbz_base_query() {
    let mut temp_files = TempFileHandler::new();

    let graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let gbz: GBZ = serialize::load_from(&graph_file).expect("Failed to load GBZ graph");
    let path_index = PathIndex::new(&gbz, GBZBase::INDEX_INTERVAL, false).expect("Failed to build path index");
    let chains_file = utils::get_test_data("micb-kir3dl1.chains");
    let chains: Chains = serialize::load_from(&chains_file).expect("Failed to load chains file");

    let gbz_base_file = build_gbz_base_truth(&mut temp_files, &graph_file, None);
    let database = GBZBase::open(&gbz_base_file).expect("Failed to open GBZ-base database");
    let mut interface = GraphInterface::new(&database).expect("Failed to create GraphInterface");

    let gfa_file = utils::get_test_data("micb-kir3dl1.gfa");
    let (queries, arg_lists) = generate_queries(&gfa_file, true);

    for (query, args) in queries.iter().zip(arg_lists.iter()) {
        for output_format in ["gfa", "json"] {
            for with_cigar in [false, true] {
                let gbz_output = run_gbz_base_query(&graph_file, Some(&chains_file), args, Some(output_format), with_cigar);
                assert!(gbz_output.status.success(), "gbz-base query with GBZ failed for query {} with format {} and CIGAR {}", query, output_format, with_cigar);
                let library_gbz = run_subgraph_query_with_gbz(&gbz, Some(&path_index), Some(&chains), query, output_format, with_cigar);
                assert!(library_gbz.is_ok(), "Subgraph query failed for query {} with format {} and CIGAR {}: {}", query, output_format, with_cigar, library_gbz.unwrap_err());
                let gbz_query_ok = library_gbz.unwrap() == gbz_output.stdout;
                assert!(gbz_query_ok, "Output mismatch for query {} with format {} and CIGAR {}", query, output_format, with_cigar);

                let db_output = run_gbz_base_query(&gbz_base_file, None, args, Some(output_format), with_cigar);
                assert!(db_output.status.success(), "gbz-base query with DB failed for query {} with format {} and CIGAR {}", query, output_format, with_cigar);
                let library_db = run_subgraph_query_with_db(&mut interface, query, output_format, with_cigar);
                assert!(library_db.is_ok(), "Subgraph query failed for query {} with format {} and CIGAR {}: {}", query, output_format, with_cigar, library_db.unwrap_err());
                let db_query_ok = library_db.unwrap() == db_output.stdout;
                assert!(db_query_ok, "Output mismatch for query {} with format {} and CIGAR {}", query, output_format, with_cigar);
            }
        }
    }
}

// Here we also test setting a safety limit.
// We already know from unit tests that the limit will trigger with all kinds of queries,
// so we don't need to repeat that here.
#[test]
fn gbz_base_query_between() {
    let mut temp_files = TempFileHandler::new();

    let graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let gbz: GBZ = serialize::load_from(&graph_file).expect("Failed to load GBZ graph");

    let gbz_base_file = build_gbz_base_truth(&mut temp_files, &graph_file, None);
    let database = GBZBase::open(&gbz_base_file).expect("Failed to open GBZ-base database");
    let mut interface = GraphInterface::new(&database).expect("Failed to create GraphInterface");

    let node_pairs = vec![
        (400, 403),   // Adjacent boundary nodes in the first component.
        (103, 255),   // Distant boundary nodes in the first component.
        (1, 659),     // Start/end of the first component.
        (1981, 1984), // Adjacent boundary nodes in the second component.
        (1169, 1550), // Distant boundary nodes in the second component.
        (660, 2891),  // Start/end of the second component.
    ];
    for (start, end) in node_pairs {
        for limit in [None, Some(10), Some(100), Some(1000)] {
            let query = SubgraphQuery::between(
                support::encode_node(start, Orientation::Forward),
                support::encode_node(end, Orientation::Forward)
            ).with_limit(limit);
            let mut args = vec![
                String::from("--between"),
                format!("{}:{}", start, end),
            ];
            if let Some(l) = limit {
                args.push(String::from("--limit"));
                args.push(l.to_string());
            }

            let gbz_output = run_gbz_base_query(&graph_file, None, &args, None, false);
            let library_gbz = run_subgraph_query_with_gbz(&gbz, None, None, &query, "gfa", false);
            assert_eq!(gbz_output.status.success(), library_gbz.is_ok(), "Query success mismatch with GBZ for query {}", query);
            if let Ok(output) = library_gbz {
                let gbz_query_ok = output == gbz_output.stdout;
                assert!(gbz_query_ok, "Output mismatch for query {}", query);
            }

            let db_output = run_gbz_base_query(&gbz_base_file, None, &args, Some("gfa"), false);
            let library_db = run_subgraph_query_with_db(&mut interface, &query, "gfa", false);
            assert_eq!(db_output.status.success(), library_db.is_ok(), "Query success mismatch with DB for query {}", query);
            if let Ok(output) = library_db {
                let db_query_ok = output == db_output.stdout;
                assert!(db_query_ok, "Output mismatch for query {}", query);
            }
        }
    }
}

// Returns GAF output.
// We assume that the subgraph query itself has been tested earlier.
fn run_gbz_gaf_base_query(
    temp_files: &mut TempFileHandler,
    graph_file: &PathBuf, gaf_base_file: &PathBuf,
    query_args: &[String], alignment_output: Option<AlignmentOutput>,
    gaf_only: bool
) -> Vec<u8> {
    let mut query_args = query_args.to_vec();
    query_args.push(String::from("--gaf-base"));
    query_args.push(gaf_base_file.to_str().unwrap().to_string());
    let gaf_output_file = if gaf_only {
        query_args.push(String::from("--gaf-only"));
        None
    } else {
        query_args.push(String::from("--gaf-output"));
        let gaf_output_file = temp_files.new_file("gaf-base-output");
        query_args.push(gaf_output_file.to_str().unwrap().to_string());
        Some(gaf_output_file)
    };
    match alignment_output {
        Some(output) => {
            query_args.push(String::from("--alignments"));
            query_args.push(output.to_string());
        },
        None => {},
    }

    let output = run_gbz_base_query(graph_file, None, &query_args, None, false);
    assert!(output.status.success(), "gbz-base query with args {:?} failed with status: {}", query_args, output.status);
    if let Some(gaf_output_file) = gaf_output_file {
        fs::read(&gaf_output_file).expect("Failed to read GAF output file")
    } else {
        output.stdout
    }
}

fn write_gaf(read_set: &ReadSet, subgraph: &Subgraph) -> Vec<u8> {
    let mut output = Vec::new();
    formats::write_gaf_file_header(&mut output, None).expect("Failed to write GAF file header");
    if let Some(graph_name) = subgraph.graph_name() {
        let header_lines = graph_name.to_gaf_header_lines();
        formats::write_header_lines(&header_lines, &mut output).expect("Failed to write GAF header lines");
    }
    read_set.to_gaf(&mut output).expect("Failed to write GAF records");
    output
}

// Returns GAF output.
// We assume that the subgraph query itself has been tested earlier.
fn run_subgraph_query_with_gbz_gaf(
    graph: &GBZ, path_index: Option<&PathIndex>, gaf_base: &GAFBase,
    query: &SubgraphQuery, output: AlignmentOutput
) -> Vec<u8> {
    let mut subgraph = Subgraph::new();
    let result = subgraph.from_gbz(graph, path_index, None, query);
    assert!(result.is_ok(), "Failed to build subgraph for query {}: {}", query, result.unwrap_err());
    let read_set = ReadSet::new(GraphReference::Gbz(graph), &subgraph, gaf_base, output);
    assert!(read_set.is_ok(), "Failed to build read set for query {} with alignment output {}: {}", query, read_set.unwrap_err(), output);
    write_gaf(&read_set.unwrap(), &subgraph)
}

// Returns GAF output.
// We assume that the subgraph query itself has been tested earlier.
fn run_subgraph_query_with_db_gaf(
    interface: &mut GraphInterface, gaf_base: &GAFBase,
    query: &SubgraphQuery, output: AlignmentOutput
) -> Vec<u8> {
    let mut subgraph = Subgraph::new();
    let result = subgraph.from_db(interface, query);
    assert!(result.is_ok(), "Failed to build subgraph for query {}: {}", query, result.unwrap_err());
    let read_set = ReadSet::new(GraphReference::Db(interface), &subgraph, gaf_base, output);
    assert!(read_set.is_ok(), "Failed to build read set for query {} with alignment output {}: {}", query, read_set.unwrap_err(), output);
    write_gaf(&read_set.unwrap(), &subgraph)
}

#[test]
fn gbz_base_query_with_gaf_base() {
    let mut temp_files = TempFileHandler::new();

    let graph_file = utils::get_test_data("micb-kir3dl1.gbz");
    let gbz: GBZ = serialize::load_from(&graph_file).expect("Failed to load GBZ graph");
    let path_index = PathIndex::new(&gbz, GBZBase::INDEX_INTERVAL, false).expect("Failed to build path index");

    let gbz_base_file = build_gbz_base_truth(&mut temp_files, &graph_file, None);
    let database = GBZBase::open(&gbz_base_file).expect("Failed to open GBZ-base database");
    let mut interface = GraphInterface::new(&database).expect("Failed to create GraphInterface");

    let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let gaf_base_file = build_gaf_base_truth(&mut temp_files, &gaf_file, None, None, &GAFBaseParams::default());
    let gaf_base = GAFBase::open(&gaf_base_file).expect("Failed to open GAF-base database");

    let gfa_file = utils::get_test_data("micb-kir3dl1.gfa");
    let (queries, arg_lists) = generate_queries(&gfa_file, false);

    for (query, args) in queries.iter().zip(arg_lists.iter()) {
        for alignment_output in [None, Some(AlignmentOutput::Overlapping), Some(AlignmentOutput::Clipped), Some(AlignmentOutput::Contained)] {

            let binary_gbz_gaf = run_gbz_gaf_base_query(&mut temp_files, &graph_file, &gaf_base_file, args, alignment_output, false);
            let output = alignment_output.unwrap_or(AlignmentOutput::Clipped);
            let library_gbz_gaf = run_subgraph_query_with_gbz_gaf(&gbz, Some(&path_index), &gaf_base, &query, output);
            let gbz_gaf_query_ok = binary_gbz_gaf == library_gbz_gaf;
            assert!(gbz_gaf_query_ok, "Output mismatch for query {} with alignment output {}", query, output);

            if alignment_output.is_none() {
                let gaf_only = run_gbz_gaf_base_query(&mut temp_files, &graph_file, &gaf_base_file, args, alignment_output, true);
                let gaf_only_ok = gaf_only == library_gbz_gaf;
                assert!(gaf_only_ok, "Output mismatch for query {} with --gaf-only", query);
            }

            let binary_db_gaf = run_gbz_gaf_base_query(&mut temp_files, &gbz_base_file, &gaf_base_file, args, alignment_output, false);
            let library_db_gaf = run_subgraph_query_with_db_gaf(&mut interface, &gaf_base, &query, output);
            let db_gaf_query_ok = binary_db_gaf == library_db_gaf;
            assert!(db_gaf_query_ok, "Output mismatch for query {} with alignment output {}", query, output);
        }
    }
}

//-----------------------------------------------------------------------------
