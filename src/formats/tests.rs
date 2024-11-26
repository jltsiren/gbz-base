use super::*;

use crate::{db, formats};

use gbwt::Metadata;

use simple_sds::ops::BitVec;
use simple_sds::serialize;

use std::io::BufRead;

//-----------------------------------------------------------------------------

// Tests for `WalkMetadata`.

#[test]
fn walk_metadata_interval() {
    let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
    let interval = 123..456;
    let metadata = WalkMetadata::path_interval(&path_name, interval.clone());

    let interval_name = FullPathName {
        sample: path_name.sample.clone(),
        contig: path_name.contig.clone(),
        haplotype: path_name.haplotype,
        fragment: path_name.fragment + interval.start
    };
    assert_eq!(metadata.name, interval_name, "Wrong path name");
    assert_eq!(metadata.end, path_name.fragment + interval.end, "Wrong end offset");
    assert!(metadata.weight.is_none(), "Weight should not be set");
    assert!(metadata.cigar.is_none(), "CIGAR string should not be set");
}

#[test]
fn walk_metadata_haplotype() {
    let filename = support::get_test_data("example.meta");
    let metadata: Metadata = serialize::load_from(&filename).unwrap();

    let len = 123;
    for path_id in 0..metadata.paths() {
        let haplotype = WalkMetadata::haplotype(&metadata, path_id, len);
        assert!(haplotype.is_some(), "Failed to create WalkMetadata for path {}", path_id);
        let haplotype = haplotype.unwrap();
        let full_name = FullPathName::from_metadata(&metadata, path_id).unwrap();
        assert_eq!(haplotype.name, full_name, "Wrong name for path {}", path_id);
        assert_eq!(haplotype.end, full_name.fragment + len, "Wrong end offset for path {}", path_id);
        assert!(haplotype.weight.is_none(), "Weight should not be set for path {}", path_id);
        assert!(haplotype.cigar.is_none(), "CIGAR string should not be set for path {}", path_id);
    }
}

#[test]
fn walk_metadata_anonymous() {
    let haplotype = 2;
    let contig = "chr19";
    let len = 456;
    let metadata = WalkMetadata::anonymous(haplotype, contig, len);

    let path_name = FullPathName {
        sample: String::from("unknown"),
        contig: String::from(contig),
        haplotype: haplotype,
        fragment: 0
    };
    assert_eq!(metadata.name, path_name, "Wrong path name");
    assert_eq!(metadata.end, path_name.fragment + len, "Wrong end offset");
    assert!(metadata.weight.is_none(), "Weight should not be set");
    assert!(metadata.cigar.is_none(), "CIGAR string should not be set");
}

#[test]
fn walk_metadata_weight_cigar() {
    let mut metadata = WalkMetadata::anonymous(2, "chr19", 456);
    assert!(metadata.weight.is_none(), "Weight should not be set");
    assert!(metadata.cigar.is_none(), "CIGAR string should not be set");

    metadata.add_weight(Some(42));
    assert_eq!(metadata.weight, Some(42), "Weight was not set correctly");
    metadata.add_weight(None);
    assert!(metadata.weight.is_none(), "Weight was not cleared");

    let cigar = "10M2D5M";
    metadata.add_cigar(Some(String::from(cigar)));
    assert_eq!(metadata.cigar, Some(String::from(cigar)), "CIGAR string was not set correctly");
    metadata.add_cigar(None);
    assert!(metadata.cigar.is_none(), "CIGAR string was not cleared");
}

//-----------------------------------------------------------------------------

// Tests for GFA writing.

#[test]
fn gfa_header() {
    {
        // No reference samples.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_header(None, &mut output);
        assert!(result.is_ok(), "Failed to write GFA header with no reference samples");
        let correct = b"H\tVN:Z:1.1\n";
        assert_eq!(&output, correct, "Wrong GFA header with no reference samples");
    }
    {
        // One reference sample.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_header(Some("GRCh38"), &mut output);
        assert!(result.is_ok(), "Failed to write GFA header with one reference sample");
        let correct = b"H\tVN:Z:1.1\tRS:Z:GRCh38\n";
        assert_eq!(&output, correct, "Wrong GFA header with one reference sample");
    }
    {
        // Two reference samples.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_header(Some("GRCh38 CHM13"), &mut output);
        assert!(result.is_ok(), "Failed to write GFA header with two reference samples");
        let correct = b"H\tVN:Z:1.1\tRS:Z:GRCh38 CHM13\n";
        assert_eq!(&output, correct, "Wrong GFA header with two reference samples");
    }
}

#[test]
fn gfa_segment() {
    {
        // Node with an integer identifier.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_node(42, b"GATTACA", &mut output);
        assert!(result.is_ok(), "Failed to write a GFA segment line for a node");
        let correct = b"S\t42\tGATTACA\n";
        assert_eq!(&output, correct, "Wrong GFA segment line for a node");
    }
    {
        // Reverse orientation.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_segment(b"42", b"GATTACA", &mut output);
        assert!(result.is_ok(), "Failed to write a GFA segment line for a segment");
        let correct = b"S\t42\tGATTACA\n";
        assert_eq!(&output, correct, "Wrong GFA segment line for a segment");
    }
}

#[test]
fn gfa_link() {
    {
        // Forward, reverse with node identifiers.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_edge(
            (42, Orientation::Forward),
            (43, Orientation::Reverse),
            &mut output
        );
        assert!(result.is_ok(), "Failed to write a GFA link line for an edge");
        let correct = b"L\t42\t+\t43\t-\t0M\n";
        assert_eq!(&output, correct, "Wrong GFA link line for an edge");
    }
    {
        // Reverse, forward with segment names.
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_link(
            (b"42", Orientation::Reverse),
            (b"43", Orientation::Forward),
            &mut output
        );
        assert!(result.is_ok(), "Failed to write a GFA link line for a link");
        let correct = b"L\t42\t-\t43\t+\t0M\n";
        assert_eq!(&output, correct, "Wrong GFA link line for a link");
    }
}

#[test]
fn gfa_walk() {
    {
        // No weight and no CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let metadata = WalkMetadata::path_interval(&path_name, interval);
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_walk(&path, &metadata, &mut output);
        assert!(result.is_ok(), "Failed to write a GFA walk line with no weight and no CIGAR string");
        let correct = b"W\tGRCh38\t0\tchr13\t1000\t1123\t>21<22>23\n";
        assert_eq!(&output, correct, "Wrong GFA walk line with no weight and no CIGAR string");
    }
    {
        // With weight but with no CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_weight(Some(42));
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_walk(&path, &metadata, &mut output);
        assert!(result.is_ok(), "Failed to write a GFA walk line with weight but no CIGAR string");
        let correct = b"W\tGRCh38\t0\tchr13\t1000\t1123\t>21<22>23\tWT:i:42\n";
        assert_eq!(&output, correct, "Wrong GFA walk line with weight but no CIGAR string");
    }
    {
        // With no weight but with a CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_walk(&path, &metadata, &mut output);
        assert!(result.is_ok(), "Failed to write a GFA walk line with no weight but with a CIGAR string");
        let correct = b"W\tGRCh38\t0\tchr13\t1000\t1123\t>21<22>23\tCG:Z:10M2D5M\n";
        assert_eq!(&output, correct, "Wrong GFA walk line with no weight but with a CIGAR string");
    }
    {
        // With weight and a CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_weight(Some(42));
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_walk(&path, &metadata, &mut output);
        assert!(result.is_ok(), "Failed to write a GFA walk line with weight and a CIGAR string");
        let correct = b"W\tGRCh38\t0\tchr13\t1000\t1123\t>21<22>23\tWT:i:42\tCG:Z:10M2D5M\n";
        assert_eq!(&output, correct, "Wrong GFA walk line with weight and a CIGAR string");
    }
}

//-----------------------------------------------------------------------------

// Tests for JSON writing.
// TODO: More tests?

#[test]
fn write_json_path() {
    {
        // No weight and no CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let metadata = WalkMetadata::path_interval(&path_name, interval);
        let json_value = json_path(&path, &metadata);
        let result = json_value.to_string();
        let correct = "{\"name\": \"GRCh38#0#chr13[1000-1123]\", \"path\": [{\"id\": \"21\", \"is_reverse\": false}, {\"id\": \"22\", \"is_reverse\": true}, {\"id\": \"23\", \"is_reverse\": false}]}";
        assert_eq!(result, correct, "Wrong JSON for a path with no weight and no CIGAR string");
    }
    {
        // With weight but with no CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_weight(Some(42));
        let json_value = json_path(&path, &metadata);
        let result = json_value.to_string();
        let correct = "{\"name\": \"GRCh38#0#chr13[1000-1123]\", \"weight\": 42, \"path\": [{\"id\": \"21\", \"is_reverse\": false}, {\"id\": \"22\", \"is_reverse\": true}, {\"id\": \"23\", \"is_reverse\": false}]}";
        assert_eq!(result, correct, "Wrong JSON for a path with weight but no CIGAR string");
    }
    {
        // With no weight but with a CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let json_value = json_path(&path, &metadata);
        let result = json_value.to_string();
        let correct = "{\"name\": \"GRCh38#0#chr13[1000-1123]\", \"cigar\": \"10M2D5M\", \"path\": [{\"id\": \"21\", \"is_reverse\": false}, {\"id\": \"22\", \"is_reverse\": true}, {\"id\": \"23\", \"is_reverse\": false}]}";
        assert_eq!(result, correct, "Wrong JSON for a path with no weight but with a CIGAR string");
    }
    {
        // With weight and a CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let mut metadata = WalkMetadata::path_interval(&path_name, interval);
        metadata.add_weight(Some(42));
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let json_value = json_path(&path, &metadata);
        let result = json_value.to_string();
        let correct = "{\"name\": \"GRCh38#0#chr13[1000-1123]\", \"weight\": 42, \"cigar\": \"10M2D5M\", \"path\": [{\"id\": \"21\", \"is_reverse\": false}, {\"id\": \"22\", \"is_reverse\": true}, {\"id\": \"23\", \"is_reverse\": false}]}";
        assert_eq!(result, correct, "Wrong JSON for a path with weight and a CIGAR string");
    }
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: parsing.

// Creates an alignment object for an unaligned sequence with no optional fields.
fn empty_alignment(name: &str, seq_len: usize) -> Alignment {
    let name = String::from(name);
    let seq_interval = 0..0;
    let path = TargetPath::Path(Vec::new());
    let path_len = 0;
    let path_interval = 0..0;
    let matches = 0;
    let edits = 0;
    let mapq = None;
    let score = None;
    let base_quality = None;
    let difference = None;
    let pair = None;
    let optional = Vec::new();
    Alignment {
        name, seq_len, seq_interval,
        path, path_len, path_interval,
        matches, edits, mapq, score,
        base_quality, difference, pair,
        optional
    }
}

// Parses the alignments from the test file.
// Returns them as a vector, unless told to expect a parse error.
fn parse_alignments(filename: &PathBuf, expect_error: bool) -> Vec<Alignment> {
    let reader = db::open_file(filename);
    assert!(reader.is_ok(), "Failed to open the test file: {}", reader.err().unwrap());
    let mut reader = reader.unwrap();

    let mut result = Vec::new();
    let mut line_num = 1;
    loop {
        let mut buf: Vec<u8> = Vec::new();
        let len = reader.read_until(b'\n', &mut buf);
        assert!(len.is_ok(), "Failed to read line {}: {}", line_num, len.err().unwrap());
        if len.unwrap() == 0 {
            break;
        }
        let alignment = Alignment::from_gaf(&buf);
        if expect_error {
            assert!(alignment.is_err(), "Expected a parse error on line {}", line_num);
        } else {
            assert!(alignment.is_ok(), "Failed to parse line {}: {}", line_num, alignment.err().unwrap());
            result.push(alignment.unwrap());
        }
        line_num += 1;
    }

    result
}

// Check that the alignment matches the truth.
// May be instructed to skip fields storing relative information and/or unknown optional fields.
fn check_alignment(alignment: &Alignment, truth: &Alignment, line: usize, skip_relative: bool, skip_optional: bool) {
    if !skip_relative {
        assert_eq!(alignment.name, truth.name, "Wrong sequence name on line {}", line);
    }
    assert_eq!(alignment.seq_len, truth.seq_len, "Wrong sequence length on line {}", line);
    assert_eq!(alignment.seq_interval, truth.seq_interval, "Wrong sequence interval on line {}", line);

    if !skip_relative {
        assert_eq!(alignment.path, truth.path, "Wrong path on line {}", line);
    }
    assert_eq!(alignment.path_len, truth.path_len, "Wrong path length on line {}", line);
    assert_eq!(alignment.path_interval, truth.path_interval, "Wrong path interval on line {}", line);

    assert_eq!(alignment.matches, truth.matches, "Wrong number of matches on line {}", line);
    assert_eq!(alignment.edits, truth.edits, "Wrong number of edits on line {}", line);
    assert_eq!(alignment.mapq, truth.mapq, "Wrong mapping quality on line {}", line);
    assert_eq!(alignment.score, truth.score, "Wrong alignment score on line {}", line);

    assert_eq!(alignment.base_quality, truth.base_quality, "Wrong base quality string on line {}", line);
    assert_eq!(alignment.difference, truth.difference, "Wrong difference string on line {}", line);
    assert_eq!(alignment.pair, truth.pair, "Wrong pair information on line {}", line);

    if !skip_optional {
        assert_eq!(alignment.optional, truth.optional, "Wrong optional fields on line {}", line);
    }
}

#[test]
fn alignment_known_good() {
    // Construct the correct alignment object.
    let name = String::from("forward");
    let seq_len = 100;
    let seq_interval = 5..90;
    let path = TargetPath::Path(vec![
        support::encode_node(10, Orientation::Forward),
        support::encode_node(12, Orientation::Forward),
        support::encode_node(14, Orientation::Reverse),
        support::encode_node(15, Orientation::Forward)
    ]);
    let path_len = 120;
    let path_interval = 10..100;
    let matches = 80;
    let edits = 13;
    let mapq = Some(60);
    let score = Some(57);
    let base_quality = Some(vec![b'?'; 100]);
    let difference = Some(vec![
        Difference::Match(20),
        Difference::Mismatch(b'C'), Difference::Mismatch(b'T'),
        Difference::Match(20),
        Difference::Deletion(8),
        Difference::Match(20),
        Difference::Insertion(b"CAT".to_vec()),
        Difference::Match(20)
    ]);
    let pair = None;
    let optional = vec![
        TypedField::String([b'x', b'x'], b"unknown".to_vec()),
    ];

    // The file contains:
    // * three variants of the same alignment
    // * three additional variants that are paired
    // * two variants of unaligned sequences
    let forward = Alignment {
        name, seq_len, seq_interval,
        path, path_len, path_interval,
        matches, edits, mapq, score,
        base_quality, difference, pair,
        optional
    };
    let mut reverse = forward.clone();
    reverse.name = String::from("reverse");
    if let TargetPath::Path(path) = &mut reverse.path {
        support::reverse_path_in_place(path);
    }
    let mut no_mapq = forward.clone();
    no_mapq.name = String::from("no_mapq");
    no_mapq.mapq = None;

    let mut first = forward.clone();
    first.name = String::from("first");
    first.pair = Some(PairedRead {
        name: String::from("second"),
        is_next: true,
        is_proper: true,
    });
    let mut second = forward.clone();
    second.name = String::from("second");
    second.pair = Some(PairedRead {
        name: String::from("first"),
        is_next: false,
        is_proper: true,
    });
    let mut same = forward.clone();
    same.name = String::from("same");
    same.pair = Some(PairedRead {
        name: String::from("same"),
        is_next: true,
        is_proper: true,
    });

    let empty = empty_alignment("empty", seq_len);
    let missing_values = empty_alignment("missing_values", seq_len);

    // And now the actual test.
    let truth = vec![forward, reverse, no_mapq, first, second, same, empty, missing_values];
    let filename = formats::get_test_data("good.gaf");
    let alignments = parse_alignments(&filename, false);
    assert_eq!(alignments.len(), truth.len(), "Wrong number of alignments in the test file");
    for i in 0..truth.len() {
        check_alignment(&alignments[i], &truth[i], i + 1, false, false);
    }
}

#[test]
fn alignment_known_bad() {
    let filename = formats::get_test_data("bad.gaf");
    let alignments = parse_alignments(&filename, true);
    assert_eq!(alignments.len(), 0, "There should be no valid alignments in the test file");
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: support functions.

#[test]
fn alignment_paths_by_sample() {
    // This is actually a GBWT of haplotypes.
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let metadata = index.metadata().unwrap();

    let true_paths: Vec<u32> = vec![0, 1, 2, 3, 4, 5];
    let true_offsets: Vec<usize> = vec![0, 2, 6];

    let result = Alignment::paths_by_sample(&metadata);
    assert!(result.is_ok(), "Failed to get a list of paths by sample: {}", result.err().unwrap());
    let (paths, offsets) = result.unwrap();

    assert_eq!(paths, true_paths, "Wrong path identifiers");
    assert_eq!(offsets.len(), paths.len() + 1, "Wrong universe size for the index");
    let offsets: Vec<usize> = offsets.one_iter().map(|(_, x)| x).collect();
    assert_eq!(offsets, true_offsets, "Wrong offsets for path id intervals by sample");
}

#[test]
fn alignment_set_relative_information() {
    let filename = support::get_test_data("example.gbwt");
    let index: GBWT = serialize::load_from(&filename).unwrap();
    let metadata = index.metadata().unwrap();

    let result = Alignment::paths_by_sample(&metadata);
    assert!(result.is_ok(), "Failed to get a list of paths by sample: {}", result.err().unwrap());
    let paths_by_sample = result.unwrap();
    let mut used_paths = RawVector::with_len(metadata.paths(), false);

    // The path we are interested in is id 3: (sample, A, 2, 0) with sample id 1.
    let name = String::from("sample");
    let seq_len = 5;
    let seq_interval = 0..5;
    let path = TargetPath::Path(vec![
        support::encode_node(11, Orientation::Forward),
        support::encode_node(13, Orientation::Forward),
        support::encode_node(14, Orientation::Forward),
        support::encode_node(16, Orientation::Forward),
        support::encode_node(17, Orientation::Forward)
    ]);
    let path_len = 5;
    let path_interval = 0..5;
    let matches = 5;
    let edits = 0;
    let mapq = None;
    let score = None;
    let base_quality = None;
    let difference = None;
    let pair = None;
    let optional = Vec::new();
    let mut alignment = Alignment {
        name, seq_len, seq_interval,
        path, path_len, path_interval,
        matches, edits,
        mapq, score,
        base_quality, difference, pair,
        optional
    };

    let original = alignment.clone();
    let result = alignment.set_relative_information(
        &index, 0, &paths_by_sample, Some(&mut used_paths)
    );
    assert!(result.is_ok(), "Failed to set relative information: {}", result.err().unwrap());

    // Only name and path fields should have changed.
    check_alignment(&alignment, &original, 0, true, false);

    // Check the relative fields.
    // This is the third path that visits (starts from) node 11 forward.
    let start_pos = Pos::new(support::encode_node(11, Orientation::Forward), 2);
    let start_pos = TargetPath::StartPosition(start_pos);
    assert_eq!(alignment.path, start_pos, "Wrong GBWT starting position for the target path");

    // Check that the path is marked as used.
    assert_eq!(used_paths.bit(3), true, "Path 3 was not marked as used");
    assert_eq!(used_paths.count_ones(), 1, "Wrong number of used paths");

    // Check that repeated calls do not change the alignment.
    let original = alignment.clone();
    let result = alignment.set_relative_information(
        &index, 0, &paths_by_sample, Some(&mut used_paths)
    );
    assert!(result.is_ok(), "Failed to set relative information again: {}", result.err().unwrap());
    assert_eq!(alignment, original, "Alignment was changed by a repeated call");
    assert_eq!(used_paths.bit(3), true, "Used path was reset by a repeated call");
    assert_eq!(used_paths.count_ones(), 1, "Used paths were changed by a repeated call");
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: encoding and decoding.

// This assumes that the alignment is relative to a GBWT index.
fn check_encode_decode(alignment: &Alignment, prefix: &str, alphabet: &[u8], line: usize) {
    let query = &alignment.name[prefix.len()..];
    let target = if let TargetPath::StartPosition(pos) = alignment.path { pos } else { unreachable!() };
    let numbers = alignment.encode_numbers();
    let quality = alignment.encode_base_quality(alphabet);
    let difference = alignment.encode_difference();
    let pair = alignment.encode_pair(prefix.len());

    let decoded = Alignment::decode(
        prefix, query, target.node, &numbers,
        quality.as_ref().map(Vec::as_slice),
        alphabet, difference.as_ref().map(Vec::as_slice),
        pair.as_ref().map(Vec::as_slice)
    );
    assert!(decoded.is_ok(), "Failed to decode alignment {}: {}", line, decoded.err().unwrap());
    let decoded = decoded.unwrap();

    // TODO: Check the optional fields.
    check_alignment(&decoded, alignment, line, false, true);
}

#[test]
fn alignment_encode_decode() {
    // Parse the alignments.
    let filename = formats::get_test_data("good.gaf");
    let mut alignments = parse_alignments(&filename, false);
    assert_eq!(alignments.len(), 8, "Unexpected number of parsed alignments");

    // Set the relative information manually, as we do not have a GBWT index.
    let pos = Pos::new(support::encode_node(10, Orientation::Forward), 0);
    let empty = Pos::new(gbwt::ENDMARKER, 0);
    for (i, alignment) in alignments.iter_mut().enumerate() {
        alignment.path = TargetPath::StartPosition(if i < 6 { pos } else { empty });
    }

    // Encode and decode the alignments.
    let prefix = "";
    let alphabet = b"?";
    for (i, aln) in alignments.iter().enumerate() {
        check_encode_decode(aln, prefix, alphabet, i + 1);
    }
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: integration.
// This is the haplotype sampling test case from vg.

fn integration_test(gaf_file: &'static str, gbwt_file: &'static str) {
    // Parse the alignments.
    let gaf_file = formats::get_test_data(gaf_file);
    let mut alignments = parse_alignments(&gaf_file, false);
    assert_eq!(alignments.len(), 12439, "Unexpected number of parsed alignments");

    // Load the GBWT index and prepare the structures.
    let gbwt_file = formats::get_test_data(gbwt_file);
    let index: GBWT = serialize::load_from(&gbwt_file).unwrap();
    let metadata = index.metadata().unwrap();
    let result = Alignment::paths_by_sample(&metadata);
    assert!(result.is_ok(), "Failed to get a list of paths by sample: {}", result.err().unwrap());
    let paths_by_sample = result.unwrap();
    let mut used_paths = RawVector::with_len(metadata.paths(), false);

    // Set relative information for all alignments.
    let prefix = "A00744:46:HV3C3DSXX:2:";
    for (i, alignment) in alignments.iter_mut().enumerate() {
        let result = alignment.set_relative_information(
            &index, prefix.len(), &paths_by_sample, Some(&mut used_paths)
        );
        assert!(result.is_ok(), "Failed to set relative information for alignment {}: {}", i + 1, result.err().unwrap());
    }

    // Check that we managed to match each alignment to a different path.
    assert_eq!(used_paths.count_ones(), used_paths.len(), "Not all paths were used");

    // Encode and decode the alignments.
    let alphabet = b"#,:F";
    for (i, aln) in alignments.iter().enumerate() {
        check_encode_decode(aln, prefix, alphabet, i + 1);
    }
}

#[test]
fn alignment_real() {
    integration_test("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt");
}

#[test]
fn alignment_real_gzipped() {
    integration_test("micb-kir3dl1_HG003.gaf.gz", "micb-kir3dl1_HG003.gbwt");
}

//-----------------------------------------------------------------------------

// Tests for `QualityEncoder`.

fn create_quality_encoder(alphabet: &[u8], dictionary: &[(usize, usize)], name: &str) -> QualityEncoder {
    let encoder = QualityEncoder::new(alphabet, dictionary);
    assert!(encoder.is_some(), "Failed to create a quality encoder for {}", name);
    encoder.unwrap()
}

fn check_quality_encoder(encoder: &QualityEncoder, sequence: &[u8], is_rle: bool, name: &str) {
    let encoded = encoder.encode(sequence);
    assert!(encoded.is_ok(), "Failed to encode the quality string for {}: {}", name, encoded.err().unwrap());
    let encoded = encoded.unwrap();
    assert_eq!(QualityEncoder::is_rle(&encoded), is_rle, "Wrong RLE status for {}", name);
    let decoded = encoder.decode(&encoded, sequence.len());
    assert!(decoded.is_ok(), "Failed to decode the quality string for {}: {}", name, decoded.err().unwrap());
    let decoded = decoded.unwrap();
    assert_eq!(decoded, sequence, "Wrong decoded quality string for {}", name);
}

fn check_invalid_sequence(encoder: &QualityEncoder, sequence: &[u8], name: &str) {
    let encoded = encoder.encode(sequence);
    assert!(encoded.is_err(), "Encoded an invalid quality string for {}", name);
}

#[test]
fn quality_empty() {
    let alphabet = b"";
    let dictionary = [];
    let name = "empty alphabet";
    let encoder = create_quality_encoder(alphabet, &dictionary, name);

    let sequence = b"";
    check_quality_encoder(&encoder, sequence, false, name);
    let invalid = b"ACGT";
    check_invalid_sequence(&encoder, invalid, name);
}

#[test]
fn quality_unary() {
    let alphabet = b"A";
    let dictionary = [(1, 1)];
    let name = "unary alphabet";
    let encoder = create_quality_encoder(alphabet, &dictionary, name);

    let empty = b"";
    check_quality_encoder(&encoder, empty, false, "unary alphabet, empty sequence");
    let sequence = b"AAAAA";
    check_quality_encoder(&encoder, sequence, false, name);
    let invalid = b"ACGT";
    check_invalid_sequence(&encoder, invalid, name);
}

#[test]
fn quality_huffman() {
    // A: 0, C: 10, G: 110, T: 111
    let alphabet = b"ACGT";
    let dictionary = [(1, 1), (2, 1), (3, 2)];
    let name = "general alphabet";
    let encoder = create_quality_encoder(alphabet, &dictionary, name);

    let empty = b"";
    check_quality_encoder(&encoder, empty, false, "general alphabet, empty sequence");
    let sequence = b"ACGTACGTACGT";
    check_quality_encoder(&encoder, sequence, false, name);
    let invalid = b"ACGTN";
    check_invalid_sequence(&encoder, invalid, name);
}

#[test]
fn quality_large_alphabet() {
    // A: 0, B: 100, C: 101, D: 1100, E: 1101, F: 111000, ...
    let alphabet = b"ABCDEFGHIJKLM";
    let dictionary = [(1, 1), (3, 2), (4, 2), (6, 8)];
    let name = "large alphabet";
    let encoder = create_quality_encoder(alphabet, &dictionary, name);

    let empty = b"";
    check_quality_encoder(&encoder, empty, false, "large alphabet, empty sequence");
    let sequence = b"ABCDEFGHIJKLMABCDEFGHIJKLM";
    check_quality_encoder(&encoder, sequence, false, name);
    let invalid = b"ACGT";
    check_invalid_sequence(&encoder, invalid, name);
}

// huffman + rle
#[test]
fn quality_with_rle() {
    let alphabet = b"ACGT";
    let dictionary = [(1, 1), (2, 1), (3, 2)];
    let name = "possible rle";
    let encoder = create_quality_encoder(alphabet, &dictionary, name);

    let empty = b"";
    check_quality_encoder(&encoder, empty, false, "possible rle, empty sequence");
    let rle = b"AAAAAAAAACAAAAAAAAG";
    check_quality_encoder(&encoder, rle, true, name);
    let wrong_symbol = b"CCCCCCCCCACCCCCCCG";
    check_quality_encoder(&encoder, wrong_symbol, false, "possible rle, wrong symbol in the runs");
    let invalid = b"ACGTN";
    check_invalid_sequence(&encoder, invalid, name);

    let balanced = [(2, 4)];
    let balanced_encoder = create_quality_encoder(alphabet, &balanced, "balanced alphabet");
    check_quality_encoder(&balanced_encoder, rle, false, "balanced alphabet, long runs");
}

//-----------------------------------------------------------------------------

// Tests for `Difference`.

fn check_difference(seq: &[u8], truth: &[Difference], name: &str, normalize: bool) {
    let result = Difference::parse(seq);
    assert!(result.is_ok(), "Failed to parse the difference string for {}: {}", name, result.err().unwrap());
    let mut result = result.unwrap();
    if normalize {
        result = Difference::normalize(result);
    }
    assert_eq!(result.len(), truth.len(), "Wrong number of operations for {}", name);
    for i in 0..truth.len() {
        assert_eq!(result[i], truth[i], "Wrong operation at position {} for {}", i, name);
    }
}

#[test]
fn difference_empty() {
    let seq = b"";
    let truth = [];
    check_difference(seq, &truth, "empty", false);
}

#[test]
fn difference_single() {
    {
        // Matching sequence.
        let seq = b"=ACGT";
        let truth = [ Difference::Match(4) ];
        check_difference(seq, &truth, "matching sequence", false);
    }
    {
        // Match length.
        let seq = b":5";
        let truth = [ Difference::Match(5) ];
        check_difference(seq, &truth, "match length", false);
    }
    {
        // Mismatch.
        let seq = b"*AC";
        let truth = [ Difference::Mismatch(b'C') ];
        check_difference(seq, &truth, "mismatch", false);
    }
    {
        // Insertion.
        let seq = b"+ACGT";
        let truth = [ Difference::Insertion(b"ACGT".to_vec()) ];
        check_difference(seq, &truth, "insertion", false);
    }
    {
        // Deletion.
        let seq = b"-ACGT";
        let truth = [ Difference::Deletion(4) ];
        check_difference(seq, &truth, "deletion", false);
    }
}

#[test]
fn difference_multi() {
    // All of these are from real reads mapped with Giraffe.
    {
        let seq = b":24*CG:78-C:35*TG:2*TG:1*TG:6";
        let truth = [
            Difference::Match(24),
            Difference::Mismatch(b'G'),
            Difference::Match(78),
            Difference::Deletion(1),
            Difference::Match(35),
            Difference::Mismatch(b'G'),
            Difference::Match(2),
            Difference::Mismatch(b'G'),
            Difference::Match(1),
            Difference::Mismatch(b'G'),
            Difference::Match(6)
        ];
        check_difference(seq, &truth, "ERR3239454.129477191 fragment 2", false);
    }
    {
        let seq = b":129+TTGAGGGGGTATAAGAGGTCG";
        let truth = [
            Difference::Match(129),
            Difference::Insertion(b"TTGAGGGGGTATAAGAGGTCG".to_vec())
        ];
        check_difference(seq, &truth, "ERR3239454.11251898 fragment 1", false);
    }
    {
        let seq = b":111*GT:20*CA:6*GT:2*GT:3*GT:3";
        let truth = [
            Difference::Match(111),
            Difference::Mismatch(b'T'),
            Difference::Match(20),
            Difference::Mismatch(b'A'),
            Difference::Match(6),
            Difference::Mismatch(b'T'),
            Difference::Match(2),
            Difference::Mismatch(b'T'),
            Difference::Match(3),
            Difference::Mismatch(b'T'),
            Difference::Match(3)
        ];
        check_difference(seq, &truth, "ERR3239454.97848632 fragment 1", false);
    }
}

#[test]
fn difference_case() {
    {
        // Lower case mismatch.
        let seq = b"*ac";
        let truth = [ Difference::Mismatch(b'C') ];
        check_difference(seq, &truth, "lower case mismatch", false);
    }
    {
        // Lower case insertion.
        let seq = b"+acgt";
        let truth = [ Difference::Insertion(b"ACGT".to_vec()) ];
        check_difference(seq, &truth, "lower case insertion", false);
    }
}

fn invalid_difference(seq: &[u8], name: &str) {
    let result = Difference::parse(seq);
    assert!(result.is_err(), "Parsed an invalid difference string: {}", name);
}

#[test]
fn difference_invalid() {
    invalid_difference(b"ABC", "invalid operation at start");
    invalid_difference(b"+AC:", "empty operation at end");
    invalid_difference(b"~AC30AC", "intron / splice operation");
    invalid_difference(b":123A", "match length is not an integer");
    invalid_difference(b"+GATTACA:", "empty match length");
    invalid_difference(b"*ACT", "mismatch is too long");
    invalid_difference(b"*A", "mismatch is too short");
}

#[test]
fn difference_len() {
    {
        let diff = Difference::Match(42);
        assert_eq!(diff.len(), 42, "Wrong length for a match");
    }
    {
        let diff = Difference::Mismatch(b'A');
        assert_eq!(diff.len(), 1, "Wrong length for a mismatch");
    }
    {
        let diff = Difference::Insertion(b"GATTACA".to_vec());
        assert_eq!(diff.len(), 7, "Wrong length for an insertion");
    }
    {
        let diff = Difference::Deletion(42);
        assert_eq!(diff.len(), 42, "Wrong length for a deletion");
    }
}

#[test]
fn difference_normalize() {
    {
        let seq = b"=ACGT:4*AC*GT:2:3=ACT*AC-ACGT-ACGT+GATTACA+CAT";
        let truth = [
            Difference::Match(8),
            Difference::Mismatch(b'C'),
            Difference::Mismatch(b'T'),
            Difference::Match(8),
            Difference::Mismatch(b'C'),
            Difference::Deletion(8),
            Difference::Insertion(b"GATTACACAT".to_vec()),
        ];
        check_difference(seq, &truth, "normalized", true);
    }
    {
        let seq = b"=:0+-*AC:30=GATTA+-:0=";
        let truth = [
            Difference::Mismatch(b'C'),
            Difference::Match(35),
        ];
        check_difference(seq, &truth, "normalized with empty operations", true);
    }
}

//-----------------------------------------------------------------------------

// Tests for `TypedField`.

#[test]
fn typed_field_empty() {
    let field = TypedField::parse(b"");
    assert!(field.is_err(), "Empty string was parsed as a typed field");
}

// This assumes that the sequence is a valid typed field.
fn check_typed_field(seq: &[u8], name: &str) {
    let field = TypedField::parse(seq);
    assert!(field.is_ok(), "Failed to parse the typed field for {}", name);
    let field = field.unwrap();
    assert_eq!(field.tag(), seq[0..2], "Wrong tag for {}", name);

    match field {
        TypedField::Char(_, value) => {
            assert_eq!(seq[3], b'A', "Parsed the type as a char for {}", name);
            assert_eq!(value, seq[5], "Wrong value for {}", name);
        }
        TypedField::String(_, value) => {
            assert_eq!(seq[3], b'Z', "Parsed the type as a string for {}", name);
            assert_eq!(value, &seq[5..], "Wrong value for {}", name);
        },
        TypedField::Int(_, value) => {
            assert_eq!(seq[3], b'i', "Parsed the type as an integer for {}", name);
            let value_seq = str::from_utf8(&seq[5..]).unwrap();
            let truth = value_seq.parse::<isize>().unwrap();
            assert_eq!(value, truth, "Wrong value for {}", name);
        },
        TypedField::Float(_, value) => {
            assert_eq!(seq[3], b'f', "Parsed the type as a float for {}", name);
            let value_seq = str::from_utf8(&seq[5..]).unwrap();
            let truth = value_seq.parse::<f64>().unwrap();
            assert_eq!(value, truth, "Wrong value for {}", name);
        },
        TypedField::Bool(_, value) => {
            assert_eq!(seq[3], b'b', "Parsed the type as a boolean for {}", name);
            let truth = match seq[5] {
                b'0' => false,
                b'1' => true,
                _ => panic!("Invalid boolean value for {}", name),
            };
            assert_eq!(value, truth, "Wrong value for {}", name);
        },
    }
}

#[test]
fn typed_field_char() {
    check_typed_field(b"tp:A:P", "char");
}

#[test]
fn typed_field_string() {
    check_typed_field(b"cg:Z:6M", "string");
    check_typed_field(b"ab:Z:", "empty string");
}

#[test]
fn typed_field_int() {
    check_typed_field(b"AS:i:150", "positive integer");
    check_typed_field(b"cd:i:-3", "negative integer");
    check_typed_field(b"ef:i:0", "zero");
}

#[test]
fn typed_field_float() {
    check_typed_field(b"gh:f:22", "positive whole number");
    check_typed_field(b"ij:f:-10", "negative whole number");
    check_typed_field(b"dv:f:0", "zero whole number");
    check_typed_field(b"kl:f:3.0", "positive whole number with decimal point");
    check_typed_field(b"mn:f:-2.0", "negative whole number with decimal point");

    check_typed_field(b"kl:f:3.14", "positive decimal number");
    check_typed_field(b"mn:f:-2.718", "negative decimal number");
    check_typed_field(b"op:f:0.0", "zero decimal number");
    check_typed_field(b"qr:f:.5", "positive decimal number without zero");
    check_typed_field(b"st:f:-.25", "negative decimal number without zero");
}

#[test]
fn typed_field_bool() {
    check_typed_field(b"pd:b:0", "false");
    check_typed_field(b"qr:b:1", "true");
}

fn invalid_typed_field(seq: &[u8], name: &str) {
    let field = TypedField::parse(seq);
    assert!(field.is_err(), "Parsed an invalid typed field for {}", name);
}

#[test]
fn typed_field_invalid() {
    // Invalid format.
    invalid_typed_field(b"AS::150", "missing type");
    invalid_typed_field(b"AS:i150", "missing separator");
    invalid_typed_field(b"ASi150", "missing separators");
    invalid_typed_field(b":i:150", "missing tag");
    invalid_typed_field(b"A:i:150", "too short tag");
    invalid_typed_field(b"ASC:i:150", "too long tag");

    // Invalid type.
    invalid_typed_field(b"AS:J:150", "JSON type (not supported)");
    invalid_typed_field(b"AS:H:150", "hexadecimal type (not supported)");
    invalid_typed_field(b"AS:B:150", "array type (not supported)");
    invalid_typed_field(b"AS:?:150", "unknown type");

    // Invalid value for char.
    invalid_typed_field(b"tp:A:", "empty char value");
    invalid_typed_field(b"tp:A:AC", "too long char value");
    invalid_typed_field(b"tp:A: A", "char value with leading space");
    invalid_typed_field(b"tp:A:A ", "char value with trailing space");

    // Invalid value for integer.
    invalid_typed_field(b"AS:i:", "empty integer value");
    invalid_typed_field(b"AS:i: 150", "integer value with leading space");
    invalid_typed_field(b"AS:i:150 ", "integer value with trailing space");
    invalid_typed_field(b"AS:i:3.14", "decimal number as integer value");
    invalid_typed_field(b"AS:i:3A", "invalid character in integer value");

    // Invalid value for float.
    invalid_typed_field(b"AS:f:", "empty float value");
    invalid_typed_field(b"AS:f: 3.14", "float value with leading space");
    invalid_typed_field(b"AS:f:3.14 ", "float value with trailing space");
    invalid_typed_field(b"AS:f:3.14A", "invalid character in float value");

    // Invalid value for boolean.
    invalid_typed_field(b"pd:b:", "empty boolean value");
    invalid_typed_field(b"pd:b:01", "too long boolean value");
    invalid_typed_field(b"pd:b:2", "wrong number in boolean value");
    invalid_typed_field(b"pd:b:?", "invalid character in boolean value");
    invalid_typed_field(b"pd:b: 0", "boolean value with leading space");
    invalid_typed_field(b"pd:b:0 ", "boolean value with trailing space");
}

//-----------------------------------------------------------------------------
