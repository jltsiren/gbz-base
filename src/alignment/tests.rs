use super::*;

use crate::{Subgraph, SubgraphQuery, GraphReference};
use crate::utils;

use gbz::GBZ;
use simple_sds::serialize;

use rand::Rng;

use std::collections::{HashMap, BTreeSet};
use std::io::BufRead;
use std::path::PathBuf;
use std::vec;

//-----------------------------------------------------------------------------

// Tests for `Alignment`: parsing.

// Creates an alignment object for an unaligned sequence with no optional fields.
fn empty_alignment(name: &str, seq_len: usize) -> Alignment {
    let mut result = Alignment::default();
    result.name = String::from(name);
    result.seq_len = seq_len;
    result
}

// Returns the GAF alignment lines from the given file.
fn read_gaf_lines(filename: &PathBuf) -> Vec<Vec<u8>> {
    let reader = utils::open_file(filename);
    assert!(reader.is_ok(), "Failed to open the test file: {}", reader.err().unwrap());
    let mut reader = reader.unwrap();

    let mut result = Vec::new();
    let mut line_num = 1;
    loop {
        let mut buf: Vec<u8> = Vec::new();
        let len = reader.read_until(b'\n', &mut buf);
        assert!(len.is_ok(), "Failed to read line {}: {}", line_num, len.err().unwrap());
        if formats::is_gaf_header_line(&buf) {
            line_num += 1;
            continue;
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.is_empty() {
            break;
        }
        result.push(buf);
        line_num += 1;
    };

    result
}

// Parses the alignments from the test file.
// Returns them as a vector, unless told to expect a parse error.
// Also validates the parsed alignments.
fn parse_alignments(filename: &PathBuf, expect_error: bool) -> Vec<Alignment> {
    let lines = read_gaf_lines(filename);

    let mut result = Vec::new();
    for (line_num, line) in lines.iter().enumerate() {
        let alignment = Alignment::from_gaf(line);
        if expect_error {
            assert!(alignment.is_err(), "Expected a parse error on line {}", line_num + 1);
        } else {
            assert!(alignment.is_ok(), "Failed to parse line {}: {}", line_num + 1, alignment.err().unwrap());
            let alignment = alignment.unwrap();
            let valid = alignment.validate();
            assert!(valid.is_ok(), "Parsed alignment on line {} is inconsistent: {}", line_num + 1, valid.err().unwrap());
            result.push(alignment);
        }
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

fn known_good_alignments() -> Vec<Alignment> {
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
    let base_quality = vec![b'?'; 100];
    let difference = vec![
        Difference::Match(20),
        Difference::Mismatch(b'C'), Difference::Mismatch(b'T'),
        Difference::Match(20),
        Difference::Deletion(8),
        Difference::Match(20),
        Difference::Insertion(b"CAT".to_vec()),
        Difference::Match(20)
    ];
    let pair = None;
    let optional = vec![
        TypedField::String([b'x', b'x'], b"unknown".to_vec()),
    ];

    // The file contains:
    // * three variants of the same alignment
    // * three additional variants that are paired
    // * three variants of unaligned sequences
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

    let mut unaligned = empty_alignment("unaligned", 10);
    unaligned.seq_interval = 0..10;
    unaligned.difference.push(Difference::Insertion(b"GATTACAGAT".to_vec()));
    unaligned.edits = 10;

    vec![forward, reverse, no_mapq, first, second, same, empty, missing_values, unaligned]
}

#[test]
fn alignment_known_good() {
    let truth = known_good_alignments();
    let filename = utils::get_test_data("good.gaf");
    let alignments = parse_alignments(&filename, false);
    assert_eq!(alignments.len(), truth.len(), "Wrong number of alignments in the test file");
    for i in 0..truth.len() {
        check_alignment(&alignments[i], &truth[i], i + 1, false, false);
    }

    // These alignments will survive the round-trip from/to GAF intact.
    // 1 is on the reverse strand, while 7 has non-canonical empty fields.
    // In 8, missing numerical fields get values from the difference string.
    let round_trip = vec![
        0, 2, 3, 4, 5, 6,
    ];
    let gaf_lines = read_gaf_lines(&filename);
    let target_sequence = b"AAAAAAAAAACCCCCCCCCCCCCCCCCCCCAGAAAAAAAAAAAAAAAAAAAAGATTACATGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAACCCCCCCCCC".to_vec();
    for index in round_trip {
        let line = alignments[index].to_gaf(&target_sequence);
        assert_eq!(line, gaf_lines[index], "Line {} did not survive the round trip", index + 1);
    }
}

#[test]
fn alignment_known_good_no_header() {
    // This is the same test as above, but the GAF file has no header lines.
    let truth = known_good_alignments();
    let filename = utils::get_test_data("no_header.gaf");
    let alignments = parse_alignments(&filename, false);
    assert_eq!(alignments.len(), truth.len(), "Wrong number of alignments in the test file");
    for i in 0..truth.len() {
        check_alignment(&alignments[i], &truth[i], i + 1, false, false);
    }
}

#[test]
fn alignment_known_bad() {
    let filename = utils::get_test_data("bad.gaf");
    let alignments = parse_alignments(&filename, true);
    assert_eq!(alignments.len(), 0, "There should be no valid alignments in the test file");
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: operations.

#[test]
fn alignment_default() {
    let default = Alignment::default();
    assert!(default.is_unaligned(), "A default alignment should be unaligned");
    assert!(!default.is_full(), "A default alignment should not be a full alignment");
    assert!(!default.is_exact(), "An empty alignment should not be exact");
    assert!(!default.has_difference_string(), "A default alignment should not have a difference string");

    assert!(default.name.is_empty(), "A default alignment should not have a name");
    assert_eq!(default.seq_len, 0, "A default alignment should have a sequence length of 0");
    assert!(default.seq_interval.is_empty(), "A default alignment should have an empty sequence interval");
    assert!(default.has_target_path(), "A default alignment should have a target path");
    assert!(!default.has_non_empty_target_path(), "A default alignment should not have a non-empty target path");
    assert!(default.target_path().unwrap().is_empty(), "A default alignment should have an empty target path");
    assert_eq!(default.path_len, 0, "A default alignment should have a target path length of 0");
    assert!(default.path_interval.is_empty(), "A default alignment should have an empty target path interval");
    assert_eq!(default.matches, 0, "A default alignment should have 0 matches");
    assert_eq!(default.edits, 0, "A default alignment should have 0 edits");
    assert!(default.mapq.is_none(), "A default alignment should not have a mapping quality");
    assert!(default.score.is_none(), "A default alignment should not have a score");
    assert!(default.base_quality.is_empty(), "A default alignment should not have base quality values");
    assert!(default.difference.is_empty(), "A default alignment should have an empty difference string");
    assert!(default.pair.is_none(), "A default alignment should not have a paired read");
    assert!(default.optional.is_empty(), "A default alignment should not have any optional fields");
}

#[test]
fn alignment_operations() {
    let mut aln = empty_alignment("test", 10);
    aln.seq_interval = 0..10;
    aln.path = TargetPath::StartPosition(Pos::new(support::encode_node(1, Orientation::Forward), 0));
    aln.path_len = 12;
    aln.path_interval = 1..11;
    aln.matches = 10;

    assert!(!aln.is_unaligned(), "The alignment should be aligned");
    assert!(aln.is_full(), "The alignment should be full");
    assert!(aln.is_exact(), "The alignment should be exact");
    assert!(!aln.has_difference_string(), "The alignment should not have a difference string");

    // Various ways of making the alignment imperfect.
    {
        let mut wrong_start = aln.clone();
        wrong_start.seq_interval.start += 1;
        wrong_start.path_interval.start += 1;
        wrong_start.matches -= 1;
        assert!(!wrong_start.is_full(), "An alignment starting after query start should not be full");
    }
    {
        let mut wrong_end = aln.clone();
        wrong_end.seq_interval.end -= 1;
        wrong_end.path_interval.end -= 1;
        wrong_end.matches -= 1;
        assert!(!wrong_end.is_full(), "An alignment ending before query end should not be full");
    }
    {
        let mut has_edits = aln.clone();
        has_edits.matches -= 1;
        has_edits.edits += 1;
        assert!(!has_edits.is_exact(), "An alignment with edits should not be exact");
    }

    // We have a GBWT start position.
    assert!(!aln.has_target_path(), "GBWT start position should not be a target path");
    assert!(!aln.has_non_empty_target_path(), "GBWT start position should not be a non-empty target path");
    assert!(aln.target_path().is_none(), "GBWT start position should not be returned as a target path");
    assert!(aln.min_handle().is_none(), "GBWT start position should not have a minimum handle");
    assert!(aln.max_handle().is_none(), "GBWT start position should not have a maximum handle");

    // Set the target path manually.
    let correct_path = vec![
        support::encode_node(1, Orientation::Forward),
        support::encode_node(2, Orientation::Forward),
    ];
    let min_handle = correct_path.iter().min().cloned();
    let max_handle = correct_path.iter().max().cloned();
    aln.set_target_path(correct_path.clone());
    assert!(aln.has_target_path(), "Target path should exist after setting it");
    assert!(aln.has_non_empty_target_path(), "Target path should be non-empty after setting it");
    assert_eq!(aln.target_path(), Some(correct_path.as_slice()), "Target path should match the set path");
    assert_eq!(aln.min_handle(), min_handle, "Minimum handle should match the set path");
    assert_eq!(aln.max_handle(), max_handle, "Maximum handle should match the set path");

    // Setting it again should update the target path.
    let new_path = vec![
        support::encode_node(1, Orientation::Forward),
        support::encode_node(3, Orientation::Forward),
    ];
    let min_handle = new_path.iter().min().cloned();
    let max_handle = new_path.iter().max().cloned();
    aln.set_target_path(new_path.clone());
    assert!(aln.has_target_path(), "Target path should still exist after setting it again");
    assert!(aln.has_non_empty_target_path(), "Target path should still be non-empty after setting it again");
    assert_eq!(aln.target_path(), Some(new_path.as_slice()), "Target path should match the new path");
    assert_eq!(aln.min_handle(), min_handle, "Minimum handle should match the new path");
    assert_eq!(aln.max_handle(), max_handle, "Maximum handle should match the new path");
}

//-----------------------------------------------------------------------------

// Returns an appropriate end for the block, ensuring that it does not mix unaligned and aligned reads.
fn block_end(alignments: &[Alignment], start: usize, block_size: usize) -> usize {
    if start >= alignments.len() {
        return start;
    }

    let mut end = start + 1;
    let unaligned = alignments[start].is_unaligned();
    while end < alignments.len() && end < start + block_size && unaligned == alignments[end].is_unaligned() {
        end += 1;
    }
    end
}

fn check_encode_decode(block: &[Alignment], index: &GBWT, first_id: usize) {
    let range = first_id..(first_id + block.len());
    let encoded = AlignmentBlock::new(block, index, first_id);
    if let Err(message) = encoded {
        panic!("Failed to encode the alignment block {}..{}: {}", range.start, range.end, message);
    }
    let encoded = encoded.unwrap();
    assert_eq!(encoded.len(), block.len(), "Wrong number of alignments in the encoded block {}..{}", range.start, range.end);

    let decoded = encoded.decode();
    assert!(decoded.is_ok(), "Failed to decode the alignment block {}..{}: {}", range.start, range.end, decoded.err().unwrap());
    let mut decoded = decoded.unwrap();
    assert_eq!(decoded.len(), block.len(), "Wrong number of alignments in the decoded block {}..{}", range.start, range.end);

    for (i, aln) in decoded.iter_mut().enumerate() {
        aln.extract_target_path(index);
        // NOTE: We do not store the true target path length in the alignment block.
        // Because we do not have the reference graph here, we simply assume that it is correct.
        assert_eq!(aln.path_len, aln.path_interval.end, "Non-zero right flank for target path of alignment {} in decoded block {}..{}", i, range.start, range.end);
        aln.path_len = block[i].path_len;
        assert_eq!(*aln, block[i], "Wrong alignment {} in the decoded block {}..{}", range.start + i, range.start, range.end);
    }
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: encode/decode integration.
// This is the haplotype sampling test case from vg.

fn integration_test(gaf_file: &'static str, gbwt_file: &'static str, block_size: usize) {
    // Parse the alignments.
    let gaf_file = utils::get_test_data(gaf_file);
    let mut alignments = parse_alignments(&gaf_file, false);
    assert_eq!(alignments.len(), 12439, "Unexpected number of parsed alignments");

    // Remove the optional fields as we do not encode them for now.
    for aln in &mut alignments {
        aln.optional.clear();
    }

    // Load the GBWT index.
    let gbwt_file = utils::get_test_data(gbwt_file);
    let index: GBWT = serialize::load_from(&gbwt_file).unwrap();

    // Encode and decode the alignments.
    let mut start = 0;
    while start < alignments.len() {
        let end = block_end(&alignments, start, block_size);
        let block = &alignments[start..end];
        check_encode_decode(block, &index, start);
        start = end;
    }
}

#[test]
fn alignment_real() {
    integration_test("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt", 1234);
}

#[test]
fn alignment_bidirectional() {
    integration_test("micb-kir3dl1_HG003.gaf", "bidirectional.gbwt", 1001);
}

#[test]
fn alignment_real_gzipped() {
    integration_test("micb-kir3dl1_HG003.gaf.gz", "micb-kir3dl1_HG003.gbwt", 789);
}

//-----------------------------------------------------------------------------

fn check_alignment_iter(aln: &Alignment, truth: &[Mapping], expect_fail: bool) {
    // This is a convenient hack. We don't need a graph when the handle is also the node length.
    let sequence_len = Arc::new(|handle| Some(handle));

    let iter = aln.iter(sequence_len);
    if expect_fail {
        assert!(iter.is_none(), "Expected no iterator for {}", aln.name);
        return;
    }

    assert!(iter.is_some(), "Failed to build the iterator for {}", aln.name);
    let mut iter = iter.unwrap();
    for (i, expected) in truth.iter().enumerate() {
        let actual = iter.next().unwrap();
        assert_eq!(actual, *expected, "Alignment mismatch at index {} of {}: {:?} != {:?}", i, aln.name, actual, expected);
    }
}

// Creates a perfect alignment that matches a single node entirely.
fn perfect_alignment(name: &str, seq_len: usize, with_path: bool, with_difference: bool) -> Alignment {
    let mut aln = empty_alignment(name, seq_len);
    aln.seq_interval = 0..seq_len;
    if with_path {
        aln.path = TargetPath::Path(vec![seq_len]);
    } else {
        aln.path = TargetPath::StartPosition(Pos::new(seq_len, 0));
    }
    aln.path_len = seq_len;
    aln.path_interval = 0..seq_len;
    aln.matches = seq_len;
    if with_difference {
        aln.difference.push(Difference::Match(seq_len));
    }
    aln
}

// TODO: Make this a public constructor for Alignment? But then we need node lengths.
fn create_alignment(
    name: &str,
    seq_len: usize, seq_start: usize,
    path: Vec<usize>, path_start: usize,
    difference: Vec<Difference>,
    sequence_len: Option<Arc<dyn Fn(usize) -> Option<usize> + '_>>
) -> Alignment {
    let mut seq_end = seq_start;
    let mut path_end = path_start;
    let mut matches = 0;
    let mut edits = 0;
    for op in difference.iter() {
        seq_end += op.query_len();
        path_end += op.target_len();
        if let Difference::Match(len) = op {
            matches += len;
        } else {
            edits += cmp::max(op.query_len(), op.target_len());
        }
    }

    let mut aln = empty_alignment(name, seq_len);
    aln.seq_interval = seq_start..seq_end;

    let sequence_len = if let Some(func) = sequence_len {
        func
    } else {
        Arc::new(|handle| Some(handle))
    };
    aln.path_len = path.iter().map(|&handle| sequence_len(handle).unwrap()).sum();

    aln.path = TargetPath::Path(path);
    aln.path_interval = path_start..path_end;
    aln.matches = matches;
    aln.edits = edits;
    aln.difference = difference;
    aln
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: Iterator.

#[test]
fn alignment_iter_special_cases() {
    let aln = empty_alignment("empty", 0);
    let truth: Vec<Mapping> = Vec::new();
    check_alignment_iter(&aln, &truth, false);

    let aln = empty_alignment("unaligned", 150);
    check_alignment_iter(&aln, &truth, false);

    let aln = perfect_alignment("no-target-path", 150, false, true);
    check_alignment_iter(&aln, &truth, true);

    let aln = perfect_alignment("no-difference-string", 150, true, false);
    check_alignment_iter(&aln, &truth, true);

    let aln = create_alignment(
        "insertion-only",
        7, 0,
        vec![], 0,
        vec![Difference::Insertion(b"GATTACA".to_vec())],
        None
    );
    let truth = vec![Mapping::unaligned(0, Difference::Insertion(b"GATTACA".to_vec()))];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "insertion-only-with-path",
        7, 0,
        vec![30], 0,
        vec![Difference::Insertion(b"GATTACA".to_vec())],
        None
    );
    let truth = vec![Mapping::new(0, 30, 30, 0, Difference::Insertion(b"GATTACA".to_vec()))];
    check_alignment_iter(&aln, &truth, false);
}

#[test]
fn alignment_iter_perfect() {
    let aln = perfect_alignment("single-node", 100, true, true);
    let truth = vec![Mapping::new(0, 100, 100, 0, Difference::Match(100))];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "multi-node",
        100, 0,
        vec![40, 60], 0,
        vec![Difference::Match(100)],
        None
    );
    let truth = vec![
        Mapping::new(0, 40, 40, 0, Difference::Match(40)),
        Mapping::new(40, 60, 60, 0, Difference::Match(60)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "contained-single-node",
        100, 0,
        vec![50, 60], 5,
        vec![Difference::Match(100)],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 5, Difference::Match(45)),
        Mapping::new(45, 60, 60, 0, Difference::Match(55)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "contained-multi-node",
        100, 0,
        vec![30, 40, 50], 10,
        vec![Difference::Match(100)],
        None
    );
    let truth = vec![
        Mapping::new(0, 30, 30, 10, Difference::Match(20)),
        Mapping::new(20, 40, 40, 0, Difference::Match(40)),
        Mapping::new(60, 50, 50, 0, Difference::Match(40)),
    ];
    check_alignment_iter(&aln, &truth, false);
}

#[test]
fn alignment_iter_full() {
    let aln = create_alignment(
        "mappings-within-nodes",
        100, 0,
        vec![30, 40, 50], 0,
        vec![
            Difference::Match(20), Difference::Deletion(10),
            Difference::Match(35), Difference::Insertion(b"AC".to_vec()), Difference::Match(5),
            Difference::Mismatch(b'T'), Difference::Match(37)
        ],
        None
    );
    let truth = vec![
        Mapping::new(0, 30, 30, 0, Difference::Match(20)),
        Mapping::new(20, 30, 30, 20, Difference::Deletion(10)),
        Mapping::new(20, 40, 40, 0, Difference::Match(35)),
        Mapping::new(55, 40, 40, 35, Difference::Insertion(b"AC".to_vec())),
        Mapping::new(57, 40, 40, 35, Difference::Match(5)),
        Mapping::new(62, 50, 50, 0, Difference::Mismatch(b'T')),
        Mapping::new(63, 50, 50, 1, Difference::Match(37)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "multi-node-mappings",
        100, 0,
        vec![10, 20, 30, 40, 50], 9,
        vec![
            Difference::Match(2),
            Difference::Mismatch(b'A'), Difference::Match(50),
            Difference::Deletion(10), Difference::Match(47)
        ],
        None
    );
    let truth = vec![
        Mapping::new(0, 10, 10, 9, Difference::Match(1)),
        Mapping::new(1, 20, 20, 0, Difference::Match(1)),
        Mapping::new(2, 20, 20, 1, Difference::Mismatch(b'A')),
        Mapping::new(3, 20, 20, 2, Difference::Match(18)),
        Mapping::new(21, 30, 30, 0, Difference::Match(30)),
        Mapping::new(51, 40, 40, 0, Difference::Match(2)),
        Mapping::new(53, 40, 40, 2, Difference::Deletion(10)),
        Mapping::new(53, 40, 40, 12, Difference::Match(28)),
        Mapping::new(81, 50, 50, 0, Difference::Match(19)),
    ];
    check_alignment_iter(&aln, &truth, false);
}

#[test]
fn alignment_iter_partial() {
    let aln = create_alignment(
        "second-to-first",
        100, 10,
        vec![80, 30], 1,
        vec![Difference::Match(40), Difference::Insertion(b"GATTACA".to_vec()), Difference::Match(40)],
        None
    );
    let truth = vec![
        Mapping::new(10, 80, 80, 1, Difference::Match(40)),
        Mapping::new(50, 80, 80, 41, Difference::Insertion(b"GATTACA".to_vec())),
        Mapping::new(57, 80, 80, 41, Difference::Match(39)),
        Mapping::new(96, 30, 30, 0, Difference::Match(1)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "middle-to-middle",
        100, 10,
        vec![30, 40, 50], 15,
        vec![Difference::Match(35), Difference::Mismatch(b'A'), Difference::Match(44)],
        None
    );
    let truth = vec![
        Mapping::new(10, 30, 30, 15, Difference::Match(15)),
        Mapping::new(25, 40, 40, 0, Difference::Match(20)),
        Mapping::new(45, 40, 40, 20, Difference::Mismatch(b'A')),
        Mapping::new(46, 40, 40, 21, Difference::Match(19)),
        Mapping::new(65, 50, 50, 0, Difference::Match(25)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "last-to-second-to-last",
        100, 10,
        vec![40, 30, 50], 39,
        vec![
            Difference::Match(45),
            Difference::Insertion(b"GATTACA".to_vec()), Difference::Deletion(6),
            Difference::Match(28)
        ],
        None
    );
    let truth = vec![
        Mapping::new(10, 40, 40, 39, Difference::Match(1)),
        Mapping::new(11, 30, 30, 0, Difference::Match(30)),
        Mapping::new(41, 50, 50, 0, Difference::Match(14)),
        Mapping::new(55, 50, 50, 14, Difference::Insertion(b"GATTACA".to_vec())),
        Mapping::new(62, 50, 50, 14, Difference::Deletion(6)),
        Mapping::new(62, 50, 50, 20, Difference::Match(28)),
    ];
    check_alignment_iter(&aln, &truth, false);

    // This is the middle-to-middle test case, but with an unused node at the start and the end.
    let aln = create_alignment(
        "unused-nodes",
        100, 10,
        vec![20, 30, 40, 50, 15], 35,
        vec![Difference::Match(35), Difference::Mismatch(b'A'), Difference::Match(44)],
        None
    );
    let truth = vec![
        Mapping::new(10, 30, 30, 15, Difference::Match(15)),
        Mapping::new(25, 40, 40, 0, Difference::Match(20)),
        Mapping::new(45, 40, 40, 20, Difference::Mismatch(b'A')),
        Mapping::new(46, 40, 40, 21, Difference::Match(19)),
        Mapping::new(65, 50, 50, 0, Difference::Match(25)),
    ];
    check_alignment_iter(&aln, &truth, false);
}

// Insertions at node boundaries are supposed to be assigned to the end of the previous node.
// Except at the start of the alignment or if the previous node is not used in the alignment.
#[test]
fn alignment_iter_insertions() {
    let aln = create_alignment(
        "at-start",
        100, 0,
        vec![50, 48], 0,
        vec![Difference::Insertion(b"AC".to_vec()), Difference::Match(98)],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 0, Difference::Insertion(b"AC".to_vec())),
        Mapping::new(2, 50, 50, 0, Difference::Match(50)),
        Mapping::new(52, 48, 48, 0, Difference::Match(48)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "at-start-with-unused",
        100, 0,
        vec![10, 50, 48], 10,
        vec![Difference::Insertion(b"AC".to_vec()), Difference::Match(98)],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 0, Difference::Insertion(b"AC".to_vec())),
        Mapping::new(2, 50, 50, 0, Difference::Match(50)),
        Mapping::new(52, 48, 48, 0, Difference::Match(48)),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "at-end",
        100, 0,
        vec![50, 48], 0,
        vec![Difference::Match(98), Difference::Insertion(b"AC".to_vec())],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 0, Difference::Match(50)),
        Mapping::new(50, 48, 48, 0, Difference::Match(48)),
        Mapping::new(98, 48, 48, 48, Difference::Insertion(b"AC".to_vec())),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "at-end-with-unused",
        100, 0,
        vec![50, 48, 10], 0,
        vec![Difference::Match(98), Difference::Insertion(b"AC".to_vec())],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 0, Difference::Match(50)),
        Mapping::new(50, 48, 48, 0, Difference::Match(48)),
        Mapping::new(98, 48, 48, 48, Difference::Insertion(b"AC".to_vec())),
    ];
    check_alignment_iter(&aln, &truth, false);

    let aln = create_alignment(
        "in-the-midle",
        100, 0,
        vec![50, 48], 0,
        vec![Difference::Match(50), Difference::Insertion(b"AC".to_vec()), Difference::Match(48)],
        None
    );
    let truth = vec![
        Mapping::new(0, 50, 50, 0, Difference::Match(50)),
        Mapping::new(50, 50, 50, 50, Difference::Insertion(b"AC".to_vec())),
        Mapping::new(52, 48, 48, 0, Difference::Match(48)),
    ];
    check_alignment_iter(&aln, &truth, false);
}

#[test]
fn alignment_iter_real() {
    // Get the alignments.
    let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let alignments = parse_alignments(&gaf_file, false);

    // Load the graph.
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    let sequence_len = Arc::new(|handle| {
        let (node_id, _) = support::decode_node(handle);
        graph.sequence_len(node_id)
    });

    // We don't know what to expect, except that we want to iterate over the entire query/target intervals.
    for (i, aln) in alignments.iter().enumerate() {
        let iter = aln.iter(sequence_len.clone());
        assert!(iter.is_some(), "Failed to build the iterator for alignment {}", i + 1);
        let iter = iter.unwrap();

        let mut seq_offset = aln.seq_interval.start;
        let mut path_offset = aln.path_interval.start;
        for mapping in iter {
            seq_offset += mapping.seq_interval().len();
            path_offset += mapping.node_interval().len();
        }
        assert_eq!(seq_offset, aln.seq_interval.end, "Failed to reach the end of the query interval in alignment {}", i + 1);
        assert_eq!(path_offset, aln.path_interval.end, "Failed to reach the end of the target interval in alignment {}", i + 1);
    }
}

//-----------------------------------------------------------------------------

fn gbz_from(filename: &'static str) -> GBZ {
    let gbz_file = support::get_test_data(filename);
    let graph = serialize::load_from(&gbz_file);
    if let Err(err) = graph {
        panic!("Failed to load GBZ graph from {}: {}", gbz_file.display(), err);
    }
    graph.unwrap()
}

fn extract_subgraph(graph: &GBZ, nodes: &[usize]) -> Subgraph {
    let mut subgraph = Subgraph::new();
    for &node_id in nodes {
        let result = subgraph.add_node_from_gbz(graph, node_id);
        if let Err(err) = result {
            panic!("Failed to add node {} to the subgraph: {}", node_id, err);
        }
    }
    subgraph
}

fn clip_alignment<'a>(aln: &Alignment, subgraph: &Subgraph, sequence_len: Arc<impl Fn(usize) -> Option<usize> + 'a>) -> Vec<Alignment> {
    let clipped = aln.clip(subgraph, sequence_len);
    if let Err(err) = clipped {
        panic!("Failed to clip alignment {}: {}", aln.name, err);
    }
    clipped.unwrap()
}

fn add_fragment_indexes(truth: &mut [Vec<Alignment>], source: &[Alignment]) {
    for (i, fragments) in truth.iter_mut().enumerate() {
        if fragments.len() == 1 && fragments[0].seq_interval == source[i].seq_interval {
            // No clipping happened, so no fragment IDs.
            continue;
        }
        for (i, aln) in fragments.iter_mut().enumerate() {
            aln.optional.push(TypedField::Int([b'f', b'i'], (i + 1) as isize));
        }
    }
}

fn clip_and_check<'a>(
    alignments: &[Alignment], truth: &[Vec<Alignment>],
    subgraph: &Subgraph, sequence_len: Arc<impl Fn(usize) -> Option<usize> + 'a>
) {
    for (aln, true_fragments) in alignments.iter().zip(truth.iter()) {
        let clipped = clip_alignment(aln, &subgraph, sequence_len.clone());
        assert_eq!(clipped.len(), true_fragments.len(), "Wrong number of fragments for {}", aln.name);
        for (i, (fragment, true_fragment)) in clipped.iter().zip(true_fragments.iter()).enumerate() {
            assert_eq!(fragment.name, true_fragment.name, "Wrong name for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.seq_len, true_fragment.seq_len, "Wrong query length for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.seq_interval, true_fragment.seq_interval, "Wrong query interval for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.path, true_fragment.path, "Wrong target path for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.path_len, true_fragment.path_len, "Wrong target length for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.path_interval, true_fragment.path_interval, "Wrong target interval for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.difference, true_fragment.difference, "Wrong difference string for fragment {} of {}", i + 1, aln.name);
            assert_eq!(fragment.optional, true_fragment.optional, "Wrong optional fields for fragment {} of {}", i + 1, aln.name);
            // assert_eq!(fragment, true_fragment, "Wrong fragment {} for {}", i + 1, aln.name);
            // The statistics are currently for the entire alignment, while the truth may have them for the fragment.
        }
    }
}

fn check_fragment(aln: &Alignment, fragment_index: usize, clipped: &[Alignment], true_subpath: &[usize]) {
    assert!(fragment_index < clipped.len(), "Missing fragment {} of alignment {}", fragment_index + 1, aln.name);
    let fragment = &clipped[fragment_index];
    let fragment_path = fragment.target_path().unwrap();
    assert_eq!(fragment_path, true_subpath, "Wrong subpath for fragment {} of alignment {}", fragment_index + 1, aln.name);
}

fn clip_and_check_real<'a>(alignments: &[Alignment], subgraph: &Subgraph, sequence_len: Arc<impl Fn(usize) -> Option<usize> + 'a>) {
    for aln in alignments.iter() {
        let clipped = clip_alignment(aln, subgraph, sequence_len.clone());
        if aln.is_unaligned() {
            assert!(clipped.is_empty(), "Clipping created fragments for unaligned read {}", aln.name);
            continue;
        }

        // We have a target path, because clipping was successful.
        let target_path = aln.target_path().unwrap();
        let mut target_offset = 0;
        let mut curr: Option<Range<usize>> = None; // Interval of `target_path` in the subgraph.
        let mut fragment_index = 0;
        for (i, &handle) in target_path.iter().enumerate() {
            let node_len = sequence_len(handle).unwrap();
            let node_is_aligned =
                aln.path_interval.start < target_offset + node_len &&
                target_offset < aln.path_interval.end;
            target_offset += node_len;
            if node_is_aligned && subgraph.has_handle(handle) {
                if let Some(range) = curr.as_mut() {
                    range.end = i + 1;
                } else {
                    curr = Some(i..(i + 1));
                }
            } else {
                if let Some(range) = curr.take() {
                    let true_subpath = &target_path[range];
                    check_fragment(aln, fragment_index, &clipped, true_subpath);
                    fragment_index += 1;
                }
            }
        }
        if let Some(range) = curr.take() {
            let true_subpath = &target_path[range];
            check_fragment(aln, fragment_index, &clipped, true_subpath);
            fragment_index += 1;
        }

        assert_eq!(fragment_index, clipped.len(), "Clipped alignment {} has fewer than expected fragments", aln.name);
    }
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: Clipping to a subgraph.

#[test]
fn alignment_clip_unaligned() {
    let graph = gbz_from("example.gbz");
    let subgraph = extract_subgraph(&graph, &[11, 14]);
    let sequence_len = Arc::new(|handle| {
        let (node_id, _) = support::decode_node(handle);
        graph.sequence_len(node_id)
    });

    let alignments = vec![
        empty_alignment("empty", 0),
        empty_alignment("unaligned", 100),
    ];

    for aln in alignments.iter() {
        let clipped = clip_alignment(aln, &subgraph, sequence_len.clone());
        assert!(clipped.is_empty(), "Clipping created fragments for {}", aln.name);
    }
}

#[test]
fn alignment_clip_aligned() {
    let graph = gbz_from("example.gbz");
    let subgraph = extract_subgraph(&graph, &[11, 12, 14, 17]);
    let sequence_len = Arc::new(|handle| {
        let (node_id, _) = support::decode_node(handle);
        graph.sequence_len(node_id)
    });

    let alignments = vec![
        create_alignment(
            "in-subgraph",
            3, 0,
            vec![
                support::encode_node(11, Orientation::Forward),
                support::encode_node(12, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
            ], 0,
            vec![Difference::Match(3)],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "starts-with-insertion",
            5, 0,
            vec![
                support::encode_node(11, Orientation::Forward),
                support::encode_node(12, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
            ], 0,
            vec![Difference::Insertion(b"AC".to_vec()), Difference::Match(3)],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "starts-with-deletion",
            1, 0,
            vec![
                support::encode_node(11, Orientation::Forward),
                support::encode_node(12, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
            ], 0,
            vec![Difference::Deletion(2), Difference::Match(1)],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "outside-subgraph",
            3, 0,
            vec![
                support::encode_node(21, Orientation::Forward),
                support::encode_node(22, Orientation::Forward),
                support::encode_node(23, Orientation::Forward),
            ], 0,
            vec![Difference::Match(3)],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "in-two-parts",
            5, 0,
            vec![
                support::encode_node(11, Orientation::Forward),
                support::encode_node(12, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
                support::encode_node(15, Orientation::Forward),
                support::encode_node(17, Orientation::Forward),
            ], 0,
            vec![Difference::Match(3), Difference::Mismatch(b'A'), Difference::Match(1)],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "missing-ends",
            3, 0,
            vec![
                support::encode_node(13, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
                support::encode_node(16, Orientation::Forward),
            ], 0,
            vec![Difference::Match(3)],
            Some(sequence_len.clone())
        ),
    ];

    let mut truth = vec![
        vec![
            create_alignment(
                "in-subgraph",
                3, 0,
                vec![
                    support::encode_node(11, Orientation::Forward),
                    support::encode_node(12, Orientation::Forward),
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Match(3)],
                Some(sequence_len.clone())
            ),
        ],
        vec![
            create_alignment(
                "starts-with-insertion",
                5, 0,
                vec![
                    support::encode_node(11, Orientation::Forward),
                    support::encode_node(12, Orientation::Forward),
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Insertion(b"AC".to_vec()), Difference::Match(3)],
                Some(sequence_len.clone())
            ),
        ],
        vec![
            create_alignment(
                "starts-with-deletion",
                1, 0,
                vec![
                    support::encode_node(11, Orientation::Forward),
                    support::encode_node(12, Orientation::Forward),
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Deletion(2), Difference::Match(1)],
                Some(sequence_len.clone())
            ),
        ],
        vec![], // "outside-subgraph"
        vec![
            create_alignment(
                "in-two-parts",
                5, 0,
                vec![
                    support::encode_node(11, Orientation::Forward),
                    support::encode_node(12, Orientation::Forward),
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Match(3)],
                Some(sequence_len.clone())
            ),
            create_alignment(
                "in-two-parts",
                5, 4,
                vec![
                    support::encode_node(17, Orientation::Forward),
                ], 0,
                vec![Difference::Match(1)],
                Some(sequence_len.clone())
            ),
        ],
        vec![
            create_alignment(
                "missing-ends",
                3, 1,
                vec![
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Match(1)],
                Some(sequence_len.clone())
            ),
        ],
    ];
    add_fragment_indexes(&mut truth, &alignments);

    clip_and_check(&alignments, &truth, &subgraph, sequence_len);
}

// with longer nodes

#[test]
fn alignment_clip_long_nodes() {
    let graph = gbz_from("example.gbz");
    let subgraph = extract_subgraph(&graph, &[11, 12, 14, 17]);

    // Alignment iteration/clipping does not use the graph. We can therefore pretend
    // that the nodes are longer than they actually are.
    let mut node_len = HashMap::new();
    for node_id in graph.node_iter() {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let handle = support::encode_node(node_id, orientation);
            node_len.insert(handle, node_id - 10);
        }
    }
    let sequence_len = Arc::new(|handle| {
        node_len.get(&handle).cloned()
    });

    let alignments = vec![
        create_alignment(
            "upper-path",
            19, 2,
            vec![
                support::encode_node(11, Orientation::Forward),
                support::encode_node(12, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
                support::encode_node(15, Orientation::Forward),
                support::encode_node(17, Orientation::Forward),
            ], 0,
            vec![
                Difference::Match(4), Difference::Deletion(2),
                Difference::Match(6), Difference::Insertion(b"GA".to_vec()), // Insertion at the end of node 15.
                Difference::Match(5),
            ],
            Some(sequence_len.clone())
        ),
        create_alignment(
            "lower-path",
            21, 1,
            vec![
                support::encode_node(13, Orientation::Forward),
                support::encode_node(14, Orientation::Forward),
                support::encode_node(16, Orientation::Forward),
                support::encode_node(17, Orientation::Forward),
            ], 2,
            vec![
                Difference::Match(5), Difference::Insertion(b"GA".to_vec()), // Insertion at the end of node 14.
                Difference::Match(7), Difference::Mismatch(b'T'),
                Difference::Match(5),
            ],
            Some(sequence_len.clone())
        ),
    ];

    let mut truth = vec![
        vec![
            create_alignment(
                "upper-path",
                19, 2,
                vec![
                    support::encode_node(11, Orientation::Forward),
                    support::encode_node(12, Orientation::Forward),
                    support::encode_node(14, Orientation::Forward),
                ], 0,
                vec![Difference::Match(4), Difference::Deletion(2), Difference::Match(1)],
                Some(sequence_len.clone())
            ),
            create_alignment(
                "upper-path",
                19, 14,
                vec![support::encode_node(17, Orientation::Forward)], 0,
                vec![Difference::Match(5)],
                Some(sequence_len.clone())
            ),
        ],
        vec![
            create_alignment(
                "lower-path",
                21, 2,
                vec![support::encode_node(14, Orientation::Forward)], 0,
                vec![Difference::Match(4), Difference::Insertion(b"GA".to_vec())],
                Some(sequence_len.clone())
            ),
            create_alignment(
                "lower-path",
                21, 14,
                vec![support::encode_node(17, Orientation::Forward)], 0,
                vec![Difference::Match(1), Difference::Mismatch(b'T'), Difference::Match(5)],
                Some(sequence_len.clone())
            ),
        ],
    ];
    add_fragment_indexes(&mut truth, &alignments);

    clip_and_check(&alignments, &truth, &subgraph, sequence_len);
}

// real reads
#[test]
fn alignment_clip_real() {
    // Get the alignments.
    let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
    let alignments = parse_alignments(&gaf_file, false);

    // Load the graph.
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    let sequence_len = Arc::new(|handle| {
        let (node_id, _) = support::decode_node(handle);
        graph.sequence_len(node_id)
    });

    // Subgraph queries for the default context around 10 random nodes.
    let node_ids: Vec<usize> = graph.node_iter().collect();
    let mut queries = Vec::new();
    let mut rng = rand::rng();
    for _ in 0..10 {
        let node_id = node_ids[rng.random_range(0..node_ids.len())];
        let mut query = BTreeSet::new();
        query.insert(node_id);
        queries.push(query);
    }

    // Execute the queries and try clipping every alignment to the resulting subgraph.
    for query in queries.iter() {
        let mut subgraph = Subgraph::new();
        let result = subgraph.around_nodes(
            GraphReference::Gbz(&graph), query, SubgraphQuery::DEFAULT_CONTEXT
        );
        if let Err(err) = result {
            panic!("Failed to extract the subgraph around nodes {:?}: {}", query, err);
        }
        clip_and_check_real(&alignments, &subgraph, sequence_len.clone());
    }
}

//-----------------------------------------------------------------------------

