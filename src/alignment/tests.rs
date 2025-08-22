use super::*;

use crate::utils;

use simple_sds::serialize;

use std::io::BufRead;
use std::path::PathBuf;

//-----------------------------------------------------------------------------

// Tests for `Alignment`: parsing.

// Creates an alignment object for an unaligned sequence with no optional fields.
fn empty_alignment(name: &str, seq_len: usize) -> Alignment {
    let mut result = Alignment::default();
    result.name = String::from(name);
    result.seq_len = seq_len;
    result
}

// Returns the GAF lines from the given file.
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
fn parse_alignments(filename: &PathBuf, expect_error: bool) -> Vec<Alignment> {
    let lines = read_gaf_lines(filename);

    let mut result = Vec::new();
    for (line_num, line) in lines.iter().enumerate() {
        let alignment = Alignment::from_gaf(line);
        if expect_error {
            assert!(alignment.is_err(), "Expected a parse error on line {}", line_num + 1);
        } else {
            assert!(alignment.is_ok(), "Failed to parse line {}: {}", line_num + 1, alignment.err().unwrap());
            result.push(alignment.unwrap());
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
    let filename = utils::get_test_data("good.gaf");
    let alignments = parse_alignments(&filename, false);
    assert_eq!(alignments.len(), truth.len(), "Wrong number of alignments in the test file");
    for i in 0..truth.len() {
        check_alignment(&alignments[i], &truth[i], i + 1, false, false);
    }

    // These alignments will survive the round-trip from/to GAF intact.
    // 1 is on the reverse strand, while 7 has non-canonical empty fields.
    let round_trip = vec![
        0, 2, 3, 4, 5, 6,
    ];
    let gaf_lines = read_gaf_lines(&filename);
    let target_sequence = b"CCCCCCCCCCCCCCCCCCCCAGAAAAAAAAAAAAAAAAAAAAGATTACATGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTT".to_vec();
    for index in round_trip {
        let line = alignments[index].to_gaf(&target_sequence);
        assert_eq!(line, gaf_lines[index], "Line {} did not survive the round trip", index + 1);
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
    assert!(!default.is_perfect(), "An empty alignment should not be perfect");

    assert!(default.name.is_empty(), "A default alignment should not have a name");
    assert_eq!(default.seq_len, 0, "A default alignment should have a sequence length of 0");
    assert!(default.seq_interval.is_empty(), "A default alignment should have an empty sequence interval");
    assert!(default.has_target_path(), "A default alignment should have a target path");
    assert!(default.target_path().unwrap().is_empty(), "A default alignment should have an empty target path");
    assert_eq!(default.path_len, 0, "A default alignment should have a target path length of 0");
    assert!(default.path_interval.is_empty(), "A default alignment should have an empty target path interval");
    assert_eq!(default.matches, 0, "A default alignment should have 0 matches");
    assert_eq!(default.edits, 0, "A default alignment should have 0 edits");
    assert!(default.mapq.is_none(), "A default alignment should not have a mapping quality");
    assert!(default.score.is_none(), "A default alignment should not have a score");
    assert!(default.base_quality.is_empty(), "A default alignment should not have base quality values");
    assert!(default.difference.is_empty(), "A default alignment should not have a difference string");
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
    assert!(aln.is_perfect(), "The alignment should be perfect");

    // Various ways of making the alignment imperfect.
    {
        let mut wrong_start = aln.clone();
        wrong_start.seq_interval.start += 1;
        wrong_start.path_interval.start += 1;
        wrong_start.matches -= 1;
        assert!(!wrong_start.is_perfect(), "An alignment starting after query start should not be perfect");
    }
    {
        let mut wrong_end = aln.clone();
        wrong_end.seq_interval.end -= 1;
        wrong_end.path_interval.end -= 1;
        wrong_end.matches -= 1;
        assert!(!wrong_end.is_perfect(), "An alignment ending before query end should not be perfect");
    }
    {
        let mut has_edits = aln.clone();
        has_edits.matches -= 1;
        has_edits.edits += 1;
        assert!(!has_edits.is_perfect(), "An alignment with edits should not be perfect");
    }

    // We have a GBWT start position.
    assert!(!aln.has_target_path(), "GBWT start position should not be a target path");
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
    assert_eq!(aln.target_path(), Some(correct_path.as_slice()), "Target path should match the set path");
    assert_eq!(aln.min_handle(), min_handle, "Minimum handle should match the set path");
    assert_eq!(aln.max_handle(), max_handle, "Maximum handle should match the set path");

    // Setting it again should not change the existing target path.
    let wrong_path = vec![
        support::encode_node(1, Orientation::Forward),
        support::encode_node(3, Orientation::Forward),
    ];
    aln.set_target_path(wrong_path);
    assert!(aln.has_target_path(), "Target path should still exist after setting it again");
    assert_eq!(aln.target_path(), Some(correct_path.as_slice()), "Target path should match the existing path");
    assert_eq!(aln.min_handle(), min_handle, "Minimum handle should match the existing path");
    assert_eq!(aln.max_handle(), max_handle, "Maximum handle should match the existing path");
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
        assert_eq!(*aln, block[i], "Wrong alignment {} in the decoded block {}..{}", range.start + i, range.start, range.end);
    }
}

//-----------------------------------------------------------------------------

// Tests for `Alignment`: integration.
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

fn check_to_bytes(ops: &[Difference], truth: &[u8], target_sequence: &[u8], name: &str) {
    let result = Difference::to_bytes(ops, target_sequence);
    assert_eq!(result, truth, "Wrong byte representation for {}", name);
}

#[test]
fn difference_empty() {
    let name = "empty";
    let seq = b"";
    let truth = [];
    check_difference(seq, &truth, name, false);
    let target_sequence = b"";
    check_to_bytes(&truth, seq, target_sequence, name);
}

#[test]
fn difference_single() {
    {
        let name = "matching sequence";
        let seq = b"=ACGT";
        let truth = [ Difference::Match(4) ];
        check_difference(seq, &truth, name, false);
        let expected = b":4"; // Converted to match length.
        let target_sequence = b"ACGT";
        check_to_bytes(&truth, expected, target_sequence, name);
    }
    {
        let name = "match length";
        let seq = b":5";
        let truth = [ Difference::Match(5) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"GATTA";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "mismatch";
        let seq = b"*AC";
        let truth = [ Difference::Mismatch(b'C') ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"A";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "insertion";
        let seq = b"+ACGT";
        let truth = [ Difference::Insertion(b"ACGT".to_vec()) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "deletion";
        let seq = b"-ACGT";
        let truth = [ Difference::Deletion(4) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"ACGT";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
}

#[test]
fn difference_multi() {
    // All of these are from real reads mapped with Giraffe.
    // Target sequences only contain the relevant bases.
    {
        let name = "ERR3239454.129477191 fragment 2";
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
        check_difference(seq, &truth, name, false);
        let mut target_sequence = vec![b'-'; 151];
        target_sequence[24] = b'C';
        target_sequence[103] = b'C';
        target_sequence[139] = b'T';
        target_sequence[142] = b'T';
        target_sequence[144] = b'T';
        check_to_bytes(&truth, seq, &target_sequence, name);
    }
    {
        let name = "ERR3239454.11251898 fragment 1";
        let seq = b":129+TTGAGGGGGTATAAGAGGTCG";
        let truth = [
            Difference::Match(129),
            Difference::Insertion(b"TTGAGGGGGTATAAGAGGTCG".to_vec())
        ];
        check_difference(seq, &truth, name, false);
        let target_sequence = vec![b'-', 129];
        check_to_bytes(&truth, seq, &target_sequence, name);
    }
    {
        let name = "ERR3239454.97848632 fragment 1";
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
        check_difference(seq, &truth, name, false);
        let mut target_sequence = vec![b'-'; 150];
        target_sequence[111] = b'G';
        target_sequence[132] = b'C';
        target_sequence[139] = b'G';
        target_sequence[142] = b'G';
        target_sequence[146] = b'G';
        check_to_bytes(&truth, seq, &target_sequence, name);
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

