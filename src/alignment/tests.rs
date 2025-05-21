use super::*;

use crate::utils;

use simple_sds::ops::BitVec;
use simple_sds::serialize;

use std::io::BufRead;
use std::path::PathBuf;

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
    let reader = utils::open_file(filename);
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
    let filename = utils::get_test_data("good.gaf");
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

// FIXME: In blocks.
// Tests for `Alignment`: encoding and decoding.

//-----------------------------------------------------------------------------

// FIXME: Use blocks.
// Tests for `Alignment`: integration.
// This is the haplotype sampling test case from vg.
/*
fn integration_test(gaf_file: &'static str, gbwt_file: &'static str) {
    // Parse the alignments.
    let gaf_file = utils::get_test_data(gaf_file);
    let mut alignments = parse_alignments(&gaf_file, false);
    assert_eq!(alignments.len(), 12439, "Unexpected number of parsed alignments");

    // Load the GBWT index and prepare the structures.
    let gbwt_file = utils::get_test_data(gbwt_file);
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
    let alphabet = b"F:#,";
    let dictionary = [(1, 1), (2, 1), (3, 2)];
    let quality_encoder = QualityEncoder::new(alphabet, &dictionary).unwrap();
    for (i, aln) in alignments.iter().enumerate() {
        check_encode_decode(aln, prefix, &quality_encoder, i + 1);
    }
}

#[test]
fn alignment_real() {
    integration_test("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt");
}

#[test]
fn alignment_real_gzipped() {
    integration_test("micb-kir3dl1_HG003.gaf.gz", "micb-kir3dl1_HG003.gbwt");
}*/

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
