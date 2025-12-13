use super::*;

use crate::utils;

use gbwt::Metadata;

use simple_sds::serialize;

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

    // Check the type and the value.
    match field.clone() {
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

    if let TypedField::Float(_, _) = field {
        // We skip floats, as they can have multiple representations for the same value.
        return;
    }

    // Conversion to a string.
    let truth_str = String::from_utf8_lossy(seq);
    assert_eq!(field.to_string(), truth_str, "Wrong string representation for {}", name);

    // Conversion to bytes.
    let mut buffer = Vec::new();
    let mut truth_bytes = seq.to_vec();
    field.append_to(&mut buffer, false);
    assert_eq!(buffer, truth_bytes, "Wrong bytes for {}", name);

    // Append as a new field.
    truth_bytes.push(b'\t');
    truth_bytes.extend_from_slice(seq);
    field.append_to(&mut buffer, true);
    assert_eq!(buffer, truth_bytes, "Wrong bytes for {} when appended as a new field", name);
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

// Tests for GAF reading.

#[test]
fn gaf_headers() {
    // Filename and expected number of header lines.
    let test_cases = vec![
        ("empty.gaf", 0),
        ("good.gaf", 2),
        ("micb-kir3dl1_HG003.gaf", 2),
        ("no_header.gaf", 0)
    ];

    for (filename, expected_lines) in test_cases {
        let filepath = utils::get_test_data(filename);
        let mut reader = utils::open_file(&filepath)
            .expect(&format!("Failed to open test GAF file {}", filename));
        let headers = read_gaf_header_lines(&mut reader);
        assert!(headers.is_ok(), "Failed to read GAF headers from file {}: {}", filename, headers.unwrap_err());
        let headers = headers.unwrap();
        assert_eq!(headers.len(), expected_lines, "Wrong number of GAF header lines in file {}", filename);
        for (i, line) in headers.iter().enumerate() {
            let is_header = is_gaf_header_line(line.as_bytes());
            assert!(is_header, "Line {} in file {} is not a valid GAF header line", i + 1, filename);
            let last = line.as_bytes().last();
            assert_ne!(last, Some(&b'\n'), "Line {} in file {} has a trailing newline", i + 1, filename);
        }
        let peek = peek_gaf_header_line(&mut reader);
        if let Ok(is_header) = peek {
            assert!(!is_header, "Additional header lines found in file {}", filename);
        }
    }
}

//-----------------------------------------------------------------------------
