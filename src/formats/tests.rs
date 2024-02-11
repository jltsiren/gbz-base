use super::*;

use gbwt::Metadata;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

#[test]
fn walk_metadata_interval() {
    let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
    let interval = 123..456;
    let metadata = WalkMetadata::path_interval(&path_name, interval.clone(), None);

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
    let metadata = WalkMetadata::anonymous(haplotype, contig, len, None);

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
    let mut metadata = WalkMetadata::anonymous(2, "chr19", 456, None);
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
        let metadata = WalkMetadata::path_interval(&path_name, interval, None);
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
        let metadata = WalkMetadata::path_interval(&path_name, interval, Some(42));
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
        let mut metadata = WalkMetadata::path_interval(&path_name, interval, None);
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
        let mut metadata = WalkMetadata::path_interval(&path_name, interval, Some(42));
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let mut output: Vec<u8> = Vec::new();
        let result = write_gfa_walk(&path, &metadata, &mut output);
        assert!(result.is_ok(), "Failed to write a GFA walk line with weight and a CIGAR string");
        let correct = b"W\tGRCh38\t0\tchr13\t1000\t1123\t>21<22>23\tWT:i:42\tCG:Z:10M2D5M\n";
        assert_eq!(&output, correct, "Wrong GFA walk line with weight and a CIGAR string");
    }
}

//-----------------------------------------------------------------------------

// TODO: More tests for JSON writing?

#[test]
fn write_json_path() {
    {
        // No weight and no CIGAR string.
        let path = [42, 45, 46];
        let path_name = FullPathName::haplotype("GRCh38", "chr13", 0, 1000);
        let interval = 0..123;
        let metadata = WalkMetadata::path_interval(&path_name, interval, None);
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
        let metadata = WalkMetadata::path_interval(&path_name, interval, Some(42));
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
        let mut metadata = WalkMetadata::path_interval(&path_name, interval, None);
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
        let mut metadata = WalkMetadata::path_interval(&path_name, interval, Some(42));
        metadata.add_cigar(Some(String::from("10M2D5M")));
        let json_value = json_path(&path, &metadata);
        let result = json_value.to_string();
        let correct = "{\"name\": \"GRCh38#0#chr13[1000-1123]\", \"weight\": 42, \"cigar\": \"10M2D5M\", \"path\": [{\"id\": \"21\", \"is_reverse\": false}, {\"id\": \"22\", \"is_reverse\": true}, {\"id\": \"23\", \"is_reverse\": false}]}";
        assert_eq!(result, correct, "Wrong JSON for a path with weight and a CIGAR string");
    }
}

//-----------------------------------------------------------------------------
