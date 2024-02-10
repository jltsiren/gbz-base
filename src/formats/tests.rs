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

// TODO: GFA writing

//-----------------------------------------------------------------------------

// TODO: JSON writing

//-----------------------------------------------------------------------------
