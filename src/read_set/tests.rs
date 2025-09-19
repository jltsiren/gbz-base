use super::*;

use crate::{Subgraph, SubgraphQuery, HaplotypeOutput};
use crate::internal;

//-----------------------------------------------------------------------------

fn gaf_base_subgraph_queries() -> Vec<(SubgraphQuery, String)> {
    vec![
        (SubgraphQuery::nodes(vec![]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("empty")),
        (SubgraphQuery::nodes(vec![1000]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("single component")),
        (SubgraphQuery::nodes(vec![500, 1500]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("two components")),
    ]
}

fn validate_read_set(read_set: &ReadSet, subgraph: &Subgraph, all_reads: &[Alignment], name: &str, contained: bool) {
    // Select the alignments that should be included in the read set.
    let mut expected_reads = Vec::new();
    for aln in all_reads {
        let target_path = aln.target_path().unwrap();
        if target_path.is_empty() {
            continue;
        }
        let mut take = false;
        for handle in target_path {
            if subgraph.has_handle(*handle) {
                take = true;
            } else if contained {
                take = false;
                break;
            }
        }
        if take {
            let mut aln = aln.clone();
            aln.optional.clear(); // We do not store unknown optional fields for the moment.
            expected_reads.push(aln);
        }
    }
    assert_eq!(read_set.len(), expected_reads.len(), "Wrong number of reads in read set for {}", name);

    // GAF-base construction and read set extraction maintain order.
    for (i, (aln, truth)) in read_set.iter().zip(expected_reads.iter()).enumerate() {
        assert_eq!(aln, truth, "Wrong read {} for {}", i, name);
    }

    // Generate the expected GAF output.
    let mut expected_gaf = Vec::new();
    for aln in expected_reads {
        let mut target_sequence = Vec::new();
        let target_path = aln.target_path().unwrap();
        for handle in target_path {
            let record = read_set.nodes.get(handle);
            assert!(record.is_some(), "Missing GBZ record for handle {} in {}", handle, name);
            let record = record.unwrap();
            target_sequence.extend_from_slice(record.sequence());
        }
        target_sequence = target_sequence[aln.path_interval.clone()].to_vec();
        let gaf_line = aln.to_gaf(&target_sequence);
        expected_gaf.extend_from_slice(&gaf_line);
        expected_gaf.push(b'\n');
    }

    // Check the GAF output.
    let mut gaf_output = Vec::new();
    let result = read_set.to_gaf(&mut gaf_output);
    assert!(result.is_ok(), "Failed to output GAF for {}: {}", name, result.unwrap_err());
    assert_eq!(gaf_output.len(), expected_gaf.len(), "Wrong GAF output length for {}", name);
    assert_eq!(gaf_output, expected_gaf, "Wrong GAF output for {}", name);
}

//-----------------------------------------------------------------------------

// Tests for ReadSet.

fn test_read_set_gbz(gbwt_part: &'static str) {
    // Load GBZ.
    let graph = internal::load_gaf_base_gbz();

    // Build and open GAF-base.
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", gbwt_part);
    let gaf_base = internal::open_gaf_base(&gaf_base_file);

    // Parse the reads as a source of truth.
    let all_reads = internal::load_gaf_base_reads(false);

    let queries = gaf_base_subgraph_queries();
    for (query, name) in queries {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_gbz(&graph, None, None, &query);
        assert!(result.is_ok(), "Failed to extract subgraph for {}: {}", name, result.unwrap_err());

        let overlapping = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &gaf_base, false);
        assert!(overlapping.is_ok(), "Failed to extract overlapping reads for {}: {}", name, overlapping.unwrap_err());
        let overlapping = overlapping.unwrap();
        validate_read_set(&overlapping, &subgraph, &all_reads, &name, false);

        let contained = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &gaf_base, true);
        assert!(contained.is_ok(), "Failed to extract fully contained reads for {}: {}", name, contained.unwrap_err());
        let contained = contained.unwrap();
        validate_read_set(&contained, &subgraph, &all_reads, &name, true);
    }

    // Cleanup.
    drop(gaf_base);
    let _ = std::fs::remove_file(&gaf_base_file);
}

#[test]
fn read_set_gbz_unidirectional() {
    test_read_set_gbz("micb-kir3dl1_HG003.gbwt");
}

#[test]
fn read_set_gbz_bidirectional() {
    test_read_set_gbz("bidirectional.gbwt");
}

fn test_read_set_db(gbwt_part: &'static str) {
    // Build and open GAF-base.
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", gbwt_part);
    let gaf_base = internal::open_gaf_base(&gaf_base_file);

    // Open GBZ-base.
    let (gbz_base, gbz_base_file) = internal::create_gaf_base_db();
    let mut graph = internal::create_graph_interface(&gbz_base);

    // Parse the reads as a source of truth.
    let all_reads = internal::load_gaf_base_reads(false);

    let queries = gaf_base_subgraph_queries();
    for (query, name) in queries {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_db(&mut graph, &query);
        assert!(result.is_ok(), "Failed to extract subgraph for {}: {}", name, result.unwrap_err());

        let overlapping = ReadSet::new(GraphReference::Db(&mut graph), &subgraph, &gaf_base, false);
        assert!(overlapping.is_ok(), "Failed to extract overlapping reads for {}: {}", name, overlapping.unwrap_err());
        let overlapping = overlapping.unwrap();
        validate_read_set(&overlapping, &subgraph, &all_reads, &name, false);

        let contained = ReadSet::new(GraphReference::Db(&mut graph), &subgraph, &gaf_base, true);
        assert!(contained.is_ok(), "Failed to extract fully contained reads for {}: {}", name, contained.unwrap_err());
        let contained = contained.unwrap();
        validate_read_set(&contained, &subgraph, &all_reads, &name, true);
    }

    // Cleanup.
    drop(graph);
    drop(gbz_base);
    let _ = std::fs::remove_file(&gbz_base_file);
    drop(gaf_base);
    let _ = std::fs::remove_file(&gaf_base_file);
}

#[test]
fn read_set_db_unidirectional() {
    test_read_set_db("micb-kir3dl1_HG003.gbwt");
}

#[test]
fn read_set_db_bidirectional() {
    test_read_set_db("bidirectional.gbwt");
}

//-----------------------------------------------------------------------------
