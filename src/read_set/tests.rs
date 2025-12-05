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

fn validate_read_set(read_set: &ReadSet, subgraph: &Subgraph, all_reads: &[Alignment], name: &str, output: AlignmentOutput) {
    // Select the alignments that should be included in the read set.
    let mut expected_reads = Vec::new();
    let mut unclipped = 0;
    for aln in all_reads {
        let target_path = aln.target_path().unwrap();
        if target_path.is_empty() {
            continue;
        }
        let mut take = false;
        for handle in target_path {
            if subgraph.has_handle(*handle) {
                take = true;
            } else if output == AlignmentOutput::Contained {
                take = false;
                break;
            }
        }
        if take {
            let mut aln = aln.clone();
            aln.optional.clear(); // We do not store unknown optional fields for the moment.
            if output == AlignmentOutput::Clipped {
                // Because we ran the test first with overlapping reads, we can assume that
                // the read set contains the necessary node records for clipping.
                let sequence_len = Arc::new(|handle| {
                    let record = read_set.nodes.get(&handle)?;
                    Some(record.sequence().len())
                });
                let clipped = aln.clip(subgraph, sequence_len);
                if let Err(err) = clipped {
                    panic!("Failed to clip an overlapping read for {}: {}", name, err);
                }
                let clipped = clipped.unwrap();
                for fragment in clipped.into_iter() {
                    expected_reads.push(fragment);
                }
            } else {
                expected_reads.push(aln);
            }
            unclipped += 1;
        }
    }
    assert_eq!(read_set.len(), expected_reads.len(), "Wrong number of {} reads in read set for {}", output, name);
    assert_eq!(read_set.unclipped(), unclipped, "Wrong number of unclipped {} reads in read set for {}", output, name);

    // GAF-base construction and read set extraction maintain order.
    for (i, (aln, truth)) in read_set.iter().zip(expected_reads.iter()).enumerate() {
        assert_eq!(aln, truth, "Wrong {} read {} for {}", output, i, name);
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
        let gaf_line = aln.to_gaf(&target_sequence);
        expected_gaf.extend_from_slice(&gaf_line);
        expected_gaf.push(b'\n');
    }

    // Check the GAF output.
    let mut gaf_output = Vec::new();
    let result = read_set.to_gaf(&mut gaf_output);
    assert!(result.is_ok(), "Failed to output {} GAF for {}: {}", output, name, result.unwrap_err());
    assert_eq!(gaf_output.len(), expected_gaf.len(), "Wrong {} GAF output length for {}", output, name);
    assert_eq!(gaf_output, expected_gaf, "Wrong {} GAF output for {}", output, name);
}

//-----------------------------------------------------------------------------

// Tests for ReadSet.

#[test]
fn read_set_default() {
    let read_set = ReadSet::default();
    assert_eq!(read_set.len(), 0, "Default ReadSet has non-zero length");
    assert_eq!(read_set.unclipped(), 0, "Default ReadSet has non-zero unclipped count");
    assert_eq!(read_set.blocks(), 0, "Default ReadSet has non-zero block count");
    assert_eq!(read_set.clusters(), 0, "Default ReadSet has non-zero cluster count");
    assert!(read_set.iter().next().is_none(), "Default ReadSet iterator is not empty");
}

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

        for output in [AlignmentOutput::Overlapping, AlignmentOutput::Clipped, AlignmentOutput::Contained] {
            let read_set = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &gaf_base, output);
            assert!(read_set.is_ok(), "Failed to extract {} reads for {}: {}", output, name, read_set.unwrap_err());
            let read_set = read_set.unwrap();
            validate_read_set(&read_set, &subgraph, &all_reads, &name, output);
        }
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

        for output in [AlignmentOutput::Overlapping, AlignmentOutput::Clipped, AlignmentOutput::Contained] {
            let read_set = ReadSet::new(GraphReference::Db(&mut graph), &subgraph, &gaf_base, output);
            assert!(read_set.is_ok(), "Failed to extract {} reads for {}: {}", output, name, read_set.unwrap_err());
            let read_set = read_set.unwrap();
            validate_read_set(&read_set, &subgraph, &all_reads, &name, output);
        }
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

#[test]
fn read_set_from_rows() {
    // Load GBZ.
    let graph = internal::load_gaf_base_gbz();

    // Build and open GAF-base.
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt");
    let gaf_base = internal::open_gaf_base(&gaf_base_file);

    // Parse the reads as a source of truth.
    let mut all_reads = internal::load_gaf_base_reads(false);

    let chunk_sizes = vec![1, 2, 5];
    for chunk_size in chunk_sizes {
        let mut found_alns = 0;
        let mut found_blocks = 0;
        let mut rowid = 1; // SQLite row ids start from 1.
        while found_alns < gaf_base.alignments() {
            let range = rowid..(rowid + chunk_size);
            let read_set = ReadSet::from_rows(&gaf_base, range.clone(), &graph);
            assert!(read_set.is_ok(), "Failed to extract reads from rows {}..{}: {}", range.start, range.end, read_set.unwrap_err());
            let read_set = read_set.unwrap();

            assert_eq!(read_set.len(), read_set.unclipped(), "Extracted clipped alignments from rows {}..{}", range.start, range.end);
            assert!(found_alns + read_set.len() <= all_reads.len(), "Extracted too many alignments from rows {}..{}", range.start, range.end);
            for (i, aln) in read_set.iter().enumerate() {
                let truth = &mut all_reads[found_alns + i];
                truth.optional.clear(); // We do not store unknown optional fields for the moment.
                assert_eq!(aln, truth, "Wrong read {} from rows {}..{}", i, range.start, range.end);
            }

            found_alns += read_set.len();
            found_blocks += read_set.blocks();
            rowid += chunk_size;
        }
        assert_eq!(found_alns, gaf_base.alignments(), "Wrong total number of alignments with chunk size {}", chunk_size);
        assert_eq!(found_blocks, gaf_base.blocks(), "Wrong total number of extracted with chunk size {}", chunk_size);
    }

    // Cleanup.
    drop(gaf_base);
    let _ = std::fs::remove_file(&gaf_base_file);
}

//-----------------------------------------------------------------------------
