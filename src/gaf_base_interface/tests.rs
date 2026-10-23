use super::*;

use crate::{GAFBaseParams, formats, internal};

use gbz::{GBZ, GBWT, Orientation, ENDMARKER};
use gbz::bwt::BWT;
use gbz::support;

use simple_sds::serialize;

use std::collections::HashSet;
use std::fs;

//-----------------------------------------------------------------------------

#[test]
fn get_record() {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);

    // We get sequences from an in-memory GBZ graph.
    let gbz: GBZ = internal::load_gaf_base_gbz();
    let mut get_sequence = |handle: usize| -> Result<Vec<u8>> {
        let (node, orientation) = support::decode_node(handle);
        let seq = gbz.sequence(node).ok_or_else(||
            Error::invalid_data(format!("Node {} does not exist", node))
        )?;
        if orientation == Orientation::Reverse {
            Ok(support::reverse_complement(seq))
        } else {
            Ok(seq.to_vec())
        }
    };

    // A GBZ record combines a GBWT record and a sequence.
    let gbwt_file = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
    let gbwt = serialize::load_from(&gbwt_file);
    assert!(gbwt.is_ok(), "Failed to load GBWT index: {}", gbwt.unwrap_err());
    let gbwt: GBWT = gbwt.unwrap();
    let bwt: &BWT = gbwt.as_ref();

    // Compare all records against the truth.
    for truth_record in bwt.iter() {
        if truth_record.id() == ENDMARKER {
            continue;
        }
        let handle = gbwt.record_to_node(truth_record.id());
        let gbz_record = interface.record(handle, &mut get_sequence);
        assert!(gbz_record.is_ok(), "Failed to get GBZ record for handle {}", handle);
        let gbz_record = gbz_record.unwrap();
        assert!(gbz_record.is_some(), "No GBZ record found for handle {}", handle);
        let gbz_record = gbz_record.unwrap();
        let truth_sequence = get_sequence(handle).unwrap();
        assert_eq!(gbz_record.sequence(), &truth_sequence, "Wrong sequence for handle {}", handle);
        // to_gbwt_record() would use the handle as the id.
        let record = gbz_record.to_gbwt_record_with_id(truth_record.id());
        assert_eq!(record, truth_record, "Wrong GBWT record for handle {}", handle);
    }

    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

//-----------------------------------------------------------------------------

fn run_direct_decompression(gbz: Option<GBZ>) {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);

    let truth = formats::read_gaf_file(utils::get_test_data("micb-kir3dl1_HG003.gaf"));
    assert!(truth.is_ok(), "Failed to read GAF file: {}", truth.unwrap_err());
    let (_, truth) = truth.unwrap();

    // Use a direct database connection for querying the rows in the Alignments table.
    let direct_connection = internal::open_connection(&gaf_base_file).unwrap();
    let statement = direct_connection.prepare("SELECT * FROM Alignments ORDER BY id");
    assert!(statement.is_ok(), "Failed to prepare statement: {}", statement.unwrap_err());
    let mut statement = statement.unwrap();

    // Return the sequence from the GBZ if available, or an empty sequence if not.
    let mut get_sequence = |handle: usize| -> Result<Vec<u8>> {
        if let Some(gbz) = gbz.as_ref() {
            let (node, orientation) = support::decode_node(handle);
            let seq = gbz.sequence(node).ok_or_else(||
                Error::invalid_data(format!("Node {} does not exist", node))
            )?;
            if orientation == Orientation::Reverse {
                Ok(support::reverse_complement(seq))
            } else {
                Ok(seq.to_vec())
            }
        } else {
            Ok(Vec::new())
        }
    };

    {
        let rows  = statement.query(());
        let mut rows = match rows {
            Ok(rows) => rows,
            Err(err) => panic!("Failed to query Alignments: {}", err),
        };
        let mut line: usize = 0;
        let mut cache = HashMap::new();
        loop {
            let row = rows.next();
            assert!(row.is_ok(), "Failed to get next row: {}", row.unwrap_err());
            let row = row.unwrap();
            if row.is_none() {
                break;
            }
            // We skip the id field.
            let block = GAFBaseInterface::decompress_block(row.unwrap(), 1);
            assert!(block.is_ok(), "Failed to decompress block: {}", block.unwrap_err());
            let mut block = block.unwrap();
            for aln in block.iter_mut() {
                let result = interface.set_target_path(aln, &mut cache, &mut get_sequence);
                assert!(result.is_ok(), "Failed to set target path: {}", result.unwrap_err());
                if gbz.is_none() {
                    // We must set target path length explicitly when `get_sequence` does not return actual sequences.
                    aln.set_target_path_len(truth[line].path_len);
                }
                assert_eq!(*aln, truth[line], "Wrong alignment at line {}", line);
                // Perform assertions or modifications on each alignment here.
                line += 1;
            }
        }
    }

    drop(statement);
    drop(direct_connection);
    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

#[test]
fn direct_decompression() {
    let gbz = None;
    run_direct_decompression(gbz);
}

#[test]
fn direct_decompression_with_gbz() {
    let gbz = Some(internal::load_gaf_base_gbz());
    run_direct_decompression(gbz);
}

//-----------------------------------------------------------------------------

// Checks that the row id / offset tags are correct and then removes them.
fn check_ids(aln: &mut Alignment, row_id: usize, row_offset: usize) {
    let row_id = row_id as isize;
    let row_offset = row_offset as isize;

    let mut found_ri = false;
    let mut found_ro = false;
    for tag in aln.optional.iter() {
        if let TypedField::Int(tag_key, value) = tag {
            if tag_key == &[b'r', b'i'] {
                assert!(!found_ri, "Duplicate row id for ({}, {})", row_id, row_offset);
                found_ri = true;
                assert_eq!(*value, row_id, "Wrong row id for ({}, {})", row_id, row_offset);
            }
            if tag_key == &[b'r', b'o'] && *value == row_offset {
                assert!(!found_ro, "Duplicate row offset for ({}, {})", row_id, row_offset);
                found_ro = true;
                assert_eq!(*value, row_offset, "Wrong row offset for ({}, {})", row_id, row_offset);
            }
        }
    }

    aln.optional.retain(|tag| {
        if let TypedField::Int(tag_key, _) = tag {
            tag_key != &[b'r', b'i'] && tag_key != &[b'r', b'o']
        } else {
            true
        }
    });
}

#[test]
fn all_alignments() {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);

    let truth = formats::read_gaf_file(utils::get_test_data("micb-kir3dl1_HG003.gaf"));
    assert!(truth.is_ok(), "Failed to read GAF file: {}", truth.unwrap_err());
    let (_, truth) = truth.unwrap();

    for include_ids in [false, true] {
        let ids_str = if include_ids { "(with ids)" } else { "(without ids)" };
        let statement = GAFBaseInterface::all_alignments(&conn);
        assert!(statement.is_ok(), "Failed to prepare statement {}: {}", ids_str, statement.unwrap_err());
        let mut statement = statement.unwrap();
        let iter = GAFBaseIter::reference_free(&mut statement, (), &mut interface, None);
        let mut iter = match iter {
            Ok(iter) => iter,
            Err(err) => panic!("Failed to create iterator {}: {}", ids_str, err),
        };
        iter.include_ids(include_ids);
        let mut i = 0;
        for (row_id, block) in iter.enumerate() {
            let mut block = block.expect(&format!("Failed to get block starting at index {} {}", i, ids_str));
            // Because we are doing reference-free iteration, we must set target path length explicitly if we need it.
            for (row_offset, aln) in block.iter_mut().enumerate() {
                aln.set_target_path_len(truth[i].path_len);
                if include_ids {
                    check_ids(aln, row_id + 1, row_offset);
                }
                assert_eq!(*aln, truth[i], "Wrong alignment at index {} {}", i, ids_str);
                i += 1;
            }
        }
    }

    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

fn check_node_cache(cache: &HashMap<usize, GBZRecord>, alignments: &[Alignment], test_case: &str) {
    let mut expected_handles: HashSet<usize> = HashSet::new();
    for aln in alignments.iter() {
        if let Some(path) = aln.target_path() {
            for handle in path.iter() {
                expected_handles.insert(*handle);
            }
        }
    }

    assert_eq!(cache.len(), expected_handles.len(), "Node cache size mismatch for test case {}", test_case);
    for handle in expected_handles {
        let record = cache.get(&handle)
            .expect(&format!("Node {} missing in cache for test case {}", handle, test_case));
        assert!(record.sequence_len() > 0, "Node {} has empty sequence in cache for test case {}", handle, test_case);
    }
}

#[test]
fn row_ranges() {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);
    let gbz = internal::load_gaf_base_gbz();

    let truth = formats::read_gaf_file(utils::get_test_data("micb-kir3dl1_HG003.gaf"));
    assert!(truth.is_ok(), "Failed to read GAF file: {}", truth.unwrap_err());
    let (_, truth) = truth.unwrap();

    // Do full decompression of the GAF-base using a reference-based iterator processing two rows at a time.
    // This tests `GAFBaseIter::reference_based` and `GAFBaseInterface::alignment_rows`.
    // It additionally checks that the node record cache contains the right nodes.
    for include_ids in [true, false] {
        let ids_str = if include_ids { "(with ids)" } else { "(without ids)" };
        let mut get_sequence = |handle: usize| -> Result<Vec<u8>> {
            let (node, orientation) = support::decode_node(handle);
            let seq = gbz.sequence(node).ok_or_else(||
                Error::invalid_data(format!("The graph does not contain a sequence for node handle {} {}", handle, ids_str))
            )?;
            match orientation {
                Orientation::Forward => Ok(seq.to_vec()),
                Orientation::Reverse => Ok(support::reverse_complement(seq)),
            }
        };
        let mut statement = GAFBaseInterface::alignment_rows(&conn)
            .expect(&format!("Failed to prepare alignment rows statement {}", ids_str));

        let mut next_alignment = 0;
        let mut next_row: usize = 1; // Row ids start from 1.
        while next_alignment < truth.len() {
            let batch_start = next_alignment;
            let mut iter = GAFBaseIter::reference_based(
                &mut statement, (next_row, next_row + 2), &mut interface, &mut get_sequence, None
            ).expect(&format!("Failed to create iterator for rows {} to {} {}", next_row, next_row + 2, ids_str));
            iter.include_ids(include_ids);
            let mut num_rows = 0;
            let mut row_id = next_row;
            while let Some(block) = iter.next() {
                let mut block = block.expect(&format!("Failed to get block starting at alignment {} {}", next_alignment, ids_str));
                for (row_offset, aln) in block.iter_mut().enumerate() {
                    if include_ids {
                        check_ids(aln, row_id, row_offset);
                    }
                    assert_eq!(*aln, truth[next_alignment], "Wrong alignment at index {} {}", next_alignment, ids_str);
                    next_alignment += 1;
                }
                num_rows += 1;
                row_id += 1;
            }
            if next_alignment >= truth.len() {
                assert!(num_rows > 0, "Iterator did not process any rows in the final batch {}", ids_str);
            } else {
                assert_eq!(num_rows, 2, "Iterator did not process rows {}..{} {}", next_row, next_row + 2, ids_str);
            }
            let cache = iter.consume();
            check_node_cache(
                &cache, &truth[batch_start..next_alignment],
                &format!("rows {}..{} {}", next_row, next_row + 2, ids_str)
            );
            next_row += 2;
        }
    }

    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

fn boundary_node_pairs() -> Vec<(usize, usize)> {
    vec![
        (400, 403),   // Adjacent boundary nodes in the first component.
        (103, 255),   // Distant boundary nodes in the first component.
        (1981, 1984), // Adjacent boundary nodes in the second component.
        (1169, 1550), // Distant boundary nodes in the second component.
        (4000, 5000), // Non-existent nodes.
    ]
}

// We take the entire block if at least one alignment is in the interval.
fn collect_blocks_in_range(alignments: &[Alignment], min_handle: usize, max_handle: usize, include_ids: bool) -> Vec<Alignment> {
    let mut result = Vec::new();
    let mut start: usize = 0;
    let mut row_id: isize = 1;
    while start < alignments.len() {
        let mut end = start + 1;
        while end < alignments.len() && end - start < GAFBaseParams::BLOCK_SIZE && !alignments[end].is_unaligned() {
            end += 1;
        }
        let mut take_this_block = false;
        for i in start..end {
            if let Some(handles) = alignments[i].target_path() {
                for handle in handles.iter() {
                    if *handle >= min_handle && *handle <= max_handle {
                        take_this_block = true;
                        break;
                    }
                }
            }
            if take_this_block {
                break;
            }
        }
        if take_this_block {
            if include_ids {
                for i in start..end {
                    let mut aln = alignments[i].clone();
                    let row_offset = (i - start) as isize;
                    aln.optional.push(TypedField::Int([b'r', b'i'], row_id));
                    aln.optional.push(TypedField::Int([b'r', b'o'], row_offset));
                    result.push(aln);
                }
            } else {
                result.extend_from_slice(&alignments[start..end]);
            }
        }
        start = end;
        row_id += 1;
    }
    result
}

#[test]
fn handle_ranges() {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);
    let (_, alignments) = formats::read_gaf_file(utils::get_test_data("micb-kir3dl1_HG003.gaf"))
        .expect("Failed to read GAF file");

    for include_ids in [true, false] {
        let ids_str = if include_ids { "(with ids)" } else { "(without ids)" };
        let mut statement = GAFBaseInterface::handle_range(&conn)
            .expect(&format!("Failed to prepare handle range statement {}", ids_str));
        // We also test reusing the cache across multiple iterations.
        let mut cache = None;
        for (min_node, max_node) in boundary_node_pairs() {
            let min_handle = support::encode_node(min_node, Orientation::Forward);
            let max_handle = support::encode_node(max_node, Orientation::Reverse);
            let truth = collect_blocks_in_range(&alignments, min_handle, max_handle, include_ids);

            // Compare iterator output to the truth.
            let mut result: Vec<Alignment> = Vec::new();
            let mut iter = GAFBaseIter::reference_free(
                &mut statement, (min_handle, max_handle), &mut interface, cache
            ).expect(&format!("Failed to create iterator for node range {} to {} {}", min_node, max_node, ids_str));
            iter.include_ids(include_ids);
            while let Some(row) = iter.next() {
                let block = row.expect(&format!("Failed to read block for node range {} to {} {}", min_node, max_node, ids_str));
                result.extend(block);
            }
            cache = Some(iter.consume());
            assert_eq!(result.len(), truth.len(), "Wrong number of alignments for node range {} to {} {}", min_node, max_node, ids_str);
            for (i, (aln, truth_aln)) in result.iter_mut().zip(truth.iter()).enumerate() {
                // Reference-free iteration cannot determine target path lengths from a reference-based GAF-base.
                aln.set_target_path_len(truth_aln.path_len);
                assert_eq!(aln, truth_aln, "Mismatch at alignment {} for node range {} to {} {}", i, min_node, max_node, ids_str);
            }
        }
    }

    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

//-----------------------------------------------------------------------------
