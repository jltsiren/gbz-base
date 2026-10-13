use super::*;

use crate::{formats, internal};

use gbz::{GBZ, GBWT, Orientation, ENDMARKER};
use gbz::bwt::BWT;
use gbz::support;

use simple_sds::serialize;

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
        let mut cache = BTreeMap::new();
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

#[test]
fn all_alignments() {
    let gaf_base_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "");
    let conn = internal::open_gaf_base(&gaf_base_file);
    let mut interface = internal::create_gaf_base_interface(&conn);

    let truth = formats::read_gaf_file(utils::get_test_data("micb-kir3dl1_HG003.gaf"));
    assert!(truth.is_ok(), "Failed to read GAF file: {}", truth.unwrap_err());
    let (_, truth) = truth.unwrap();

    {
        let statement = GAFBaseInterface::all_alignments(&conn);
        assert!(statement.is_ok(), "Failed to prepare statement: {}", statement.unwrap_err());
        let mut statement = statement.unwrap();
        let iter = GAFBaseIter::reference_free(&mut statement, &mut interface);
        let iter = match iter {
            Ok(iter) => iter,
            Err(err) => panic!("Failed to create iterator: {}", err),
        };
        for (i, aln) in iter.enumerate() {
            assert!(aln.is_ok(), "Failed to get alignment at index {}: {}", i, aln.unwrap_err());
            let mut aln = aln.unwrap();
            // Because we are doing reference-free iteration, we must set target path length explicitly if we need it.
            aln.set_target_path_len(truth[i].path_len);
            assert_eq!(aln, truth[i], "Wrong alignment at index {}", i);
        }
    }

    drop(interface);
    drop(conn);
    let _ = fs::remove_file(&gaf_base_file);
}

//-----------------------------------------------------------------------------
