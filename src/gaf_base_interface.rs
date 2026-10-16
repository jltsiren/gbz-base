//! GAF-base query interface.

use crate::{GAFBase, GBZRecord, Alignment, AlignmentBlock};
use crate::{Result, Error};
use crate::alignment::{TargetPath, Flags};

use crate::utils;

use gbz::bwt::Record;

use rusqlite::{Statement, Row, Rows, OptionalExtension, Params};

use std::collections::HashMap;
use std::iter::FusedIterator;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// [`GAFBase`] query interface.
///
/// This structure provides a low-level query interface to a GAF-base database.
/// It is mostly intended to be used with [`GAFBaseIter`].
/// See [`crate::ReadSet`] for a higher-level interface for extracting alignments to a subgraph.
///
/// # Examples
///
/// ```
/// use gbz_base::{GAFBase, GAFBaseParams, GraphReference};
/// use gbz_base::{GAFBaseInterface, GAFBaseIter};
/// use gbz_base::{Result, Error};
/// use gbz_base::{formats, utils};
/// use gbz::{GBZ, Orientation};
/// use gbz::support;
/// use simple_sds::serialize;
///
/// // Build a GAF-base.
/// let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
/// let db_file = serialize::temp_file_name("gaf-base");
/// GAFBase::create_from_files(
///     &gaf_file,
///     None, // Build a new GBWT index.
///     &db_file,
///     GraphReference::None, // Build a reference-based GAF-base.
///     &GAFBaseParams::default()
/// ).expect("Failed to create GAF-base");
/// let conn = GAFBase::open(&db_file)
///     .expect("Failed to open GAF-base");
/// let mut interface = GAFBaseInterface::new(&conn)
///     .expect("Failed to create GAF-base interface");
///
/// // We get sequences from an in-memory GBZ graph.
/// let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
/// let gbz: GBZ = serialize::load_from(&gbz_file)
///     .expect("Failed to load GBZ file");
/// let mut get_sequence = |handle: usize| -> Result<Vec<u8>> {
///     let (node, orientation) = support::decode_node(handle);
///     let seq = gbz.sequence(node).ok_or_else(||
///         Error::invalid_data(format!("Node {} does not exist", node))
///     )?;
///     if orientation == Orientation::Reverse {
///         Ok(support::reverse_complement(seq))
///     } else {
///         Ok(seq.to_vec())
///     }
/// };
///
/// // Iterate over the GAF-base and compare to the original file.
/// let (_headers, alignments) = formats::read_gaf_file(&gaf_file)
///     .expect("Failed to read GAF file");
/// {
///     let mut statement = GAFBaseInterface::all_alignments(&conn)
///         .expect("Failed to prepare all alignments statement");
///     let iter = GAFBaseIter::reference_free(&mut statement, (), &mut interface, None)
///         .expect("Failed to create GAF-base iterator");
///     let mut index = 0;
///     for block in iter {
///         let mut block = block.expect("Failed to retrieve next row");
///         for aln in block.iter_mut() {
///             // Reference-free iteration does not determine target path lengths.
///             aln.set_target_path_len(alignments[index].path_len);
///             assert_eq!(*aln, alignments[index], "Wrong alignment {}", index);
///             index += 1;
///         }
///     }
/// }
///
/// // Iterate over rows 4..6 that corresponds to alignments 3000..5000.
/// {
///     let mut statement = GAFBaseInterface::alignment_rows(&conn)
///         .expect("Failed to prepare alignment rows statement");
///     let mut iter = GAFBaseIter::reference_based(
///         &mut statement, (4, 6), &mut interface, &mut get_sequence, None
///     ).expect("Failed to create iterator for row 4");
///     let mut index = 3000;
///     while let Some(block) = iter.next() {
///         let mut block = block.expect("Failed to retrieve next row");
///         for aln in block.iter_mut() {
///             let aln = aln;
///             assert_eq!(*aln, alignments[index], "Wrong alignment {}", index);
///             index += 1;
///         }
///     }
/// }
///
/// drop(interface);
/// drop(conn);
/// let _ = std::fs::remove_file(&db_file);
/// ```
#[derive(Debug)]
pub struct GAFBaseInterface<'a> {
    get_record: Statement<'a>,
}

impl<'a> GAFBaseInterface<'a> {
    /// Returns a new interface to the given database.
    ///
    /// Passes through any [`ErrorKind::Database`](crate::ErrorKind::Database) errors.
    pub fn new(database: &'a GAFBase) -> Result<Self> {
        // We prefer not preparing other statements here, because we want to keep creating new interfaces cheap.
        let get_record = database.prepare("SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?")?;
        Ok(Self {
            get_record,
        })
    }

    // Returns the node record for the given handle.
    // If the GAF-base is reference-based, the `get_sequence` function will be used to retrieve the sequence.
    // The function is mutable to allow using a `GraphInterface`.
    // If the sequence is not necessary, this can be a trivial function.
    pub(crate) fn record(
        &mut self, handle: usize,
        get_sequence: &mut dyn FnMut(usize) -> Result<Vec<u8>>
    ) -> Result<Option<GBZRecord>> {
        let result = self.get_record.query_row(
            (handle,),
            |row| {
                let edge_bytes: Vec<u8> = row.get(0)?;
                let (edges, _) = Record::decompress_edges(&edge_bytes).unwrap();
                let bwt: Vec<u8> = row.get(1)?;
                let encoded_sequence: Vec<u8> = row.get(2)?;
                let sequence: Vec<u8> = utils::decode_sequence(&encoded_sequence);
                Ok((edges, bwt, sequence))
            }
        ).optional()?;
        if let Some((edges, bwt, mut sequence)) = result {
            if sequence.is_empty() {
                sequence = get_sequence(handle)?;
            }
            let record = unsafe {
                GBZRecord::from_raw_parts(handle, edges, bwt, sequence, None)
            };
            Ok(Some(record))
        } else {
            Ok(None)
        }
    }

    /// Returns a prepared statement for retrieving all alignments from table `Alignments`.
    ///
    /// The statement can be queried without any parameters.
    /// Alignments will be returned in the original order.
    /// This can be used for creating a [`GAFBaseIter`].
    pub fn all_alignments<'conn>(database: &'conn GAFBase) -> Result<Statement<'conn>> {
        database.prepare(
            "SELECT * FROM Alignments"
        )
    }

    /// Returns a prepared statement for retrieving the given range of rows from table `Alignments`.
    ///
    /// The statement takes two integer parameters that define a half-open interval of row ids.
    /// Alignments will be returned in the original order.
    /// This can be used for creating a [`GAFBaseIter`].
    pub fn alignment_rows<'conn>(database: &'conn GAFBase) -> Result<Statement<'conn>> {
        database.prepare(
            "SELECT * FROM Alignments WHERE id >= ? AND id < ?"
        )
    }

    /// Returns a prepared statement for retrieving rows overlapping the given range of handles from table `Alignments`.
    ///
    /// The statement takes two integer parameters: the minimum handle and the maximum handle.
    /// Alignment order may differ from the original if multiple rows share the same `min_handle` value.
    /// This can be used for creating a [`GAFBaseIter`].
    pub fn handle_range<'conn>(database: &'conn GAFBase) -> Result<Statement<'conn>> {
        database.prepare(
            "SELECT * FROM Alignments WHERE max_handle >= ? AND min_handle <= ?"
        )
    }

    // Decompresses an alignment block from a row, starting from index.
    // Fields `from_idx` to `from_idx + 10` must correspond to columns `min_handle` to `optional` in table `Alignments`.
    pub(crate) fn decompress_block(row: &Row, from_idx: usize) -> Result<Vec<Alignment>> {
        let min_handle: Option<usize> = row.get(from_idx + 0)?;
        let max_handle: Option<usize> = row.get(from_idx + 1)?;
        let alignments: usize = row.get(from_idx + 2)?;
        let read_length: Option<usize> = row.get(from_idx + 3)?;
        let gbwt_starts: Vec<u8> = row.get(from_idx + 4)?;
        let names: Vec<u8> = row.get(from_idx + 5)?;
        let quality_strings: Vec<u8> = row.get(from_idx + 6)?;
        let difference_strings: Vec<u8> = row.get(from_idx + 7)?;
        let flags: Vec<u8> = row.get(from_idx + 8)?;
        let numbers: Vec<u8> = row.get(from_idx + 9)?;
        let optional: Vec<u8> = row.get(from_idx + 10)?;
        let block = AlignmentBlock {
            min_handle, max_handle, alignments, read_length,
            gbwt_starts, names,
            quality_strings, difference_strings,
            flags: Flags::from(flags), numbers,
            optional,
        };
        block.decode()
    }

    // Sets the target path for an alignment.
    // See `Self::record` for assumptions on `get_sequence`.
    // If the sequences are not set, target path length will be 0.
    pub(crate) fn set_target_path(
        &mut self, aln: &mut Alignment,
        cache: &mut HashMap<usize, GBZRecord>,
        get_sequence: &mut dyn FnMut(usize) -> Result<Vec<u8>>
    ) -> Result<()> {
        let mut pos = match &aln.path {
            TargetPath::StartPosition(p) => Some(*p),
            _ => return Ok(()),
        };
        let mut path = Vec::new();
        let mut path_len = 0;
        while let Some(p) = pos {
            path.push(p.node);
            if !cache.contains_key(&p.node) {
                let record = self.record(p.node, get_sequence)?.ok_or(
                    Error::invalid_data(format!("GAF-base does not contain a record for node handle {}", p.node))
                )?;
                cache.insert(p.node, record);
            }
            let record = cache.get(&p.node).unwrap();
            path_len += record.sequence_len();
            pos = record.to_gbwt_record().lf(p.offset);
        }
        aln.set_target_path(path);
        aln.set_target_path_len(path_len);
        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// An iterator over alignments in a GAF-base.
///
/// This iterator can do either reference-free or reference-based iteration.
/// Reference-free iteration uses node records directly from the GAF-base.
/// If the GAF-base is reference-based, target path lengths will not be set in the returned alignments.
/// They can be set with [`Alignment::set_target_path_len`] if needed.
/// Reference-based iteration avoids this by using another source for node sequences.
pub struct GAFBaseIter<'reference, 'interface> {
    // Query results from table `Alignments`.
    rows: Rows<'reference>,
    // GAF-base query interface for retrieving node records.
    record_interface: &'reference mut GAFBaseInterface<'interface>,
    // Optional function for getting node sequences from another source.
    get_sequence: Option<&'reference mut dyn FnMut(usize) -> Result<Vec<u8>>>,
    // Cached GBZ records for the nodes used by the alignments so far.
    cache: HashMap<usize, GBZRecord>,
}

impl<'reference, 'interface> Iterator for GAFBaseIter<'reference, 'interface> {
    type Item = Result<Vec<Alignment>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_next_row() {
            Ok(block) => {
                if block.is_empty() {
                    None
                } else {
                    Some(Ok(block))
                }
            },
            Err(e) => Some(Err(e))
        }
    }
}

impl<'reference, 'interface> FusedIterator for GAFBaseIter<'reference, 'interface> {}

impl<'reference, 'interface> GAFBaseIter<'reference, 'interface> {
    // Column index of `min_handle` in query results with all columns included.
    const FROM_IDX: usize = 1;

    /// Returns a reference-free iterator over the alignments returned by the given statement.
    ///
    /// As the iterator does reference-free iteration, target path lengths are not set in the alignments.
    /// If needed, they can be set with [`Alignment::set_target_path_len`].
    /// See [`GAFBaseInterface`] for an example.
    ///
    /// # Arguments
    ///
    /// * `statement`: Prepared SQL statement that returns rows with all columns from table `Alignments`.
    /// * `params`: Parameters to bind to the SQL statement.
    /// * `record_interface`: The GAF-base interface for retrieving node records.
    /// * `cache`: Optional node record cache that is expected to overlap with the iterated alignments.
    pub fn reference_free<P: Params>(
        statement: &'reference mut Statement<'_>, params: P,
        record_interface: &'reference mut GAFBaseInterface<'interface>,
        cache: Option<HashMap<usize, GBZRecord>>,
    ) -> Result<Self> {
        let rows = statement.query(params)?;
        Ok(Self {
            rows,
            record_interface,
            get_sequence: None,
            cache: cache.unwrap_or_default(),
        })
    }

    /// Returns a reference-based iterator over the alignments returned by the given statement.
    ///
    /// If the GAF-base is reference-based, this version uses the given function to determine node sequences.
    /// After iteration has finished, node records with sequences can be obtained with [`Self::consume`].
    /// See [`GAFBaseInterface`] for an example.
    ///
    /// # Arguments
    ///
    /// * `statement`: Prepared SQL statement that returns rows with all columns from table `Alignments`.
    /// * `record_interface`: The GAF-base interface for retrieving node records.
    /// * `get_sequence`: Function for retrieving node sequences.
    /// * `cache`: Optional node record cache that is expected to overlap with the iterated alignments.
    pub fn reference_based<P: Params>(
        statement: &'reference mut Statement<'_>, params: P,
        record_interface: &'reference mut GAFBaseInterface<'interface>,
        get_sequence: &'reference mut dyn FnMut(usize) -> Result<Vec<u8>>,
        cache: Option<HashMap<usize, GBZRecord>>,
    ) -> Result<Self> {
        let rows = statement.query(params)?;
        Ok(Self {
            rows,
            record_interface,
            get_sequence: Some(get_sequence),
            cache: cache.unwrap_or_default(),
        })
    }

    /// Returns the number of node records in the cache.
    pub fn cache_size(&self) -> usize {
        self.cache.len()
    }

    /// Consumes the iterator and returns the node record cache.
    ///
    /// This can be useful for GAF output after reference-based iteration.
    /// Some difference string operations require target bases, which can be determined from the cache.
    pub fn consume(self) -> HashMap<usize, GBZRecord> {
        self.cache
    }

    // Returns the next row as a block of alignments, or an empty block if there are no more rows.
    fn read_next_row(&mut self) -> Result<Vec<Alignment>> {
        let result = self.rows.next().map_err(Error::from)?;
        if let Some(row) = result {
            let mut block = GAFBaseInterface::decompress_block(row, Self::FROM_IDX)?;
            for aln in block.iter_mut() {
                if let Some(get_sequence) = &mut self.get_sequence {
                    self.record_interface.set_target_path(aln, &mut self.cache, get_sequence)?;
                } else {
                    self.record_interface.set_target_path(aln, &mut self.cache, &mut |_| Ok(Vec::new()))?;
                }
            }
            Ok(block)
        } else {
            Ok(Vec::new())
        }
    }
}

//-----------------------------------------------------------------------------
