//! GAF-base query interface.

use crate::{GAFBase, GBZRecord, Alignment, AlignmentBlock};
use crate::{Result, Error};
use crate::alignment::{TargetPath, Flags};

use crate::utils;

use gbz::bwt::Record;

use rusqlite::{Statement, Row, Rows, OptionalExtension};

use std::collections::BTreeMap;
use std::iter::FusedIterator;

//#[cfg(test)]
//mod tests;

//-----------------------------------------------------------------------------

// FIXME: test
/// [`GAFBase`] query interface.
///
/// This structure provides a low-level query interface to a GAF-base database.
/// See [`crate::ReadSet`] for a higher-level interface for extracting alignments to a subgraph.
///
/// # Examples
///
/// ```
/// use gbz_base::{GAFBase, GAFBaseParams, GraphReference};
/// use gbz_base::{GAFBaseInterface, GAFBaseIter};
/// use gbz_base::{Alignment, Result, Error};
/// use gbz_base::{formats, utils};
/// use gbz::{GBZ, GBWT, Orientation};
/// use gbz::bwt::BWT;
/// use gbz::support;
/// use simple_sds::serialize;
///
/// // Build a GAF-base.
/// let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
/// let db_file = serialize::temp_file_name("gaf-base");
/// let db = GAFBase::create_from_files(
///     &gaf_file,
///     None, // Build a new GBWT index.
///     &db_file,
///     GraphReference::None, // Build a reference-based GAF-base.
///     &GAFBaseParams::default()
/// );
/// assert!(db.is_ok(), "Failed to build GAF-base");
/// let conn = GAFBase::open(&db_file);
/// assert!(conn.is_ok(), "Failed to open GAF-base");
/// let conn = conn.unwrap();
/// let interface = GAFBaseInterface::new(&conn);
/// assert!(interface.is_ok(), "Failed to create GAF-base interface");
/// let mut interface = interface.unwrap();
///
/// // We get sequences from an in-memory GBZ graph.
/// let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
/// let gbz: GBZ = serialize::load_from(&gbz_file).unwrap();
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
/// // GBZ records combine a GBWT record and a sequence.
/// let gbwt_file = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
/// let gbwt: GBWT = serialize::load_from(&gbwt_file).unwrap();
/// let (node, handle) = (1, support::encode_node(1, Orientation::Forward));
/// let gbz_record = interface.record(handle, &mut get_sequence);
/// assert!(gbz_record.is_ok(), "Failed to query GBZ record");
/// let gbz_record = gbz_record.unwrap();
/// assert!(gbz_record.is_some(), "GBZ record does not exist");
/// let gbz_record = gbz_record.unwrap();
/// assert_eq!(gbz_record.handle(), handle);
/// let sequence = gbz_record.sequence();
/// assert_eq!(Some(sequence), gbz.sequence(node));
/// let bwt: &BWT = gbwt.as_ref();
/// let record_id = gbwt.node_to_record(handle);
/// let from_gbwt = bwt.record(record_id);
/// let from_db = gbz_record.to_gbwt_record_with_id(record_id);
/// assert_eq!(Some(from_db), from_gbwt);
///
/// // Iterate over the GAF-base using a second connection.
/// let truth = formats::parse_gaf_file(&gaf_file);
/// assert!(
///     truth.is_ok(), "Failed to parse GAF file {}: {}",
///     gaf_file.display(), truth.unwrap_err()
/// );
/// let (_headers, truth) = truth.unwrap();
/// let iter_conn = GAFBase::open(&db_file);
/// assert!(iter_conn.is_ok(), "Failed to open GAF-base for iteration");
/// let iter_conn = iter_conn.unwrap();
/// let iter_interface = GAFBaseInterface::new(&iter_conn);
/// assert!(iter_interface.is_ok(), "Failed to create GAF-base interface for iteration");
/// let mut iter_interface = iter_interface.unwrap();
/// let iter = GAFBaseIter::all_alignments(&mut iter_interface, &mut interface);
/// assert!(iter.is_ok(), "Failed to create GAF-base iterator");
/// let iter = iter.unwrap();
/// for (i, aln) in iter.enumerate() {
///     assert!(aln.is_ok(), "Failed to retrieve alignment");
///     let mut aln = aln.unwrap();
///     aln.set_target_path_len(truth[i].path_len);
///     assert_eq!(aln, truth[i], "Wrong alignment {}", i);
/// }
///
/// drop(iter_interface);
/// drop(iter_conn);
/// drop(interface);
/// drop(conn);
/// let _ = std::fs::remove_file(&db_file);
/// ```
#[derive(Debug)]
pub struct GAFBaseInterface<'a> {
    get_record: Statement<'a>,
    get_all_alignments: Statement<'a>,
}

impl<'a> GAFBaseInterface<'a> {
    /// Returns a new interface to the given database.
    ///
    /// Passes through any [`ErrorKind::Database`](crate::ErrorKind::Database) errors.
    pub fn new(database: &'a GAFBase) -> Result<Self> {
        let get_record = database.connection.prepare("SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?")?;
        let get_all_alignments = database.connection.prepare(
            "SELECT * FROM Alignments ORDER BY id"
        )?;
        Ok(Self {
            get_record,
            get_all_alignments,
        })
    }

    /// Returns the node record for the given handle, or [`None`] if the node does not exist.
    ///
    /// If the GAF-base is reference-based, the `get_sequence` function will be used to retrieve the sequence.
    /// The function is mutable to allow using a [`crate::GraphInterface`].
    /// If the sequence is not necessary, this can be a trivial function.
    pub fn record(
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
        cache: &mut BTreeMap<usize, GBZRecord>,
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

// FIXME: tests
// Note: target path length is missing the right flank
/// An iterator over alignments in a GAF-base.
///
/// The iterator uses two database connetions: one for the alignments and another for node records.
/// Alignments will be iterated in the order they appear in the database.
pub struct GAFBaseIter<'reference, 'conn> {
    // Query results from table `Alignments` in GAF-base.
    // The results must include all fields from `min_handle` to `optional` contiguously.
    rows: Rows<'reference>,
    // Field number corresponding to `min_handle`.
    from_idx: usize,
    // GAF-base query interface for retrieving node records.
    record_interface: &'reference mut GAFBaseInterface<'conn>,
    // Cached GBZ records for the nodes used by the alignments so far.
    cache: BTreeMap<usize, GBZRecord>,
    // Alignments decoded from the current row.
    buffer: Vec<Alignment>,
    // Index of the next alignment to return from the buffer.
    buffer_index: usize,
    // Have we reached the end?
    finished: bool,
}

impl<'reference, 'conn> Iterator for GAFBaseIter<'reference, 'conn> {
    type Item = Result<Alignment>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        if self.buffer_index >= self.buffer.len() {
            let result = self.read_next_row();
            if let Err(e) = result {
                return Some(Err(e));
            }
            if self.finished {
                return None;
            }
        }
        let alignment = self.buffer[self.buffer_index].clone();
        self.buffer_index += 1;
        Some(Ok(alignment))
    }
}

impl<'reference, 'conn> FusedIterator for GAFBaseIter<'reference, 'conn> {}

impl<'reference, 'conn> GAFBaseIter<'reference, 'conn> {
    /// Returns an iterator over all alignments in the GAF-base.
    ///
    /// Because this does reference-free iteration, target path lengths are unset in the alignments.
    /// If needed, they can be set with [`Alignment::set_target_path_len`].
    /// See [`GAFBaseInterface`] for an example.
    pub fn all_alignments(iter_interface: &'reference mut GAFBaseInterface<'_>, record_interface: &'reference mut GAFBaseInterface<'conn>) -> Result<Self> {
        let rows = iter_interface.get_all_alignments.query(())?;
        Ok(Self {
            rows,
            from_idx: 1,
            record_interface,
            cache: BTreeMap::new(),
            buffer: Vec::new(),
            buffer_index: 0,
            finished: false,
        })
    }

    // Reads the next row into the buffer.
    // Passes through database errors.
    fn read_next_row(&mut self) -> Result<()> {
        let result = self.rows.next().map_err(Error::from)?;
        if let Some(row) = result {
            let mut buffer = GAFBaseInterface::decompress_block(&row, self.from_idx)?;
            for aln in &mut buffer {
                self.record_interface.set_target_path(aln, &mut self.cache, &mut |_| Ok(Vec::new()))?;
            }
            self.buffer = buffer;
            self.buffer_index = 0;
        } else {
            self.finished = true;
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------
