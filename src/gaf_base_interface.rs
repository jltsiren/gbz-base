//! GAF-base query interface.
// FIXME: document; connections to ReadSet

use crate::{GBZRecord, Alignment, AlignmentBlock};
use crate::{Result, Error};
use crate::alignment::{TargetPath, Flags};

use crate::utils;

use gbz::bwt::Record;

use rusqlite::{Connection, Statement, Row, Rows, OptionalExtension};

use std::collections::BTreeMap;
use std::iter::FusedIterator;

//#[cfg(test)]
//mod tests;

//-----------------------------------------------------------------------------

// FIXME: implement, document, test
#[derive(Debug)]
pub struct GAFBaseInterface<'a> {
    get_record: Statement<'a>,
}

impl<'a> GAFBaseInterface<'a> {
    /// Returns a new interface to the given database.
    ///
    /// Passes through any [`ErrorKind::Database`](crate::ErrorKind::Database) errors.
    pub fn new(connection: &'a Connection) -> Result<Self> {
        let get_record = connection.prepare("SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?")?;
        Ok(Self {
            get_record,
        })
    }

    // FIXME: inconvenient type for get_sequence
    /// Returns the node record for the given handle, or [`None`] if the node does not exist.
    ///
    /// If the GAF-base is reference-based, the `get_sequence` function will be used to retrieve the sequence.
    /// If the sequence is not necessary, this can be a trivial function.
    pub fn get_record(
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

    /// Decompresses an alignment block from a row, starting from index.
    ///
    /// Fields `from_idx` to `from_idx + 10` must correspond to columns `min_handle` to `optional` in table `Alignments`.
    pub fn decompress_block(row: &Row, from_idx: usize) -> Result<Vec<Alignment>> {
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

    // FIXME: document
    pub fn set_target_path(
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
                let record = self.get_record(p.node, get_sequence)?.ok_or(
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

// FIXME: document, tests
// Note: target path length is missing the right flank
pub struct GAFBaseIter<'a> {
    // Query results from table `Alignments` in GAF-base.
    // FIXME: which fields do we assume?
    rows: Rows<'a>,
    // GAF-base query interface for retrieving node records.
    interface: &'a mut GAFBaseInterface<'a>,
    // Cached GBZ records for the nodes used by the alignments so far.
    cache: BTreeMap<usize, GBZRecord>,
    // Alignments decoded from the current row.
    buffer: Vec<Alignment>,
    // Index of the next alignment to return from the buffer.
    buffer_index: usize,
    // Have we reached the end?
    finished: bool,
}

impl<'a> Iterator for GAFBaseIter<'a> {
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

impl<'a> FusedIterator for GAFBaseIter<'a> {}

impl<'a> GAFBaseIter<'a> {
    // FIXME: document
    pub fn new(rows: Rows<'a>, interface: &'a mut GAFBaseInterface<'a>) -> Result<Self> {
        Ok(Self {
            rows,
            interface,
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
            let mut buffer = GAFBaseInterface::decompress_block(&row, 0)?; // FIXME: which fields do we assume?
            for aln in &mut buffer {
                self.interface.set_target_path(aln, &mut self.cache, &mut |_| Ok(Vec::new()))?;
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
