//! A set of reads extracted from a GAF-base.

use crate::{GAFBase, GBZRecord, GraphReference, Subgraph, Alignment, AlignmentBlock};
use crate::alignment::{Flags, TargetPath};

use gbwt::Pos;
use gbwt::bwt::Record;

use rusqlite::{Row, OptionalExtension};

use std::collections::BTreeMap;
use std::fmt::Display;
use std::io::Write;
use std::sync::Arc;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Output options for alignments in a subgraph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AlignmentOutput {
    /// Reads overlapping with the subgraph.
    Overlapping,
    /// Overlapping reads clipped to the subgraph.
    Clipped,
    /// Reads fully contained within the subgraph.
    Contained,
}

impl Display for AlignmentOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            AlignmentOutput::Overlapping => write!(f, "overlapping"),
            AlignmentOutput::Clipped => write!(f, "clipped"),
            AlignmentOutput::Contained => write!(f, "contained"),
        }
    }
}

//-----------------------------------------------------------------------------

/// A set of reads extracted from [`GAFBase`].
///
/// This is a counterpart to [`Subgraph`].
/// Sets of reads fully contained in a subgraph or overlapping with it can be created using [`ReadSet::new`].
/// The reads can be iterated over with [`ReadSet::iter`] and converted to GAF lines with [`ReadSet::to_gaf`].
/// The reads will appear in the same order as in the database.
///
/// # Examples
///
/// ```
/// use gbz_base::{Subgraph, SubgraphQuery, HaplotypeOutput};
/// use gbz_base::{GAFBase, GAFBaseParams, ReadSet, GraphReference};
/// use gbz_base::utils;
/// use gbwt::GBZ;
/// use simple_sds::serialize;
///
/// // Get an in-memory graph.
/// let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
/// let graph = serialize::load_from(&gbz_file).unwrap();
///
/// // Extract a 100 bp subgraph around node 150.
/// let nodes = vec![150];
/// let query = SubgraphQuery::nodes(nodes).with_output(HaplotypeOutput::Distinct);
/// let mut subgraph = Subgraph::new();
/// let _ = subgraph.from_gbz(&graph, None, None, &query).unwrap();
///
/// // Create a database of reads aligned to the graph.
/// let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf");
/// let gbwt_file = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
/// let db_file = serialize::temp_file_name("gaf-base");
/// let params = GAFBaseParams::default();
/// let db = GAFBase::create_from_files(&gaf_file, &gbwt_file, &db_file, &params);
/// assert!(db.is_ok());
///
/// // Extract all reads fully within the subgraph.
/// let db = GAFBase::open(&db_file);
/// assert!(db.is_ok());
/// let db = db.unwrap();
/// let read_set = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &db, true);
/// assert!(read_set.is_ok());
/// let read_set = read_set.unwrap();
/// assert_eq!(read_set.len(), 148);
///
/// // The extracted reads are aligned and fully within the subgraph.
/// for aln in read_set.iter() {
///     for handle in aln.target_path().unwrap() {
///         assert!(subgraph.has_handle(*handle));
///     }
/// }
///
/// drop(db);
/// let _ = std::fs::remove_file(&db_file);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ReadSet {
    // GBZ records from the GAF GBWT, with sequence from the graph/subgraph.
    nodes: BTreeMap<usize, GBZRecord>,
    reads: Vec<Alignment>,
    // Number of alignments before clipping.
    unclipped: usize,
    blocks: usize,
}

impl ReadSet {
    // Decompresses an alignment block from a row.
    fn decompress_block(row: &Row) -> Result<Vec<Alignment>, String> {
        let min_handle: Option<usize> = row.get(0).map_err(|x| x.to_string())?;
        let max_handle: Option<usize> = row.get(1).map_err(|x| x.to_string())?;
        let alignments: usize = row.get(2).map_err(|x| x.to_string())?;
        let read_length: Option<usize> = row.get(3).map_err(|x| x.to_string())?;
        let gbwt_starts: Vec<u8> = row.get(4).map_err(|x| x.to_string())?;
        let names: Vec<u8> = row.get(5).map_err(|x| x.to_string())?;
        let quality_strings: Vec<u8> = row.get(6).map_err(|x| x.to_string())?;
        let difference_strings: Vec<u8> = row.get(7).map_err(|x| x.to_string())?;
        let flags: Vec<u8> = row.get(8).map_err(|x| x.to_string())?;
        let numbers: Vec<u8> = row.get(9).map_err(|x| x.to_string())?;

        let block = AlignmentBlock {
            min_handle, max_handle, alignments, read_length,
            gbwt_starts, names,
            quality_strings, difference_strings,
            flags: Flags::from(flags), numbers
        };
        block.decode()
    }

    // Replaces the GBWT starting position of the alignment with the path.
    // Requires that the path overlaps with / is fully contained in the subgraph.
    // If the path is valid, inserts all missing node records into the read set.
    fn set_target_path(
        &mut self, alignment: &mut Alignment, subgraph: &Subgraph,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>,
        contained: bool
    ) -> Result<(), String> {
        let mut pos = match alignment.path {
            TargetPath::Path(_) => return Ok(()),
            TargetPath::StartPosition(pos) => Some(pos),
        };

        let mut path = Vec::new();
        let mut overlap = false;
        while let Some(p) = pos {
            if subgraph.has_handle(p.node) {
                overlap = true;
            } else if contained {
                return Ok(()); // Not fully contained in the subgraph.
            }

            // Now get the record for the node.
            let mut record = self.nodes.get(&p.node);
            if record.is_none() {
                let result = get_record(p.node)?;
                self.nodes.insert(p.node, result);
                record = self.nodes.get(&p.node);
            }

            // Navigate to the next position.
            path.push(p.node);
            let record = record.unwrap();
            pos = record.to_gbwt_record().lf(p.offset);
        }

        // Set the target path in the alignment.
        if overlap {
            alignment.set_target_path(path);
        }
        Ok(())
    }

    /// Extracts a set of reads overlapping with the subgraph.
    ///
    /// The extracted reads will be in the same order as in the database.
    /// That corresponds to the order in the original GAF file.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ-compatible graph.
    /// * `subgraph`: The subgraph used as the query region.
    /// * `database`: A database storing reads aligned to the graph.
    /// * `output`: Which reads to include in the read set.
    ///
    /// # Errors
    ///
    /// Passes through any database errors.
    /// Returns an error if an alignment cannot be decompressed.
    pub fn new(graph: GraphReference<'_, '_>, subgraph: &Subgraph, database: &GAFBase, output: AlignmentOutput) -> Result<Self, String> {
        let mut read_set = ReadSet {
            nodes: BTreeMap::new(),
            reads: Vec::new(),
            unclipped: 0,
            blocks: 0,
        };

        // Build a record from the databases.
        let mut get_node = database.connection.prepare(
            "SELECT edges, bwt FROM Nodes WHERE handle = ?1"
        ).map_err(|x| x.to_string())?;
        let mut graph = graph;
        let mut get_record = |handle: usize| -> Result<GBZRecord, String> {
            // Get the edges and the BWT fragment from the GAF-base.
            let gaf_result = get_node.query_row(
                (handle,),
                |row: &Row<'_>| -> rusqlite::Result<(Vec<Pos>, Vec<u8>)> {
                    let edge_bytes: Vec<u8> = row.get(0)?;
                    let (edges, _) = Record::decompress_edges(&edge_bytes).unwrap();
                    let bwt: Vec<u8> = row.get(1)?;
                    Ok((edges, bwt))
                }
            ).optional().map_err(|x| x.to_string())?;
            if gaf_result.is_none() {
                return Err(format!("Could not find the record for handle {} in GAF-base", handle));
            }
            let (edges, bwt) = gaf_result.unwrap();

            // Get the sequence from the subgraph or from the GBZ-base.
            let sequence = subgraph.sequence(handle);
            let sequence = match sequence {
                Some(seq) => seq.to_vec(),
                None => {
                    let gbz_record = graph.gbz_record(handle)?;
                    gbz_record.sequence().to_vec()
                }
            };

            unsafe {
                Ok(GBZRecord::from_raw_parts(handle, edges, bwt, sequence, None))
            }
        };

        // Get the reads that may overlap with the subgraph.
        let min_handle = subgraph.min_handle();
        let max_handle = subgraph.max_handle();
        let mut get_reads = database.connection.prepare(
            "SELECT min_handle, max_handle, alignments, read_length, gbwt_starts, names, quality_strings, difference_strings, flags, numbers
            FROM Alignments
            WHERE min_handle <= ?1 AND max_handle >= ?2"
        ).map_err(|x| x.to_string())?;
        let mut rows = get_reads.query((max_handle, min_handle)).map_err(|x| x.to_string())?;
        while let Some(row) = rows.next().map_err(|x| x.to_string())? {
            let block = Self::decompress_block(row)?;
            for mut alignment in block {
                read_set.set_target_path(&mut alignment, subgraph, &mut get_record, output == AlignmentOutput::Contained)?;
                if alignment.has_target_path() {
                    if output == AlignmentOutput::Clipped {
                        let sequence_len = Arc::new(|handle| {
                            let record = read_set.nodes.get(&handle)?;
                            Some(record.sequence().len())
                        });
                        let clipped = alignment.clip(subgraph, sequence_len)?;
                        for aln in clipped.into_iter() {
                            read_set.reads.push(aln);
                        }
                    } else {
                        read_set.reads.push(alignment);
                    }
                    read_set.unclipped += 1;
                }
            }
            read_set.blocks += 1;
        }

        Ok(read_set)
    }

    /// Returns the number of alignment fragments in the set.
    #[inline]
    pub fn len(&self) -> usize {
        self.reads.len()
    }

    /// Returns `true` if the set is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    /// Returns the original number of alignments (before clipping) in the set.
    #[inline]
    pub fn unclipped(&self) -> usize {
        self.unclipped
    }

    /// Returns the number of alignment blocks decompressed when creating the read set.
    #[inline]
    pub fn blocks(&self) -> usize {
        self.blocks
    }

    /// Returns the number of node records in the read set.
    ///
    /// Each record corresponds to an oriented node, and the opposite orientation may not be present.
    /// This includes all node records encountered while tracing the alignments, even when the alignment was not included in the read set.
    #[inline]
    pub fn node_records(&self) -> usize {
        self.nodes.len()
    }

    /// Returns an iterator over the reads in the set.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &Alignment> {
        self.reads.iter()
    }

    // Extracts the target sequence for the given alignment.
    fn target_sequence(&self, alignment: &Alignment) -> Result<Vec<u8>, String> {
        let target_path = alignment.target_path();
        if target_path.is_none() {
            return Ok(Vec::new());
        }
        let target_path = target_path.unwrap();

        let mut sequence = Vec::new();
        for handle in target_path{
            let record = self.nodes.get(handle);
            if record.is_none() {
                return Err(format!("Read {}: Missing record for node handle {}", alignment.name, handle));
            }
            let record = record.unwrap();
            sequence.extend_from_slice(record.sequence());
        }

        if sequence.len() != alignment.path_len {
            return Err(format!(
                "Read {}: Target path length {} does not match the expected length {}",
                alignment.name, sequence.len(), alignment.path_len
            ));
        }
        sequence = sequence[alignment.path_interval.clone()].to_vec();
        Ok(sequence)
    }

    /// Serializes the read set in the GAF format.
    ///
    /// Returns an error if the target sequence for a read is invalid or cannot be determined.
    /// Passes through any I/O errors.
    pub fn to_gaf<W: Write>(&self, writer: &mut W) -> Result<(), String> {
        for alignment in self.reads.iter() {
            let target_sequence = self.target_sequence(alignment)?;
            let mut line = alignment.to_gaf(&target_sequence);
            line.push(b'\n');
            writer.write_all(&line).map_err(|x| x.to_string())?;
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------
