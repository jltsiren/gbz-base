//! A set of reads extracted from a GAF-base.

use crate::{GAFBase, GBZRecord, GraphReference, Subgraph, Alignment, GAFBaseInterface};
use crate::alignment::TargetPath;
use crate::error::{Error, Result};
use crate::utils;

use gbz::{Orientation, Pos, GBZ};
use gbz::bwt::Record;
use gbz::support;

use rusqlite::{Row, OptionalExtension};

use std::collections::{BTreeMap, HashSet};
use std::fmt::Display;
use std::io::Write;
use std::ops::Range;
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
/// use gbz_base::{GAFBase, GAFBaseParams, ReadSet, AlignmentOutput, GraphReference};
/// use gbz_base::utils;
/// use gbz::GBZ;
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
/// let gbwt_file = None; // Build a new GBWT index.
/// let db_file = serialize::temp_file_name("gaf-base");
/// let graph_ref = GraphReference::None; // Do not store sequences in the database.
/// let params = GAFBaseParams::default();
/// let db = GAFBase::create_from_files(&gaf_file, gbwt_file, &db_file, graph_ref, &params);
/// assert!(db.is_ok());
///
/// // Extract all reads fully within the subgraph.
/// let db = GAFBase::open(&db_file);
/// assert!(db.is_ok());
/// let db = db.unwrap();
/// let read_set = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &db, AlignmentOutput::Contained);
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
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReadSet {
    // GBZ records from the GAF GBWT, with sequence from the graph/subgraph.
    nodes: BTreeMap<usize, GBZRecord>,
    reads: Vec<Alignment>,
    // Number of alignments before clipping.
    unclipped: usize,
    // Number of alignment blocks decompressed.
    blocks: usize,
    // Number of candidate alignments in the decompressed blocks.
    candidates: usize,
    // Number of node id clusters in the subgraph.
    clusters: usize,
}

impl ReadSet {
    // TODO: Should this be configurable?
    /// Gap length threshold for clustering node ids.
    pub const CLUSTER_GAP_THRESHOLD: usize = 1000;

    // Returns the row id from the first column.
    fn get_row_id(row: &Row) -> Result<usize> {
        let row_id: usize = row.get(0)?;
        Ok(row_id)
    }

    // Replaces the GBWT starting position of the alignment with the path and sets the true target path length.
    // Requires that the path overlaps with / is fully contained in the subgraph.
    // If the path is valid, inserts all missing node records into the read set.
    fn set_target_path(
        &mut self, alignment: &mut Alignment, subgraph: &Subgraph,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord>,
        contained: bool
    ) -> Result<()> {
        let mut pos = match alignment.path {
            TargetPath::Path(_) => return Ok(()),
            TargetPath::StartPosition(pos) => Some(pos),
        };

        let mut path = Vec::new();
        let mut overlap = false;
        let mut target_path_len = 0;
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
            target_path_len += record.sequence_len();
            pos = record.to_gbwt_record().lf(p.offset);
        }

        // Set the target path in the alignment.
        if overlap {
            alignment.set_target_path(path);
            alignment.set_target_path_len(target_path_len);
        }
        Ok(())
    }

    // Replaces the GBWT starting position of the alignment with the path and sets the true target path length.
    // Inserts all missing node records into the read set.
    fn set_target_path_simple(
        &mut self, alignment: &mut Alignment,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord>,
    ) -> Result<()> {
        let mut pos = match alignment.path {
            TargetPath::Path(_) => return Ok(()),
            TargetPath::StartPosition(pos) => Some(pos),
        };

        let mut path = Vec::new();
        let mut target_path_len = 0;
        while let Some(p) = pos {
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
            target_path_len += record.sequence_len();
            pos = record.to_gbwt_record().lf(p.offset);
        }

        // Set the target path in the alignment.
        alignment.set_target_path(path);
        alignment.set_target_path_len(target_path_len);
        Ok(())
    }

    // FIXME: use GAFBaseInterface
    /// Extracts a set of reads overlapping with the subgraph.
    ///
    /// The extracted reads will be in the same order as in the database.
    /// That corresponds to the order in the original GAF file.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ-compatible graph for querying a reference-based GAF-base, or no graph for a reference-free one.
    /// * `subgraph`: The subgraph used as the query region.
    /// * `database`: A database storing reads aligned to the graph.
    /// * `output`: Which reads to include in the read set.
    ///
    /// # Errors
    ///
    /// Returns an [`ErrorKind::InvalidData`](crate::ErrorKind::InvalidData) error if the database is missing a node record that an alignment refers to.
    /// Passes through any [`ErrorKind::Database`](crate::ErrorKind::Database) errors.
    /// Returns an error if an alignment cannot be decompressed.
    pub fn new(graph: GraphReference<'_, '_>, subgraph: &Subgraph, database: &GAFBase, output: AlignmentOutput) -> Result<Self> {
        let mut read_set = ReadSet::default();

        // Build a record from the databases.
        let mut get_node = database.connection.prepare(
            "SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?1"
        )?;
        let mut graph = graph;
        let mut get_record = |handle: usize| -> Result<GBZRecord> {
            // Get the node record from the GAF-base.
            let gaf_result = get_node.query_row(
                (handle,),
                |row: &Row<'_>| -> rusqlite::Result<(Vec<Pos>, Vec<u8>, Vec<u8>)> {
                    let edge_bytes: Vec<u8> = row.get(0)?;
                    let (edges, _) = Record::decompress_edges(&edge_bytes).unwrap();
                    let bwt: Vec<u8> = row.get(1)?;
                    let encoded_sequence: Vec<u8> = row.get(2)?;
                    let sequence = utils::decode_sequence(&encoded_sequence);
                    Ok((edges, bwt, sequence))
                }
            ).optional()?;
            if gaf_result.is_none() {
                return Err(Error::invalid_data(format!("GAF-base does not contain a record for node handle {}", handle)));
            }
            let (edges, bwt, mut sequence) = gaf_result.unwrap();

            // If we have a reference-based database, try to get the sequence from the
            // subgraph or the original graph.
            if sequence.is_empty() {
                let seq = subgraph.sequence_for_handle(handle);
                sequence = match seq {
                    Some(seq) => seq.to_vec(),
                    None => {
                        let gbz_record = graph.gbz_record(handle)?;
                        gbz_record.sequence().to_vec()
                    }
                };
            }

            unsafe {
                Ok(GBZRecord::from_raw_parts(handle, edges, bwt, sequence, None))
            }
        };

        // Cluster the handles in the subgraph into reasonable intervals. Because node
        // i corresponds to handles 2i and 2i+1, it is easier to work with node ids.
        let node_ids: Vec<usize> = subgraph.node_iter().collect();
        let clusters = utils::cluster_node_ids(node_ids, Self::CLUSTER_GAP_THRESHOLD);
        let clusters: Vec<(usize, usize)> = clusters.into_iter()
            .map(|r| (support::encode_node(*r.start(), Orientation::Forward), support::encode_node(*r.end(), Orientation::Reverse)))
            .collect();
        read_set.clusters = clusters.len();

        // Get the reads that may overlap with the subgraph. We keep track of row ids
        // we have encountered to avoid duplicating reads that overlap multiple clusters.
        let mut row_ids: HashSet<usize> = HashSet::new();
        let mut get_reads = database.connection.prepare(
            "SELECT id, min_handle, max_handle, alignments, read_length, gbwt_starts, names, quality_strings, difference_strings, flags, numbers, optional
            FROM Alignments
            WHERE min_handle <= ?1 AND max_handle >= ?2"
        )?;
        for (min_handle, max_handle) in clusters.into_iter() {
            let mut rows = get_reads.query((max_handle, min_handle))?;
            while let Some(row) = rows.next()? {
                let row_id = Self::get_row_id(row)?;
                if row_ids.contains(&row_id) {
                    continue;
                }
                row_ids.insert(row_id);
                let block = GAFBaseInterface::decompress_block(row, 1)?;
                let block_size = block.len();
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
                read_set.candidates += block_size;
            }
        }

        Ok(read_set)
    }

    // FIXME: use GAFBaseInterface
    /// Extracts all reads from the given range of row ids.
    ///
    /// The extracted reads will be in the same order as in the database.
    /// That corresponds to the order in the original GAF file.
    ///
    /// # Arguments
    ///
    /// * `database`: A database storing reads aligned to the graph.
    /// * `row_range`: The range of row ids to extract.
    /// * `graph`: A GBZ graph if the database is reference-based, or [`None`] for a reference-free one.
    ///
    /// # Errors
    ///
    /// Returns an [`ErrorKind::InvalidData`](crate::ErrorKind::InvalidData) error if the database is missing a node record that an alignment refers to.
    /// Passes through any [`ErrorKind::Database`](crate::ErrorKind::Database) errors.
    /// Returns an error if an alignment cannot be decompressed.
    pub fn from_rows(database: &GAFBase, row_range: Range<usize>, graph: Option<&GBZ>) -> Result<Self> {
        let mut read_set = ReadSet { clusters: 1, ..Default::default() };

        // Build a record from the GAF-base, with the sequence possibly from the GBZ graph.
        let mut get_node = database.connection.prepare(
            "SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?1"
        )?;
        let mut get_record = |handle: usize| -> Result<GBZRecord> {
            // Get the edges and the BWT fragment from the GAF-base.
            let gaf_result = get_node.query_row(
                (handle,),
                |row: &Row<'_>| -> rusqlite::Result<(Vec<Pos>, Vec<u8>, Vec<u8>)> {
                    let edge_bytes: Vec<u8> = row.get(0)?;
                    let (edges, _) = Record::decompress_edges(&edge_bytes).unwrap();
                    let bwt: Vec<u8> = row.get(1)?;
                    let encoded_sequence: Vec<u8> = row.get(2)?;
                    let sequence = utils::decode_sequence(&encoded_sequence);
                    Ok((edges, bwt, sequence))
                }
            ).optional()?;
            if gaf_result.is_none() {
                return Err(Error::invalid_data(format!("GAF-base does not contain a record for node handle {}", handle)));
            }
            let (edges, bwt, mut sequence) = gaf_result.unwrap();
            if sequence.is_empty() {
                if let Some(graph) = graph {
                    let seq = graph.sequence(support::node_id(handle)).ok_or_else(||
                        Error::invalid_data(format!("The graph does not contain a sequence for node handle {}", handle))
                    )?;
                    sequence = seq.to_vec();
                } else {
                    return Err(Error::invalid_query("No reference provided for a reference-based GAF-base"));
                }
                if support::node_orientation(handle) == Orientation::Reverse {
                    sequence = support::reverse_complement(&sequence);
                }
            }

            unsafe {
                Ok(GBZRecord::from_raw_parts(handle, edges, bwt, sequence, None))
            }
        };

        // Get the reads in the given range of row ids.
        let mut get_reads = database.connection.prepare(
            "SELECT min_handle, max_handle, alignments, read_length, gbwt_starts, names, quality_strings, difference_strings, flags, numbers, optional
            FROM Alignments
            WHERE id >= ?1 AND id < ?2"
        )?;
        let mut rows = get_reads.query((row_range.start, row_range.end))?;
        while let Some(row) = rows.next()? {
            let block = GAFBaseInterface::decompress_block(row, 0)?;
            let block_size = block.len();
            for mut alignment in block {
                read_set.set_target_path_simple(&mut alignment, &mut get_record)?;
                if alignment.has_target_path() {
                    read_set.reads.push(alignment);
                    read_set.unclipped += 1;
                }
            }
            read_set.blocks += 1;
            read_set.candidates += block_size;
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

    /// Returns the number of candidate alignments in the decompressed blocks.
    #[inline]
    pub fn candidates(&self) -> usize {
        self.candidates
    }

    /// Returns the number of node records in the read set.
    ///
    /// Each record corresponds to an oriented node, and the opposite orientation may not be present.
    /// This includes all node records encountered while tracing the alignments, even when the alignment was not included in the read set.
    #[inline]
    pub fn node_records(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of node id clusters in the subgraph.
    #[inline]
    pub fn clusters(&self) -> usize {
        self.clusters
    }

    /// Returns an iterator over the reads in the set.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &Alignment> {
        self.reads.iter()
    }

    // Extracts the target sequence for the given alignment.
    fn target_sequence(&self, alignment: &Alignment) -> Result<Vec<u8>> {
        let target_path = alignment.target_path();
        if target_path.is_none() {
            return Ok(Vec::new());
        }
        let target_path = target_path.unwrap();

        let mut sequence = Vec::new();
        for handle in target_path{
            let record = self.nodes.get(handle);
            if record.is_none() {
                return Err(Error::invalid_data(format!("Read {}: Missing record for node handle {}", alignment.name, handle)));
            }
            let record = record.unwrap();
            sequence.extend_from_slice(record.sequence());
        }

        if sequence.len() != alignment.path_len {
            return Err(Error::invalid_data(format!(
                "Read {}: Target path length {} does not match the expected length {}",
                alignment.name, sequence.len(), alignment.path_len
            )));
        }
        Ok(sequence)
    }

    /// Serializes the read set in the GAF format.
    ///
    /// The output does not include any header lines, as the GAF file may consist of multiple read sets.
    /// Returns an error if the target sequence for a read is invalid or cannot be determined.
    /// Passes through any [`ErrorKind::Io`](crate::ErrorKind::Io) errors.
    pub fn to_gaf<W: Write>(&self, writer: &mut W) -> Result<()> {
        for alignment in self.reads.iter() {
            let target_sequence = self.target_sequence(alignment)?;
            let mut line = alignment.to_gaf(&target_sequence);
            line.push(b'\n');
            writer.write_all(&line)?;
        }
        Ok(())
    }
}

//-----------------------------------------------------------------------------
