//! A subgraph in a GBZ graph.
//!
//! This module provides functionality for extracting a subgraph around a specific position or interval of a specific path.
//! The subgraph contains all nodes within a given context and all edges between them.
//! All other paths within the subgraph can also be extracted, but they will not have any true metadata associated with them.

use crate::{GBZRecord, GBZPath, GraphInterface};
use crate::formats::{self, WalkMetadata, JSONValue};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap, HashMap};
use std::fmt::Display;
use std::io::{self, Write};
use std::ops::Range;

use gbwt::ENDMARKER;
use gbwt::{GBZ, GBWT, Orientation, Pos, FullPathName};

use gbwt::support;

use simple_sds::sparse_vector::{SparseVector, SparseBuilder};
use simple_sds::ops::PredSucc;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Output options for the haplotypes in the subgraph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HaplotypeOutput {
    /// Output all haplotypes as separate paths.
    All,

    /// Output only distinct haplotypes with the number of duplicates stored in the weight field.
    Distinct,

    /// Output only the reference path.
    ReferenceOnly,
}

impl Display for HaplotypeOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            HaplotypeOutput::All => write!(f, "all"),
            HaplotypeOutput::Distinct => write!(f, "distinct"),
            HaplotypeOutput::ReferenceOnly => write!(f, "reference only"),
        }
    }
}

// FIXME: There should be an enum (Offet<usize>, Interval<Range<usize>>, Node<usize>),
// but we first need to implement context extraction based on an interval.
/// Arguments for extracting a subgraph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubgraphQuery {
    // Name of the path to use as the reference path.
    path_name: Option<FullPathName>,

    // Position in the reference path (in bp).
    // Fields `path_name` and `offset` must be both `Some` or both `None`.
    offset: Option<usize>,

    // Node to be used as the reference position.
    // Fields `path_name` and `node` must be both `None` if this field is in use.
    node_id: Option<usize>,

    // Context size around the reference position (in bp).
    context: usize,

    // How to output the haplotypes.
    output: HaplotypeOutput,
}

impl SubgraphQuery {
    /// Creates a query that retrieves a subgraph around a path offset.
    ///
    /// # Arguments
    ///
    /// * `path_name`: Name of the reference path.
    /// * `offset`: Position in the reference path (in bp).
    /// * `context`: Context length around the reference position (in bp).
    /// * `output`: How to output the haplotypes.
    pub fn path_offset(path_name: &FullPathName, offset: usize, context: usize, output: HaplotypeOutput) -> Self {
        SubgraphQuery {
            path_name: Some(path_name.clone()),
            offset: Some(offset),
            node_id: None,
            context,
            output
        }
    }

    // FIXME: If interval length is odd, the right context is one bp too short.
    /// Cretes a query that retrieves a subgraph around a path interval.
    ///
    /// # Arguments
    ///
    /// * `path_name`: Name of the reference path.
    /// * `interval`: Interval of the reference path (in bp).
    /// * `context`: Context length around the reference interval (in bp).
    /// * `output`: How to output the haplotypes.
    pub fn path_interval(path_name: &FullPathName, interval: Range<usize>, context: usize, output: HaplotypeOutput) -> Self {
        let radius = interval.len() / 2;
        let midpoint = interval.start + radius;
        SubgraphQuery {
            path_name: Some(path_name.clone()),
            offset: Some(midpoint),
            node_id: None,
            context: context + radius,
            output
        }
    }

    /// Creates a query that retrieves a subgraph around a node.
    ///
    /// # Arguments
    ///
    /// * `node_id`: Identifier of the reference node.
    /// * `context`: Context length around the reference node (in bp).
    /// * `output`: How to output the haplotypes.
    ///
    /// # Panics
    ///
    /// Panics if `output` is [`HaplotypeOutput::ReferenceOnly`].
    pub fn node(node_id: usize, context: usize, output: HaplotypeOutput) -> Self {
        if output == HaplotypeOutput::ReferenceOnly {
            panic!("Cannot output a reference path in a node-based query");
        }
        SubgraphQuery {
            path_name: None,
            offset: None,
            node_id: Some(node_id),
            context,
            output
        }
    }

    /// Returns `true` if this query is based on a path offset.
    pub fn is_path_based(&self) -> bool {
        self.path_name.is_some() && self.offset.is_some()
    }

    // Returns `true` if this query is based on a node.
    pub fn is_node_based(&self) -> bool {
        self.node_id.is_some()
    }

    /// Returns the path name for this query, or [`None`] if the query does not use a reference path.
    pub fn path_name(&self) -> Option<&FullPathName> {
        self.path_name.as_ref()
    }

    /// Returns the offset for this query, or [`None`] if the query does not use a reference path.
    pub fn offset(&self) -> Option<usize> {
        self.offset
    }

    /// Returns the node identifier for this query, or [`None`] if the query does not use a node.
    pub fn node_id(&self) -> Option<usize> {
        self.node_id
    }

    /// Returns the context length (in bp) for the query.
    pub fn context(&self) -> usize {
        self.context
    }

    /// Returns the output format for the query.
    pub fn output(&self) -> HaplotypeOutput {
        self.output
    }
}

impl Display for SubgraphQuery {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_path_based() {
            write!(f, "(path {}, offset {}, context {}, {})", self.path_name().unwrap(), self.offset.unwrap(), self.context, self.output)
        } else if self.is_node_based() {
            write!(f, "(node {}, context {}, {})", self.node_id.unwrap(), self.context, self.output)
        } else {
            write!(f, "(invalid query)")
        }
    }
}

//-----------------------------------------------------------------------------

/// An index for random access to reference and generic paths in a GBZ graph.
///
/// Indexed paths are identified by their offsets in the index.
/// The offsets range from 0 to `path_count() - 1`.
///
/// The combination of [`GBZ`] and [`PathIndex`] is functionally similar to [`GraphInterface`] but tens of times faster.
/// An in-memory graph is better for batch operations, where the loading time is a negligible fraction of the total time.
/// The database is better for interactive applications, where the user works with relatively small subgraphs.
/// For a human graph, the database should be faster than the in-memory graph for subgraphs of up to 1 Mbp.
///
/// # Examples
///
/// ```
/// use gbz_base::PathIndex;
/// use gbwt::{GBZ, FullPathName, Pos, REF_SAMPLE};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let graph: GBZ = serialize::load_from(&filename).unwrap();
///
/// // Create a path index with 3 bp intervals.
/// let path_index = PathIndex::new(&graph, 3, false);
/// assert!(path_index.is_ok());
/// let path_index = path_index.unwrap();
///
/// // We have two components with one generic path in each.
/// assert_eq!(path_index.path_count(), 2);
/// assert_eq!(path_index.path_length(0), Some(5));
/// assert_eq!(path_index.path_length(1), Some(4));
///
/// // Consider the generic path in component A.
/// let path_name = FullPathName::generic("A");
/// let index_offset = path_index.find_path(&graph, &path_name);
/// assert_eq!(index_offset, Some(0));
/// let index_offset = index_offset.unwrap();
///
/// // There should be two indexed positions for the path.
/// let first_sample = path_index.indexed_position(index_offset, 2);
/// assert_eq!(first_sample, Some((0, Pos::new(22, 0))));
/// let second_sample = path_index.indexed_position(index_offset, 5);
/// assert_eq!(second_sample, Some((3, Pos::new(30, 0))));
/// let next_sample = path_index.indexed_position(index_offset, 100);
/// assert_eq!(next_sample, second_sample);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathIndex {
    // Maps path identifiers to index offsets.
    path_to_offset: HashMap<usize, usize>,

    // Maps index offsets to path identifiers.
    offset_to_path: HashMap<usize, usize>,

    // Sequence lengths for each path in bp.
    path_lengths: Vec<usize>,

    // Indexed sequence positions for each path.
    sequence_positions: Vec<SparseVector>,

    // GBWT positions corresponding to the indexed sequence positions.
    gbwt_positions: Vec<Vec<Pos>>,
}

impl PathIndex {
    /// Creates a new path index for the given GBZ graph.
    ///
    /// The index is built for all reference and generic paths.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ graph.
    /// * `interval`: Approximate distance between indexed positions (in bp).
    /// * `verbose`: Print progress information to stderr.
    pub fn new(graph: &GBZ, interval: usize, verbose: bool) -> Result<Self, String> {
        if verbose {
            eprintln!("Building path index");
        }

        let reference_paths = graph.reference_positions(interval, verbose);
        if reference_paths.is_empty() {
            return Err(String::from("No reference paths to index"));
        }

        let mut path_to_offset: HashMap<usize, usize> = HashMap::with_capacity(reference_paths.len());
        let mut offset_to_path: HashMap<usize, usize> = HashMap::with_capacity(reference_paths.len());
        let mut path_lengths: Vec<usize> = Vec::with_capacity(reference_paths.len());
        let mut sequence_positions = Vec::with_capacity(reference_paths.len());
        let mut gbwt_positions = Vec::with_capacity(reference_paths.len());
        for ref_path in reference_paths.iter() {
            path_to_offset.insert(ref_path.id, path_to_offset.len());
            offset_to_path.insert(offset_to_path.len(), ref_path.id);
            path_lengths.push(ref_path.len);
            let mut sequence = SparseBuilder::new(ref_path.len, ref_path.positions.len()).map_err(
                String::from
            )?;
            let mut gbwt = Vec::with_capacity(ref_path.positions.len());
            for (sequence_pos, gbwt_pos) in ref_path.positions.iter() {
                sequence.set(*sequence_pos);
                gbwt.push(*gbwt_pos);
            }
            sequence_positions.push(SparseVector::try_from(sequence).map_err(String::from)?);
            gbwt_positions.push(gbwt);
        }

        Ok(PathIndex { path_to_offset, offset_to_path, path_lengths, sequence_positions, gbwt_positions })
    }

    /// Returns the number of indexed paths.
    pub fn path_count(&self) -> usize {
        self.sequence_positions.len()
    }

    /// Returns the index offset for the indexed path with the given identifier.
    ///
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    pub fn path_to_offset(&self, path_id: usize) -> Option<usize> {
        self.path_to_offset.get(&path_id).cloned()
    }

    /// Returns the path identifier for the indexed path with the given index offset.
    ///
    /// Returns [`None`] if the path does not exist.
    pub fn offset_to_path(&self, index_offset: usize) -> Option<usize> {
        self.offset_to_path.get(&index_offset).cloned()
    }

    /// Returns the index offset for the path with the given metadata.
    ///
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    pub fn find_path(&self, graph: &GBZ, path_name: &FullPathName) -> Option<usize> {
        let metadata = graph.metadata()?;
        let path_id = metadata.find_path(path_name)?;
        self.path_to_offset(path_id)
    }

    /// Returns the length of the indexed path with the given index offset.
    ///
    /// Returns [`None`] if the path does not exist.
    pub fn path_length(&self, index_offset: usize) -> Option<usize> {
        self.path_lengths.get(index_offset).cloned()
    }

    /// Returns the last indexed position at or before `offset` on the path with name `path_name`.
    ///
    /// The return value consists of a sequence offset and a GBWT position.
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    /// This is similar to [`GraphInterface::find_path`] followed by [`GraphInterface::indexed_position`].
    ///
    /// # Arguments
    ///
    /// * `index_offset`: Offset of the path in the index.
    /// * `query_offset`: Sequence position in the path (in bp).
    pub fn indexed_position(&self, index_offset: usize, query_offset: usize) -> Option<(usize, Pos)> {
        let mut iter = self.sequence_positions[index_offset].predecessor(query_offset);
        if let Some((sample_offset, sequence_offset)) = iter.next() {
            let gbwt_pos = self.gbwt_positions[index_offset][sample_offset];
            Some((sequence_offset, gbwt_pos))
        } else {
            None
        }
    }
}

//-----------------------------------------------------------------------------

/// A subgraph based on a context around a graph position.
///
/// The position can be based on a path offset, a path interval, or a node identifier.
/// The path used for extracting the subgraph becomes the reference path for it.
/// Non-reference haplotypes do not have any metadata associated with them, as we cannot determine the identifier of a path from its GBWT position efficiently.
pub struct Subgraph {
    // Node records for the subgraph.
    records: BTreeMap<usize, GBZRecord>,

    // Paths in the subgraph.
    paths: Vec<PathInfo>,

    // Offset in `paths` for the reference path, if any.
    ref_id: Option<usize>,

    // Metadata for the reference path, if any.
    ref_path: Option<GBZPath>,

    // Interval of the reference path that is present in the subgraph, if any.
    ref_interval: Option<Range<usize>>,
}

// TODO: This could implement an interface similar to the node/edge part of GBZ.
impl Subgraph {
    /// Extracts a subgraph around the given query position.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ graph.
    /// * `path_index`: A path index for the graph, if the query is path-based.
    /// * `query`: Arguments for extracting the subgraph.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{PathIndex, Subgraph, SubgraphQuery, HaplotypeOutput};
    /// use gbwt::{GBZ, FullPathName};
    /// use gbwt::support;
    /// use simple_sds::serialize;
    ///
    /// // Get the graph.
    /// let gbz_file = support::get_test_data("example.gbz");
    /// let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    ///
    /// // Create a path index with 3 bp intervals.
    /// let path_index = PathIndex::new(&graph, 3, false).unwrap();
    ///
    /// // Extract a subgraph that contains an 1 bp context around path A offset 2.
    /// let path_name = FullPathName::generic("A");
    /// let query = SubgraphQuery::path_offset(&path_name, 2, 1, HaplotypeOutput::All);
    /// let subgraph = Subgraph::from_gbz(&graph, Some(&path_index), &query);
    /// assert!(subgraph.is_ok());
    /// let subgraph = subgraph.unwrap();
    ///
    /// // The subgraph should be centered around 1 bp node 14 of degree 4.
    /// assert_eq!(subgraph.node_count(), 5);
    /// assert_eq!(subgraph.path_count(), 3);
    ///
    /// // We get the same result using a node id.
    /// let query = SubgraphQuery::node(14, 1, HaplotypeOutput::All);
    /// let subgraph = Subgraph::from_gbz(&graph, None, &query);
    /// assert!(subgraph.is_ok());
    /// let subgraph = subgraph.unwrap();
    /// assert_eq!(subgraph.node_count(), 5);
    /// assert_eq!(subgraph.path_count(), 3);
    /// ```
    pub fn from_gbz(graph: &GBZ, path_index: Option<&PathIndex>, query: &SubgraphQuery) -> Result<Self, String> {
        // FIXME: Do interval-based and node-based properly.
        if query.is_path_based() {
            let path_index = path_index.ok_or(
                String::from("Path index is required for path-based queries")
            )?;
            let ref_path = GBZPath::with_name(graph, query.path_name().unwrap()).ok_or(
                format!("Cannot find path {}", query.path_name().unwrap())
            )?;
            let (query_pos, gbwt_pos) = query_position_from_gbz(graph, path_index, &ref_path, query.offset.unwrap())?;
            let records = extract_context_from_gbz(graph, query_pos, query.context)?;
            Self::path_based(records, ref_path, query.offset().unwrap(), query_pos, gbwt_pos, query.output())
        } else if query.is_node_based() {
            let node_id = query.node_id().unwrap();
            let node_len = graph.sequence_len(node_id).ok_or(
                format!("Cannot find node {}", node_id)
            )?;
            // FIXME: If node length is even, the right context is one bp too long.
            // What about odd length?
            let query_pos = GraphPosition {
                node: node_id,
                orientation: Orientation::Forward,
                offset: node_len / 2
            };
            let context = query.context + node_len / 2;
            let records = extract_context_from_gbz(graph, query_pos, context)?;
            Self::node_based(records, query.output())
        } else {
            Err(format!("Invalid query: {}", query))
        }
    }

    /// Extracts a subgraph around the given query position.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZBase, GraphInterface, Subgraph, SubgraphQuery, HaplotypeOutput};
    /// use gbwt::FullPathName;
    /// use gbwt::support;
    /// use simple_sds::serialize;
    /// use std::fs;
    ///
    /// // Create the database.
    /// let gbz_file = support::get_test_data("example.gbz");
    /// let db_file = serialize::temp_file_name("subgraph");
    /// let result = GBZBase::create_from_file(&gbz_file, &db_file);
    /// assert!(result.is_ok());
    ///
    /// // Open the database and create a graph interface.
    /// let database = GBZBase::open(&db_file).unwrap();
    /// let mut interface = GraphInterface::new(&database).unwrap();
    ///
    /// // Extract a subgraph that contains an 1 bp context around path A offset 2.
    /// let path_name = FullPathName::generic("A");
    /// let query = SubgraphQuery::path_offset(&path_name, 2, 1, HaplotypeOutput::All);
    /// let subgraph = Subgraph::from_db(&mut interface, &query);
    /// assert!(subgraph.is_ok());
    /// let subgraph = subgraph.unwrap();
    ///
    /// // The subgraph should be centered around 1 bp node 14 of degree 4.
    /// assert_eq!(subgraph.node_count(), 5);
    /// assert_eq!(subgraph.path_count(), 3);
    ///
    /// // Clean up.
    /// drop(interface);
    /// drop(database);
    /// fs::remove_file(&db_file).unwrap();
    /// ```
    pub fn from_db(graph: &mut GraphInterface, query: &SubgraphQuery) -> Result<Self, String> {
        // FIXME: Do interval-based and node-based properly.
        if query.is_path_based() {
            let ref_path = graph.find_path(query.path_name().unwrap())?;
            let ref_path = ref_path.ok_or(format!("Cannot find path {}", query.path_name().unwrap()))?;
            if !ref_path.is_indexed {
                return Err(format!("Path {} has not been indexed for random access", query.path_name().unwrap()));
            }
            let (query_pos, gbwt_pos) = query_position_from_db(graph, &ref_path, query.offset.unwrap())?;
            let records = extract_context_from_db(graph, query_pos, query.context)?;
            Self::path_based(records, ref_path, query.offset().unwrap(), query_pos, gbwt_pos, query.output())
        } else if query.is_node_based() {
            let node_id = query.node_id().unwrap();
            let record = graph.get_record(support::encode_node(node_id, Orientation::Forward))?;
            let record = record.ok_or(
                format!("Cannot find node {}", node_id)
            )?;
            // FIXME: If node length is odd, the right context is one bp too short.
            let node_len = record.sequence_len();
            let query_pos = GraphPosition {
                node: node_id,
                orientation: Orientation::Forward,
                offset: node_len / 2
            };
            let context = query.context + node_len / 2;
            let records = extract_context_from_db(graph, query_pos, context)?;
            Self::node_based(records, query.output())
        } else {
            Err(format!("Invalid query: {}", query))
        }
    }

    // Shared internal constructor for path-based queries.
    fn path_based(
        records: BTreeMap<usize, GBZRecord>,
        ref_path: GBZPath,
        offset: usize,
        query_pos: GraphPosition,
        gbwt_pos: Pos,
        output: HaplotypeOutput
    ) -> Result<Self, String> {
        let (mut paths, ref_id_offset) = extract_paths(&records, Some(gbwt_pos))?;
        let (mut ref_id, path_offset) = ref_id_offset.unwrap();
        let ref_start = offset - distance_to(&records, &paths[ref_id].path, path_offset, query_pos.offset);
        let ref_interval = ref_start..ref_start + paths[ref_id].len;

        if output == HaplotypeOutput::Distinct {
            // TODO: We should use bidirectional search in GBWT to find the distinct paths directly.
            let res = make_distinct(paths, Some(ref_id));
            paths = res.0;
            ref_id = res.1.unwrap();
        } else if output == HaplotypeOutput::ReferenceOnly {
            paths = vec![paths[ref_id].clone()];
            ref_id = 0;
        }

        Ok(Subgraph {
            records,
            paths,
            ref_id: Some(ref_id),
            ref_path: Some(ref_path),
            ref_interval: Some(ref_interval),
        })
    }

    // Shared internal constructor for node-based queries.
    fn node_based(
        records: BTreeMap<usize, GBZRecord>,
        output: HaplotypeOutput
    ) -> Result<Self, String> {
        let (mut paths, _) = extract_paths(&records, None)?;
        if output == HaplotypeOutput::Distinct {
            let res = make_distinct(paths, None);
            paths = res.0;
        }
        Ok(Subgraph {
            records,
            paths,
            ref_id: None,
            ref_path: None,
            ref_interval: None,
        })
    }

    /// Returns the number of nodes in the subgraph.
    pub fn node_count(&self) -> usize {
        self.records.len() / 2
    }

    /// Returns the number of paths in the subgraph.
    pub fn path_count(&self) -> usize {
        self.paths.len()
    }

    // Returns the total length of the nodes in the given path interval.
    fn path_len(&self, path: &[usize], interval: Range<usize>) -> usize {
        let mut result = 0;
        for i in interval {
            let record = self.records.get(&path[i]).unwrap();
            result += record.sequence_len();
        }
        result
    }

    // Appends a new edit or extends the previous one.
    fn append_edit(edits: &mut Vec<(EditOperation, usize)>, op: EditOperation, len: usize) {
        if let Some((prev_op, prev_len)) = edits.last_mut() {
            if *prev_op == op {
                *prev_len += len;
                return;
            }
        }
        edits.push((op, len));
    }

    // Appends the relevant edits for two path intervals, which are assumed to be diverging.
    fn append_edits(&self, edits: &mut Vec<(EditOperation, usize)>, path: &[usize], path_interval: Range<usize>, ref_path: &[usize], ref_interval: Range<usize>) {
        let path_len = self.path_len(path, path_interval);
        let ref_len = self.path_len(ref_path, ref_interval);
        if path_len == ref_len && path_len > 0 && path_len < 5 {
            Self::append_edit(edits, EditOperation::Match, path_len);
        } else {
            if path_len > 0 {
                Self::append_edit(edits, EditOperation::Insertion, path_len);
            }
            if ref_len > 0 {
                Self::append_edit(edits, EditOperation::Deletion, ref_len);
            }
        }
    }

    // Returns the CIGAR string for the given path, aligned to the reference path.
    // Takes the alignment from the longest common subsequence of the node sequences.
    // Diverging parts that are of equal length and at most 4 bp are represented as
    // matches. Otherwise they are represented as an insertion and/or a deletion.
    fn align_to_ref(&self, path_id: usize) -> Option<String> {
        let ref_id = self.ref_id?;
        if path_id == ref_id || path_id >= self.paths.len() {
            return None;
        }

        // Fill in the LCS matrix.
        // TODO: Something more efficient?
        let ref_path = &self.paths[ref_id].path;
        let path = &self.paths[path_id].path;
        let mut dp_matrix = vec![vec![0; ref_path.len() + 1]; path.len() + 1];
        for (path_offset, path_value) in path.iter().enumerate() {
            for (ref_offset, ref_value) in ref_path.iter().enumerate() {
                if path_value == ref_value {
                    dp_matrix[path_offset + 1][ref_offset + 1] = dp_matrix[path_offset][ref_offset] + 1;
                } else {
                    dp_matrix[path_offset + 1][ref_offset + 1] = dp_matrix[path_offset][ref_offset + 1].max(dp_matrix[path_offset + 1][ref_offset]);
                }
            }
        }

        // Trace back the LCS.
        let mut lcs: Vec<(usize, usize)> = Vec::new();
        let mut path_offset = path.len();
        let mut ref_offset = ref_path.len();
        while path_offset > 0 && ref_offset > 0 {
            if path[path_offset - 1] == ref_path[ref_offset - 1] {
                lcs.push((path_offset - 1, ref_offset - 1));
                path_offset -= 1;
                ref_offset -= 1;
            } else if dp_matrix[path_offset - 1][ref_offset] > dp_matrix[path_offset][ref_offset - 1] {
                path_offset -= 1;
            } else {
                ref_offset -= 1;
            }
        }
        lcs.reverse();

        // Convert the LCS to a sequence of edit operations
        let mut edits: Vec<(EditOperation, usize)> = Vec::new();
        path_offset = 0;
        ref_offset = 0;
        for (next_path_offset, next_ref_offset) in lcs.iter() {
            self.append_edits(&mut edits, path, path_offset..*next_path_offset, ref_path, ref_offset..*next_ref_offset);
            let node_len = self.records.get(&path[*next_path_offset]).unwrap().sequence_len();
            Self::append_edit(&mut edits, EditOperation::Match, node_len);
            path_offset = next_path_offset + 1;
            ref_offset = next_ref_offset + 1;
        }
        self.append_edits(&mut edits, path, path_offset..path.len(), ref_path, ref_offset..ref_path.len());

        // Convert the edits to a CIGAR string.
        let mut result = String::new();
        for (op, len) in edits.iter() {
            result.push_str(&format!("{}{}", len, op));
        }
        Some(result)
    }

    /// Writes the subgraph in the GFA format to the given output.
    ///
    /// If `cigar` is true, the CIGAR strings for the non-reference haplotypes are included in the output.
    pub fn write_gfa<T: Write>(&self, output: &mut T, cigar: bool) -> io::Result<()> {
        // Header.
        let reference_samples = self.ref_path.as_ref().map(|path| path.name.sample.as_ref());
        formats::write_gfa_header(reference_samples, output)?;

        // Segments.
        for (handle, record) in self.records.iter() {
            if support::node_orientation(*handle) == Orientation::Forward {
                formats::write_gfa_node(record.id(), record.sequence(), output)?;
            }
        }

        // Links.
        for (handle, record) in self.records.iter() {
            let from = support::decode_node(*handle);
            for successor in record.successors() {
                let to = support::decode_node(successor);
                if self.records.contains_key(&successor) && support::edge_is_canonical(from, to) {
                    formats::write_gfa_link(
                        (from.0.to_string().as_bytes(), from.1),
                        (to.0.to_string().as_bytes(), to.1),
                        output
                    )?;
                }
            }
        }

        // Paths.
        if let Some((metadata, ref_id)) = self.ref_metadata() {
            formats::write_gfa_walk(&self.paths[ref_id].path, &metadata, output)?;
        }
        let mut haplotype = 1;
        let contig_name = self.contig_name();
        for (id, path_info) in self.paths.iter().enumerate() {
            if Some(id) == self.ref_id {
                continue;
            }
            let mut metadata = WalkMetadata::anonymous(haplotype, &contig_name, path_info.len);
            metadata.add_weight(path_info.weight);
            if cigar {
                metadata.add_cigar(self.align_to_ref(id));
            }
            formats::write_gfa_walk(&path_info.path, &metadata, output)?;
            haplotype += 1;
        }

        Ok(())
    }

    /// Writes the subgraph in the JSON format to the given output.
    ///
    /// If `cigar` is true, the CIGAR strings for the non-reference haplotypes are included in the output.
    pub fn write_json<T: Write>(&self, output: &mut T, cigar: bool) -> io::Result<()> {
        // Nodes.
        let mut nodes: Vec<JSONValue> = Vec::new();
        for (_, record) in self.records.iter() {
            let (id, orientation) = support::decode_node(record.handle());
            if orientation == Orientation::Reverse {
                continue;
            }
            let node = JSONValue::Object(vec![
                ("id".to_string(), JSONValue::String(id.to_string())),
                ("sequence".to_string(), JSONValue::String(String::from_utf8_lossy(record.sequence()).to_string())),
            ]);
            nodes.push(node);
        }

        // Edges.
        let mut edges: Vec<JSONValue> = Vec::new();
        for (handle, record) in self.records.iter() {
            let from = support::decode_node(*handle);
            for successor in record.successors() {
                let to = support::decode_node(successor);
                if self.records.contains_key(&successor) && support::edge_is_canonical(from, to) {
                    let edge = JSONValue::Object(vec![
                        ("from".to_string(), JSONValue::String(from.0.to_string())),
                        ("from_is_reverse".to_string(), JSONValue::Boolean(from.1 == Orientation::Reverse)),
                        ("to".to_string(), JSONValue::String(to.0.to_string())),
                        ("to_is_reverse".to_string(), JSONValue::Boolean(to.1 == Orientation::Reverse)),
                    ]);
                    edges.push(edge);
                }
            }
        }

        // Paths.
        let mut paths: Vec<JSONValue> = Vec::new();
        if let Some((metadata, ref_id)) = self.ref_metadata() {
            let ref_path = formats::json_path(&self.paths[ref_id].path, &metadata);
            paths.push(ref_path);
        }
        let mut haplotype = 1;
        let contig_name = self.contig_name();
        for (id, path_info) in self.paths.iter().enumerate() {
            if Some(id) == self.ref_id {
                continue;
            }
            let mut metadata = WalkMetadata::anonymous(haplotype, &contig_name, path_info.len);
            metadata.add_weight(path_info.weight);
            if cigar {
                metadata.add_cigar(self.align_to_ref(id));
            }
            let path = formats::json_path(&path_info.path, &metadata);
            paths.push(path);
            haplotype += 1;
        }

        let result = JSONValue::Object(vec![
            ("nodes".to_string(), JSONValue::Array(nodes)),
            ("edges".to_string(), JSONValue::Array(edges)),
            ("paths".to_string(), JSONValue::Array(paths)),
        ]);
        output.write_all(result.to_string().as_bytes())?;

        Ok(())
    }

    // Builds metadata for the reference path.
    fn ref_metadata(&self) -> Option<(WalkMetadata, usize)> {
        let ref_id = self.ref_id?;
        let ref_path = self.ref_path.as_ref()?;
        let interval = self.ref_interval.as_ref()?.clone();
        let mut metadata = WalkMetadata::path_interval(ref_path.as_ref(), interval);
        metadata.add_weight(self.paths[ref_id].weight);
        Some((metadata, ref_id))
    }

    // Determines a contig name for the subgraph.
    fn contig_name(&self) -> String {
        if let Some(ref_path) = self.ref_path.as_ref() {
            ref_path.name.contig.clone()
        } else {
            String::from("unknown")
        }
    }
}

//-----------------------------------------------------------------------------

// Returns the graph position and the GBWT position for the given offset.
fn query_position_from_db(graph: &mut GraphInterface, path: &GBZPath, query_offset: usize) -> Result<(GraphPosition, Pos), String> {
    let result = graph.indexed_position(path.handle, query_offset)?;
    let (mut path_offset, mut pos) = result.ok_or(
        format!("Path {} has not been indexed for random access", path.name())
    )?;

    let mut graph_pos: Option<GraphPosition> = None;
    let mut gbwt_pos: Option<Pos> = None;
    loop {
        let handle = pos.node;
        let record = graph.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset + record.sequence_len() > query_offset {
            graph_pos = Some(GraphPosition {
                node: support::node_id(handle),
                orientation: support::node_orientation(handle),
                offset: query_offset - path_offset,
            });
            gbwt_pos = Some(pos);
            break;
        }
        path_offset += record.sequence_len();
        let gbwt_record = record.to_gbwt_record();
        let next = gbwt_record.lf(pos.offset);
        if next.is_none() {
            break;
        }
        pos = next.unwrap();
    }

    let graph_pos = graph_pos.ok_or(
        format!("Path {} does not contain offset {}", path.name(), query_offset)
    )?;
    let gbwt_pos = gbwt_pos.unwrap();
    Ok((graph_pos, gbwt_pos))
}


// Returns the graph position and the GBWT position for the given offset.
fn query_position_from_gbz(graph: &GBZ, path_index: &PathIndex, path: &GBZPath, query_offset: usize) -> Result<(GraphPosition, Pos), String> {
    // Path id to an indexed position.
    let index_offset = path_index.path_to_offset(path.handle).ok_or(
        format!("Path {} has not been indexed for random access", path.name())
    )?;
    let (mut path_offset, mut pos) = path_index.indexed_position(index_offset, query_offset).unwrap();

    let mut graph_pos: Option<GraphPosition> = None;
    let mut gbwt_pos: Option<Pos> = None;
    let index: &GBWT = graph.as_ref();
    loop {
        let node_id = support::node_id(pos.node);
        let node_len = graph.sequence_len(node_id).unwrap();
        if path_offset + node_len > query_offset {
            graph_pos = Some(GraphPosition {
                node: node_id,
                orientation: support::node_orientation(pos.node),
                offset: query_offset - path_offset,
            });
            gbwt_pos = Some(pos);
            break;
        }
        path_offset += node_len;
        let next = index.forward(pos);
        if next.is_none() {
            break;
        }
        pos = next.unwrap();
    }

    let graph_pos = graph_pos.ok_or(
        format!("Path {} does not contain offset {}", path.name(), query_offset)
    )?;
    let gbwt_pos = gbwt_pos.unwrap();
    Ok((graph_pos, gbwt_pos))
}

//-----------------------------------------------------------------------------

fn distance_to_end(record: &GBZRecord, from: GraphPosition, orientation: Orientation) -> usize {
    if orientation == from.orientation {
        record.sequence_len() - from.offset
    } else {
        from.offset + 1
    }
}

fn extract_context(
    from: GraphPosition,
    context: usize,
    mut get_record: impl FnMut(usize) -> Result<GBZRecord, String>
) -> Result<BTreeMap<usize, GBZRecord>, String> {
    // Start graph traversal from the initial node.
    let mut active: BinaryHeap<Reverse<(usize, GraphPosition)>> = BinaryHeap::new(); // (distance, node id)
    active.push(Reverse((0, from)));

    // Traverse in both directions.
    let mut selected: BTreeMap<usize, GBZRecord> = BTreeMap::new();
    while !active.is_empty() {
        let (distance, position) = active.pop().unwrap().0;
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let handle = support::encode_node(position.node, orientation);
            if selected.contains_key(&handle) {
                continue;
            }
            let record = get_record(handle)?;
            let next_distance = distance + distance_to_end(&record, position, orientation);
            if next_distance <= context {
                for successor in record.successors() {
                    if !selected.contains_key(&successor) {
                        let next = GraphPosition {
                            node: support::node_id(successor),
                            orientation: support::node_orientation(successor),
                            offset: 0,
                        };
                        active.push(Reverse((next_distance, next)));
                    }
                }
            }
            selected.insert(handle, record);
        }
    }

    Ok(selected)
}

// Returns all node records within `context` bp around the given position.
fn extract_context_from_gbz(
    graph: &GBZ,
    from: GraphPosition,
    context: usize
) -> Result<BTreeMap<usize, GBZRecord>, String> {
    extract_context(
        from, context,
        |handle| GBZRecord::from_gbz(graph, handle).ok_or(
            format!("The graph does not contain handle {}", handle)
        )
    )
}

// Returns all node records within `context` bp around the given position.
fn extract_context_from_db(
    graph: &mut GraphInterface,
    from: GraphPosition,
    context: usize
) -> Result<BTreeMap<usize, GBZRecord>, String> {
    extract_context(
        from, context,
        |handle| {
            let record = graph.get_record(handle)?;
            record.ok_or(format!("The graph does not contain handle {}", handle))
        }
    )
}

//-----------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PathInfo {
    path: Vec<usize>,
    len: usize,
    weight: Option<usize>,
}

impl PathInfo {
    fn new(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: None }
    }

    fn weighted(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: Some(1) }
    }

    fn increment_weight(&mut self) {
        if let Some(weight) = self.weight {
            self.weight = Some(weight + 1);
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum EditOperation {
    Match,
    Insertion,
    Deletion,
}

impl Display for EditOperation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EditOperation::Match => write!(f, "M"),
            EditOperation::Insertion => write!(f, "I"),
            EditOperation::Deletion => write!(f, "D"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct GraphPosition {
    node: usize,
    orientation: Orientation,
    offset: usize,
}

//-----------------------------------------------------------------------------

fn next_pos(pos: Pos, successors: &BTreeMap<usize, Vec<(Pos, bool)>>) -> Option<Pos> {
    if let Some(v) = successors.get(&pos.node) {
        let (next, _) = v[pos.offset];
        if next.node == ENDMARKER || !successors.contains_key(&next.node) {
            None
        } else {
            Some(next)
        }
    } else {
        None
    }
}

// Extract all paths in the subgraph. If a GBWT position for the reference path is
// specified, the second return value is (offset in result, offset on that path) for
// the handle corresponding to `ref_pos`.
fn extract_paths(
    subgraph: &BTreeMap<usize, GBZRecord>,
    ref_pos: Option<Pos>
) -> Result<(Vec<PathInfo>, Option<(usize, usize)>), String> {
    // Decompress the GBWT node records for the subgraph.
    let mut keys: Vec<usize> = Vec::new();
    let mut successors: BTreeMap<usize, Vec<(Pos, bool)>> = BTreeMap::new();
    for (handle, gbz_record) in subgraph.iter() {
        let gbwt_record = gbz_record.to_gbwt_record();
        let decompressed: Vec<(Pos, bool)> = gbwt_record.decompress().into_iter().map(|x| (x, false)).collect();
        keys.push(*handle);
        successors.insert(*handle, decompressed);
    }

    // Mark the positions that have predecessors in the subgraph.
    for handle in keys.iter() {
        let decompressed = successors.get(handle).unwrap().clone();
        for (pos, _) in decompressed.iter() {
            if let Some(v) = successors.get_mut(&pos.node) {
                v[pos.offset].1 = true;
            }
        }
    }

    // FIXME: Check for infinite loops.
    // Extract all paths and note if one of them passes through `ref_pos`.
    let mut result: Vec<PathInfo> = Vec::new();
    let mut ref_id_offset: Option<(usize, usize)> = None;
    for (handle, positions) in successors.iter() {
        for (offset, (_, has_predecessor)) in positions.iter().enumerate() {
            if *has_predecessor {
                continue;
            }
            let mut curr = Some(Pos::new(*handle, offset));
            let mut is_ref = false;
            let mut path: Vec<usize> = Vec::new();
            let mut len = 0;
            while let Some(pos) = curr {
                if Some(pos) == ref_pos {
                    ref_id_offset = Some((result.len(), path.len()));
                    is_ref = true;
                }
                path.push(pos.node);
                len += subgraph.get(&pos.node).unwrap().sequence_len();
                curr = next_pos(pos, &successors);
            }
            if is_ref {
                if !support::encoded_path_is_canonical(&path) {
                    eprintln!("Warning: the reference path is not in canonical orientation");
                }
                result.push(PathInfo::new(path, len));
            } else if support::encoded_path_is_canonical(&path) {
                result.push(PathInfo::new(path, len));
            }
        }
    }

    if ref_pos.is_some() && ref_id_offset.is_none() {
        return Err(String::from("Could not find the reference path"));
    }
    Ok((result, ref_id_offset))
}

fn distance_to(subgraph: &BTreeMap<usize, GBZRecord>, path: &[usize], path_offset: usize, node_offset: usize) -> usize {
    let mut result = node_offset;
    for handle in path.iter().take(path_offset) {
        result += subgraph.get(handle).unwrap().sequence_len();
    }
    result
}

// Returns all distinct paths and uses the weight field for storing their counts.
// Also updates `ref_id` if given.
// TODO: Use hashing?
fn make_distinct(
    paths: Vec<PathInfo>,
    ref_id: Option<usize>
) -> (Vec<PathInfo>, Option<usize>) {
    let ref_path = ref_id.map(|id| paths[id].path.clone());
    let mut paths = paths;
    paths.sort_unstable();

    let mut new_paths: Vec<PathInfo> = Vec::new();
    let mut ref_id = None;
    for info in paths.iter() {
        if new_paths.is_empty() || new_paths.last().unwrap().path != info.path {
            if let Some(ref_path) = &ref_path {
                if info.path == *ref_path {
                    ref_id = Some(new_paths.len());
                }
            }
            new_paths.push(PathInfo::weighted(info.path.clone(), info.len));
        } else {
            new_paths.last_mut().unwrap().increment_weight();
        }
    }

    (new_paths, ref_id)
}

//-----------------------------------------------------------------------------
