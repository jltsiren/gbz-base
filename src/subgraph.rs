//! A subgraph in a GBZ graph.
//!
//! This module provides functionality for extracting a subgraph around a specific position or interval of a specific path.
//! The subgraph contains all nodes within a given context and all edges between them.
//! All other paths within the subgraph can also be extracted, but they will not have any true metadata associated with them.

use crate::{GBZRecord, GBZPath, GraphInterface};
use crate::formats::{self, WalkMetadata, JSONValue};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap, BTreeSet, HashMap};
use std::fmt::Display;
use std::io::{self, Write};
use std::ops::Range;
use std::cmp;

use gbwt::ENDMARKER;
use gbwt::{GBZ, GBWT, GraphPosition, Orientation, Pos, FullPathName};

use gbwt::{algorithms, support};

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

// FIXME: There should be an enum (Offet(usize), Interval(Range<usize>), Node(usize)),
// FIXME: Or maybe Nodes(Vec<usize>)
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
    /// The reference path should be specified by using a sample name, a contig name, and optionally a haplotype number.
    /// The fragment field should not be used.
    /// If the reference haplotype is fragmented, the query will try to find the right fragment.
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
    /// The reference path should be specified by using a sample name, a contig name, and optionally a haplotype number.
    /// The fragment field should not be used.
    /// If the reference haplotype is fragmented, the query will try to find the right fragment.
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

    // Returns a path name containing the query offset, assuming that this is a path-based query.
    fn path_name_with_offset(&self) -> FullPathName {
        let mut result = self.path_name().unwrap().clone();
        result.fragment = self.offset().unwrap();
        result
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

        // NOTE: For fragmented haplotypes, the reference positions are relative
        // to the start of the fragment.
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

// FIXME: Tests and examples for the new functionality.

// FIXME: Document the intended use of this struct, for example as a sliding window over a reference path.
/// A subgraph based on a context around a graph position.
///
/// The position can be based on a path offset, a path interval, or a node identifier.
/// The path used for extracting the subgraph becomes the reference path for it.
/// Non-reference haplotypes do not have any metadata associated with them, as we cannot determine the identifier of a path from its GBWT position efficiently.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Subgraph {
    // Node records for the subgraph.
    records: BTreeMap<usize, GBZRecord>,

    // Paths in the subgraph.
    paths: Vec<PathInfo>,

    // Offset in `paths` for the reference path, if any.
    ref_id: Option<usize>,

    // Metadata for the reference path, if any.
    ref_path: Option<FullPathName>,

    // Interval of the reference path that is present in the subgraph, if any.
    ref_interval: Option<Range<usize>>,
}

//-----------------------------------------------------------------------------

// TODO: This could implement an interface similar to the node/edge part of GBZ.
/// Construction.
impl Subgraph {
    /// Creates a new empty subgraph.
    pub fn new() -> Self {
        Subgraph::default()
    }

    // FIXME: also for a node and a path with an interval
    // FIXME: Make this report the number of inserted and removed nodes.
    /// Updates the subgraph to a context around the given graph position.
    ///
    /// Reuses existing records when possible.
    /// Removes all path information from the subgraph.
    ///
    /// # Arguments
    ///
    /// * `pos`: The reference position for the subgraph.
    /// * `context`: The context length around the reference position (in bp).
    /// * `get_record`: A function that returns the record for the given node handle.
    ///
    /// # Errors
    ///
    /// Passes through any errors from `get_record`.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZRecord, Subgraph, HaplotypeOutput};
    /// use gbwt::{GBZ, FullPathName, GraphPosition, Orientation};
    /// use gbwt::support;
    /// use simple_sds::serialize;
    ///
    /// // Get the graph.
    /// let gbz_file = support::get_test_data("example.gbz");
    /// let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    ///
    /// // Extract a subgraph that contains an 1 bp context around node 14.
    /// let pos = GraphPosition::new(14, Orientation::Forward, 0);
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.with_context(pos, 1, &mut |handle| {
    ///     GBZRecord::from_gbz(&graph, handle).ok_or(
    ///         format!("The graph does not contain handle {}", handle)
    ///     )
    /// });
    /// assert!(result.is_ok());
    ///
    /// // The subgraph should be centered around 1 bp node 14 of degree 4.
    /// let true_nodes = [12, 13, 14, 15, 16];
    /// assert_eq!(subgraph.node_count(), true_nodes.len());
    /// let node_ids = subgraph.node_ids();
    /// assert!(node_ids.iter().eq(true_nodes.iter()));
    ///
    /// // But there are no paths.
    /// assert_eq!(subgraph.path_count(), 0);
    /// ```
    pub fn with_context(
        &mut self,
        pos: GraphPosition,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        // Start graph traversal from the initial node.
        let mut active: BinaryHeap<Reverse<(usize, GraphPosition)>> = BinaryHeap::new(); // (distance, node id)
        active.push(Reverse((0, pos)));

        self.insert_context(active, context, get_record)
    }

    // Returns the distance from the given position to the successors of the node in the given orientation.
    fn distance_to_successors(record: &GBZRecord, from: GraphPosition, orientation: Orientation) -> usize {
        if orientation == from.orientation {
            record.sequence_len() - from.offset
        } else {
            from.offset + 1
        }
    }

    // FIXME: check that if we call this twice, we do not add or remove any records.
    // Inserts all nodes within the context around the active positions into the subgraph.
    // Reuses existing records if possible.
    // Removes existing nodes that are not in the new context.
    // Clears all path information.
    fn insert_context(
        &mut self,
        active: BinaryHeap<Reverse<(usize, GraphPosition)>>,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        self.clear_paths();

        let mut active = active;
        let mut visited: BTreeSet<usize> = BTreeSet::new();
        let mut to_remove = self.node_ids();

        while !active.is_empty() {
            let (distance, position) = active.pop().unwrap().0;
            if visited.contains(&position.node) {
                continue;
            }
            visited.insert(position.node);
            to_remove.remove(&position.node);
            if !self.has_node(position.node) {
                self.add_node_internal(position.node, get_record)?;
            }
            for orientation in [Orientation::Forward, Orientation::Reverse] {
                let handle = support::encode_node(position.node, orientation);
                let record = self.record(handle).unwrap();
                let next_distance = distance + Self::distance_to_successors(record, position, orientation);
                if next_distance <= context {
                    for successor in record.successors() {
                        if !visited.contains(&support::node_id(successor)) {
                            let next = GraphPosition::new(
                                support::node_id(successor),
                                support::node_orientation(successor),
                                0
                            );
                            active.push(Reverse((next_distance, next)));
                        }
                    }
                }
            }
        }

        for node_id in to_remove {
            self.remove_node_internal(node_id);
        }

        Ok(())
    }

    // FIXME: from_gbz and from_db are very similar.
    /// Extracts a subgraph around the given query position.
    ///
    /// Reuses existing records when possible.
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
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.from_gbz(&graph, Some(&path_index), &query);
    /// assert!(result.is_ok());
    ///
    /// // The subgraph should be centered around 1 bp node 14 of degree 4.
    /// assert_eq!(subgraph.node_count(), 5);
    /// assert_eq!(subgraph.path_count(), 3);
    ///
    /// // We get the same result using a node id.
    /// let query = SubgraphQuery::node(14, 1, HaplotypeOutput::All);
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.from_gbz(&graph, None, &query);
    /// assert!(result.is_ok());
    /// assert_eq!(subgraph.node_count(), 5);
    /// assert_eq!(subgraph.path_count(), 3);
    /// ```
    pub fn from_gbz(&mut self, graph: &GBZ, path_index: Option<&PathIndex>, query: &SubgraphQuery) -> Result<(), String> {
        // These depend on the query type.
        let query_pos;
        let mut context = query.context;
        let mut reference_path = None;

        // FIXME: Do interval-based and node-based properly.
        if query.is_path_based() {
            let path_index = path_index.ok_or(
                String::from("Path index is required for path-based queries")
            )?;
            let path_name = query.path_name_with_offset();
            let ref_path = GBZPath::with_name(graph, &path_name).ok_or(
                format!("Cannot find a path covering {}", path_name)
            )?;
            // Transform the offset relative to the haplotype to the offset relative to the path.
            let query_offset = path_name.fragment - ref_path.name.fragment;
            let ref_pos = query_position_from_gbz(graph, path_index, &ref_path, query_offset)?;
            query_pos = ref_pos.graph_pos;
            reference_path = Some((ref_pos, ref_path.name));
        } else if query.is_node_based() {
            if query.output() == HaplotypeOutput::ReferenceOnly {
                return Err(String::from("Cannot output a reference path in a node-based query"));
            }
            let node_id = query.node_id().unwrap();
            let node_len = graph.sequence_len(node_id).ok_or(
                format!("Cannot find node {}", node_id)
            )?;
            // FIXME: If node length is even, the right context is one bp too long.
            // What about odd length?
            query_pos = GraphPosition::new(node_id, Orientation::Forward, node_len / 2);
            context = query.context + node_len / 2;
        } else {
            return Err(format!("Invalid query: {}", query));
        }

        // Now we can build the subgraph.
        self.with_context(query_pos, context, &mut |handle| {
            GBZRecord::from_gbz(graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
        })?;
        self.extract_paths(reference_path, query.output())?;

        Ok(())
    }

    /// Extracts a subgraph around the given query position.
    ///
    /// Reuses existing records when possible.
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
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.from_db(&mut interface, &query);
    /// assert!(result.is_ok());
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
    pub fn from_db(&mut self, graph: &mut GraphInterface, query: &SubgraphQuery) -> Result<(), String> {
        // These depend on the query type.
        let query_pos;
        let mut context = query.context;
        let mut reference_path = None;

        // FIXME: Do interval-based and node-based properly.
        if query.is_path_based() {
            let path_name = query.path_name_with_offset();
            let ref_path = graph.find_path(&path_name)?;
            let ref_path = ref_path.ok_or(format!("Cannot find a path covering {}", path_name))?;
            if !ref_path.is_indexed {
                return Err(format!("Path {} has not been indexed for random access", ref_path.name));
            }
            // Transform the offset relative to the haplotype to the offset relative to the path.
            let query_offset = path_name.fragment - ref_path.name.fragment;
            let ref_pos = query_position_from_db(graph, &ref_path, query_offset)?;
            query_pos = ref_pos.graph_pos;
            reference_path = Some((ref_pos, ref_path.name));
        } else if query.is_node_based() {
            if query.output() == HaplotypeOutput::ReferenceOnly {
                return Err(String::from("Cannot output a reference path in a node-based query"));
            }
            let node_id = query.node_id().unwrap();
            let record = graph.get_record(support::encode_node(node_id, Orientation::Forward))?;
            let record = record.ok_or(
                format!("Cannot find node {}", node_id)
            )?;
            // FIXME: If node length is odd, the right context is one bp too short.
            let node_len = record.sequence_len();
            query_pos = GraphPosition::new(node_id, Orientation::Forward, node_len / 2);
            context = query.context + node_len / 2;
        } else {
            return Err(format!("Invalid query: {}", query));
        }

        // Now we can build the subgraph.
        self.with_context(query_pos, context, &mut |handle| {
            let record = graph.get_record(handle)?;
            record.ok_or(format!("The graph does not contain handle {}", handle))
        })?;
        self.extract_paths(reference_path, query.output())?;

        Ok(())
    }

    // Returns the successor position for the given GBWT position, if it is in the subgraph.
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

    /// Extracts all paths in the subgraph.
    ///
    /// Clears the current paths in the subgraph.
    /// If a reference path is given, one of the paths is assumed to visit the position.
    /// This will then set all information related to the reference path.
    ///
    /// # Arguments
    ///
    /// * `reference_path`: Position on the reference path and the name of the path.
    /// * `output`: How to output the haplotypes.
    ///
    /// # Errors
    ///
    /// Returns an error if no path visits the reference position.
    /// Returns an error if reference-only output is requested without a reference path.
    /// Clears all path information in case of an error.
    pub fn extract_paths(
        &mut self,
        reference_path: Option<(PathPosition, FullPathName)>,
        output: HaplotypeOutput
    ) -> Result<(), String> {
        self.clear_paths();

        let ref_pos;
        if let Some((position, name)) = reference_path {
            ref_pos = Some(position);
            self.ref_path = Some(name);
        } else {
            ref_pos = None;
        }

        // Decompress the GBWT node records for the subgraph.
        let mut keys: Vec<usize> = Vec::new();
        let mut successors: BTreeMap<usize, Vec<(Pos, bool)>> = BTreeMap::new();
        for (handle, gbz_record) in self.records.iter() {
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
        // `ref_offset` is the offset of the node containing `ref_pos`.
        let mut ref_offset: Option<usize> = None;
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
                    if let Some(position) = ref_pos.as_ref() {
                        if pos == position.gbwt_pos {
                            self.ref_id = Some(self.paths.len());
                            ref_offset = Some(path.len());
                            is_ref = true;
                        }
                    }
                    path.push(pos.node);
                    len += self.records.get(&pos.node).unwrap().sequence_len();
                    curr = Self::next_pos(pos, &successors);
                }
                if is_ref {
                    if !support::encoded_path_is_canonical(&path) {
                        eprintln!("Warning: the reference path is not in canonical orientation");
                    }
                    self.paths.push(PathInfo::new(path, len));
                } else if support::encoded_path_is_canonical(&path) {
                    self.paths.push(PathInfo::new(path, len));
                }
            }
        }

        // Now we can set the reference interval.
        if ref_pos.is_some() {
            if let Some(offset) = ref_offset {
                let ref_info = &self.paths[self.ref_id.unwrap()];
                self.ref_interval = Some(ref_info.path_interval(self, offset, ref_pos.as_ref().unwrap()));
            } else {
                self.clear_paths();
                return Err(String::from("Could not find the reference path"));
            }
        }

        // Haplotype output.
        if output == HaplotypeOutput::Distinct {
            self.distinct_paths();
        } else if output == HaplotypeOutput::ReferenceOnly {
            self.reference_only()?;
        }

        Ok(())
    }

    // Sorts the paths and merges duplicates, storing the count in the weight field.
    fn distinct_paths(&mut self) {
        let ref_path = self.ref_id.map(|id| self.paths[id].path.clone());
        self.paths.sort_unstable();

        let mut new_paths: Vec<PathInfo> = Vec::new();
        let mut ref_id = None;
        for info in self.paths.iter() {
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

        self.paths = new_paths;
        self.ref_id = ref_id;
    }

    // Removes all paths except the reference path.
    fn reference_only(&mut self) -> Result<(), String> {
        if self.ref_id.is_none() {
            return Err(String::from("Reference path is required for reference-only output"));
        }
        let ref_id = self.ref_id.unwrap();
        let ref_info = self.paths[ref_id].clone();
        self.paths = vec![ref_info];
        self.ref_id = Some(0);
        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// Operations.
impl Subgraph {
    /// Returns the number of nodes in the subgraph.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.records.len() / 2
    }

    /// Returns `true` if the subgraph contains the given node.
    #[inline]
    pub fn has_node(&self, node_id: usize) -> bool {
        self.records.contains_key(&support::encode_node(node_id, Orientation::Forward))
    }

    /// Returns `true` if the subgraph contains the given handle.
    #[inline]
    pub fn has_handle(&self, handle: usize) -> bool {
        self.records.contains_key(&handle)
    }

    /// Returns the set of node identifiers in the subgraph.
    pub fn node_ids(&self) -> BTreeSet<usize> {
        self.records.keys().map(|&handle| support::node_id(handle)).collect()
    }

    /// Returns the set of handles in the subgraph.
    pub fn handles(&self) -> BTreeSet<usize> {
        self.records.keys().cloned().collect()
    }

    /// Returns the record for the given node handle, or [`None`] if the node is not in the subgraph.
    #[inline]
    pub fn record(&self, handle: usize) -> Option<&GBZRecord> {
        self.records.get(&handle)
    }

    /// Returns the record for the given node in the given orientation, or [`None`] if the node is not in the subgraph.
    #[inline]
    pub fn record_by_id(&self, node_id: usize, orientation: Orientation) -> Option<&GBZRecord> {
        self.records.get(&support::encode_node(node_id, orientation))
    }

    // Adds a new node to the subgraph.
    fn add_node_internal(&mut self, node_id: usize, get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>) -> Result<(), String> {
        let forward = get_record(support::encode_node(node_id, Orientation::Forward))?;
        let reverse = get_record(support::encode_node(node_id, Orientation::Reverse))?;
        self.records.insert(forward.handle(), forward);
        self.records.insert(reverse.handle(), reverse);
        Ok(())
    }

    /// Inserts the given node into the subgraph.
    ///
    /// No effect if the node is already in the subgraph.
    /// Clears all path information from the subgraph.
    ///
    /// # Arguments
    ///
    /// * `node_id`: Identifier of the node to insert.
    /// * `get_record`: A function that returns the record for the given node handle.
    ///
    /// # Errors
    ///
    /// Passes through any errors from `get_record`.
    pub fn add_node(&mut self, node_id: usize, get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>) -> Result<(), String> {
        if self.has_node(node_id) {
            return Ok(());
        }
        self.clear_paths();
        self.add_node_internal(node_id, get_record)
    }

    // Removes the given node from the subgraph.
    fn remove_node_internal(&mut self, node_id: usize) {
        self.records.remove(&support::encode_node(node_id, Orientation::Forward));
        self.records.remove(&support::encode_node(node_id, Orientation::Reverse));
    }

    /// Removes the given node from the subgraph.
    ///
    /// No effect if the node is not in the subgraph.
    /// Clears all path information from the subgraph.
    pub fn remove_node(&mut self, node_id: usize) {
        if !self.has_node(node_id) {
            return;
        }
        self.clear_paths();
        self.remove_node_internal(node_id);
    }

    /// Returns the number of paths in the subgraph.
    #[inline]
    pub fn path_count(&self) -> usize {
        self.paths.len()
    }

    /// Clears all path information from the subgraph.
    pub fn clear_paths(&mut self) {
        self.paths.clear();
        self.ref_id = None;
        self.ref_path = None;
        self.ref_interval = None;
    }
}

//-----------------------------------------------------------------------------

/// Alignment to the reference.
impl Subgraph {
    // Returns the total length of the nodes.
    fn path_len(&self, path: &[usize]) -> usize {
        let mut result = 0;
        for handle in path.iter() {
            let record = self.records.get(handle).unwrap();
            result += record.sequence_len();
        }
        result
    }

    // Appends a new edit or extends the previous one. No effect if `len` is zero.
    fn append_edit(edits: &mut Vec<(EditOperation, usize)>, op: EditOperation, len: usize) {
        if len == 0 {
            return;
        }
        if let Some((prev_op, prev_len)) = edits.last_mut() {
            if *prev_op == op {
                *prev_len += len;
                return;
            }
        }
        edits.push((op, len));
    }

    // Returns the number of matching bases at the start of the two paths.
    fn prefix_matches(&self, path: &[usize], ref_path: &[usize]) -> usize {
        let mut result = 0;
        let mut path_offset = 0;
        let mut ref_offset = 0;
        let mut path_base = 0;
        let mut ref_base = 0;
        while path_offset < path.len() && ref_offset < ref_path.len() {
            let a = self.records.get(&path[path_offset]).unwrap().sequence();
            let b = self.records.get(&ref_path[ref_offset]).unwrap().sequence();
            while path_base < a.len() && ref_base < b.len() {
                if a[path_base] != b[ref_base] {
                    return result;
                }
                path_base += 1;
                ref_base += 1;
                result += 1;
            }
            if path_base == a.len() {
                path_offset += 1;
                path_base = 0;
            }
            if ref_base == b.len() {
                ref_offset += 1;
                ref_base = 0;
            }
        }
        result
    }

    // Returns the number of matching bases at the end of the two paths.
    fn suffix_matches(&self, path: &[usize], ref_path: &[usize]) -> usize {
        let mut result = 0;
        let mut path_offset = 0;
        let mut ref_offset = 0;
        let mut path_base = 0;
        let mut ref_base = 0;
        while path_offset < path.len() && ref_offset < ref_path.len() {
            let a = self.records.get(&path[path.len() - path_offset - 1]).unwrap().sequence();
            let b = self.records.get(&ref_path[ref_path.len() - ref_offset - 1]).unwrap().sequence();
            while path_base < a.len() && ref_base < b.len() {
                if a[a.len() - path_base - 1] != b[b.len() - ref_base - 1] {
                    return result;
                }
                path_base += 1;
                ref_base += 1;
                result += 1;
            }
            if path_base == a.len() {
                path_offset += 1;
                path_base = 0;
            }
            if ref_base == b.len() {
                ref_offset += 1;
                ref_base = 0;
            }
        }
        result
    }

    // Returns the penalty for a mismatch of the given length.
    fn mismatch_penalty(len: usize) -> usize {
        4 * len
    }

    // Returns the penalty for a gap of the given length.
    fn gap_penalty(len: usize) -> usize {
        if len == 0 {
            0
        } else {
            6 + (len - 1)
        }
    }

    // Appends edits corresponding to the alignment of the given (sub)paths.
    // The paths are assumed to be diverging, but there may be base-level matches at the start/end.
    // The middle part is either mismatch + gap or insertion + deletion, with gap length possibly 0.
    // Alignment scoring follows the parameters used in vg:
    // +1 for a match, -4 for a mismatch, -6 for gap open, -1 for gap extend.
    //
    // NOTE: It is possible that some of the mismatches are actually matches.
    // But we ignore this possibility, as the paths are assumed to be diverging.
    fn align(&self, path: &[usize], ref_path: &[usize], edits: &mut Vec<(EditOperation, usize)>) {
        let path_len = self.path_len(path);
        let ref_len = self.path_len(ref_path);
        let prefix = self.prefix_matches(path, ref_path);
        let mut suffix = self.suffix_matches(path, ref_path);
        if prefix + suffix > path_len {
            suffix = path_len - prefix;
        }
        if prefix + suffix > ref_len {
            suffix = ref_len - prefix;
        }


        Self::append_edit(edits, EditOperation::Match, prefix);
        let path_middle = path_len - prefix - suffix;
        let ref_middle = ref_len - prefix - suffix;
        if path_middle == 0 {
            Self::append_edit(edits, EditOperation::Deletion, ref_middle);
        } else if ref_middle == 0 {
            Self::append_edit(edits, EditOperation::Insertion, path_middle);
        } else {
            let mismatch = cmp::min(path_middle, ref_middle);
            let mismatch_indel =
                Self::mismatch_penalty(mismatch) +
                Self::gap_penalty(path_middle - mismatch) +
                Self::gap_penalty(ref_middle - mismatch);
            let insertion_deletion = Self::gap_penalty(path_middle) + Self::gap_penalty(ref_middle);
            if mismatch_indel <= insertion_deletion {
                Self::append_edit(edits, EditOperation::Match, mismatch);
                Self::append_edit(edits, EditOperation::Insertion, path_middle - mismatch);
                Self::append_edit(edits, EditOperation::Deletion, ref_middle - mismatch);
            } else {
                Self::append_edit(edits, EditOperation::Insertion, path_middle);
                Self::append_edit(edits, EditOperation::Deletion, ref_middle);
            }
        }
        Self::append_edit(edits, EditOperation::Match, suffix);
    }

    // Returns the CIGAR string for the given path, aligned to the reference path.
    // Takes the alignment from the LCS of the paths weighted by node lengths.
    // Diverging parts are aligned using `align()`.
    fn align_to_ref(&self, path_id: usize) -> Option<String> {
        let ref_id = self.ref_id?;
        if path_id == ref_id || path_id >= self.paths.len() {
            return None;
        }

        // Find the LCS of the paths weighted by node lengths.
        let weight = &|handle: usize| -> usize {
            self.records.get(&handle).unwrap().sequence_len()
        };
        let path = &self.paths[path_id].path;
        let ref_path = &self.paths[ref_id].path;
        let (lcs, _) = algorithms::fast_weighted_lcs(path, ref_path, weight);

        // Convert the LCS to a sequence of edit operations
        let mut edits: Vec<(EditOperation, usize)> = Vec::new();
        let mut path_offset = 0;
        let mut ref_offset = 0;
        for (next_path_offset, next_ref_offset) in lcs.iter() {
            let path_interval = &path[path_offset..*next_path_offset];
            let ref_interval = &ref_path[ref_offset..*next_ref_offset];
            self.align(path_interval, ref_interval, &mut edits);
            let node_len = self.records.get(&path[*next_path_offset]).unwrap().sequence_len();
            Self::append_edit(&mut edits, EditOperation::Match, node_len);
            path_offset = next_path_offset + 1;
            ref_offset = next_ref_offset + 1;
        }
        self.align(&path[path_offset..], &ref_path[ref_offset..], &mut edits);

        // Convert the edits to a CIGAR string.
        let mut result = String::new();
        for (op, len) in edits.iter() {
            result.push_str(&format!("{}{}", len, op));
        }
        Some(result)
    }
}

//-----------------------------------------------------------------------------

/// Output formats.
impl Subgraph {
    /// Writes the subgraph in the GFA format to the given output.
    ///
    /// If `cigar` is true, the CIGAR strings for the non-reference haplotypes are included in the output.
    pub fn write_gfa<T: Write>(&self, output: &mut T, cigar: bool) -> io::Result<()> {
        // Header.
        let reference_samples = self.ref_path.as_ref().map(|path| path.sample.as_ref());
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
        let mut metadata = WalkMetadata::path_interval(ref_path, interval);
        metadata.add_weight(self.paths[ref_id].weight);
        Some((metadata, ref_id))
    }

    // Determines a contig name for the subgraph.
    fn contig_name(&self) -> String {
        if let Some(ref_path) = self.ref_path.as_ref() {
            ref_path.contig.clone()
        } else {
            String::from("unknown")
        }
    }
}

//-----------------------------------------------------------------------------

// Returns the graph position and the GBWT position for the given offset.
fn query_position_from_db(graph: &mut GraphInterface, path: &GBZPath, query_offset: usize) -> Result<PathPosition, String> {
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
            graph_pos = Some(GraphPosition::new(
                support::node_id(handle),
                support::node_orientation(handle),
                query_offset - path_offset
            ));
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
    Ok(PathPosition {
        seq_offset: query_offset,
        graph_pos,
        gbwt_pos,
    })
}

// Returns the graph position and the GBWT position for the given offset.
fn query_position_from_gbz(graph: &GBZ, path_index: &PathIndex, path: &GBZPath, query_offset: usize) -> Result<PathPosition, String> {
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
            graph_pos = Some(GraphPosition::new(
                node_id,
                support::node_orientation(pos.node),
                query_offset - path_offset,
            ));
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
    Ok(PathPosition {
        seq_offset: query_offset,
        graph_pos,
        gbwt_pos,
    })
}

//-----------------------------------------------------------------------------

/// A path position represented in multiple ways.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathPosition {
    /// Offset in the path in bp.
    pub seq_offset: usize,
    /// Graph position corresponding to the offset.
    pub graph_pos: GraphPosition,
    /// GBWT position for the oriented node containing the position.
    pub gbwt_pos: Pos,
}

// TODO: Add hash of the path for faster comparisons?
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

    // Returns the sequence interval (in bp) for this path, assuming that `path[path_offset]` matches `path_pos`.
    // The sequence interval is for the part of the path contained in the subgraph.
    fn path_interval(&self, subgraph: &Subgraph, path_offset: usize, path_pos: &PathPosition) -> Range<usize> {
        // Distance from the start of the path to `path_pos`.
        let mut ref_pos = path_pos.graph_pos.offset;
        for handle in self.path.iter().take(path_offset) {
            ref_pos += subgraph.records.get(handle).unwrap().sequence_len();
        }

        let start = path_pos.seq_offset - ref_pos;
        let end = start + self.len;
        start..end
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

//-----------------------------------------------------------------------------
