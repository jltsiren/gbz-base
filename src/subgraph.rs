//! A subgraph in a GBZ graph.
//!
//! This module provides functionality for extracting a subgraph around a specific position or interval of a specific path.
//! The subgraph contains all nodes within a given context and all edges between them.
//! All other paths within the subgraph can also be extracted, but they will not have any true metadata associated with them.

// TODO: Could we just provide the GBZ / DB versions of get_record() somewhere?

use crate::{GBZRecord, GBZPath, GraphInterface, PathIndex};
use crate::formats::{self, WalkMetadata, JSONValue};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap, BTreeSet};
use std::fmt::Display;
use std::io::{self, Write};
use std::iter::FusedIterator;
use std::ops::Range;
use std::cmp;

use gbwt::ENDMARKER;
use gbwt::{GBZ, GraphPosition, Orientation, Pos, FullPathName};
use gbwt::bwt::Record;

use gbwt::{algorithms, support};

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

#[derive(Clone, Debug, PartialEq, Eq)]
enum QueryType {
    // Path name and offset in bp stored in the fragment field.
    PathOffset(FullPathName),
    // Starting position as in `PathOffset` and length in bp.
    PathInterval((FullPathName, usize)),
    // Set of node identifiers.
    Nodes(BTreeSet<usize>),
}

// but we first need to implement context extraction based on an interval.
/// Arguments for extracting a subgraph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubgraphQuery {
    query_type: QueryType,

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
        let mut path_name = path_name.clone();
        path_name.fragment = offset;
        SubgraphQuery {
            query_type: QueryType::PathOffset(path_name),
            context,
            output,
        }
    }

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
        let mut path_name = path_name.clone();
        path_name.fragment = interval.start;
        SubgraphQuery {
            query_type: QueryType::PathInterval((path_name, interval.len())),
            context,
            output,
        }
    }

    /// Creates a query that retrieves a subgraph around a node.
    ///
    /// # Arguments
    ///
    /// * `nodes`: Set of node identifiers.
    /// * `context`: Context length around the reference node (in bp).
    /// * `output`: How to output the haplotypes.
    ///
    /// # Panics
    ///
    /// Panics if `output` is [`HaplotypeOutput::ReferenceOnly`].
    pub fn nodes(nodes: impl IntoIterator<Item = usize>, context: usize, output: HaplotypeOutput) -> Self {
        if output == HaplotypeOutput::ReferenceOnly {
            panic!("Cannot output a reference path in a node-based query");
        }
        SubgraphQuery {
            query_type: QueryType::Nodes(nodes.into_iter().collect()),
            context,
            output,
        }
    }

    fn query_type(&self) -> &QueryType {
        &self.query_type
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
        match self.query_type() {
            QueryType::PathOffset(path_name) => write!(f, "(path {}, context {}, {})", path_name, self.context, self.output),
            QueryType::PathInterval((path_name, len)) => write!(f, "(path {}, len {}, context {}, {})", path_name, len, self.context, self.output),
            QueryType::Nodes(nodes) => write!(f, "(nodes {:#?}, context {}, {})", nodes, self.context, self.output),
        }
    }
}

//-----------------------------------------------------------------------------

/// A subgraph induced by a set of node identifiers.
///
/// The subgraph contains the specified nodes and the edges between them.
/// Paths in the current subgraph can be extracted using [`Subgraph::extract_paths`].
/// Only the provided reference path will have proper metadata.
/// Other paths remain anonymous, as we cannot identify them efficiently using the GBWT.
/// Any changes to the subgraph will remove the extracted paths.
/// When a subgraph is changed, the operations will try to reuse already extracted node records.
/// This makes it viable as a sliding window over a reference path.
///
/// Node identifiers can be selected using:
/// * [`Subgraph::add_node`] and [`Subgraph::remove_node`] with individual ids.
/// * [`Subgraph::around_position`] for a context around a graph position.
/// * [`Subgraph::around_interval`] for a context around path interval.
/// * [`Subgraph::around_nodes`] for a context around a set of nodes.
///
/// [`Subgraph::from_gbz`] and [`Subgraph::from_db`] are integrated methods for extracting a subgraph with paths using a [`SubgraphQuery`].
///
/// `Subgraph` implements a similar graph interface to the node/edge operations of [`GBZ`].
/// It can also be serialized in GFA and JSON formats using [`Subgraph::write_gfa`] and [`Subgraph::write_json`].
///
/// # Examples
///
/// The following example replicates an offset-based [`Subgraph::from_gbz`] query using lower-level functions.
///
/// ```
/// use gbz_base::{GBZRecord, PathIndex, Subgraph, HaplotypeOutput};
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
/// // Find the position for path A offset 2.
/// let mut query_pos = FullPathName::generic("A");
/// query_pos.fragment = 2;
/// let mut subgraph = Subgraph::new();
/// let result = subgraph.path_pos_from_gbz(&graph, &path_index, &query_pos);
/// assert!(result.is_ok());
/// let (path_pos, path_name) = result.unwrap();
///
/// // Extract a 1 bp context around the position.
/// let result = subgraph.around_position(path_pos.graph_pos(), 1, &mut |handle| {
///    GBZRecord::from_gbz(&graph, handle).ok_or(
///        format!("The graph does not contain handle {}", handle)
///    )
/// });
/// assert!(result.is_ok());
///
/// // Extract all paths in the subgraph.
/// let result = subgraph.extract_paths(Some((path_pos, path_name)), HaplotypeOutput::All);
/// assert!(result.is_ok());
///
/// // The subgraph should be centered around 1 bp node 14 of degree 4.
/// assert_eq!(subgraph.nodes(), 5);
/// assert_eq!(subgraph.paths(), 3);
/// ```
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

/// Construction.
impl Subgraph {
    /// Creates a new empty subgraph.
    pub fn new() -> Self {
        Subgraph::default()
    }

    /// Returns the path position for the haplotype offset represented by the query position.
    ///
    /// `query_pos.fragment` is used as an offset in the haplotype.
    /// The return value consists of the position and metadata for the path covering the position.
    /// Updates the subgraph to include the path from the nearest indexed position to the query position.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ graph.
    /// * `path_index`: A path index for the graph.
    /// * `query_pos`: Query position.
    ///
    /// # Errors
    ///
    /// Returns an error if there is no path covering the given position or the path has not been indexed for random access.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZRecord, Subgraph, PathIndex};
    /// use gbwt::{GBZ, FullPathName, Orientation};
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
    /// // Query for path A offset 2.
    /// let mut query_pos = FullPathName::generic("A");
    /// query_pos.fragment = 2;
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.path_pos_from_gbz(&graph, &path_index, &query_pos);
    /// assert!(result.is_ok());
    /// let (path_pos, path_name) = result.unwrap();
    /// assert_eq!(path_pos.seq_offset(), query_pos.fragment);
    ///
    /// // The query position is at the start of node 14 in forward orientation.
    /// assert_eq!(path_pos.node_id(), 14);
    /// assert_eq!(path_pos.orientation(), Orientation::Forward);
    /// assert_eq!(path_pos.node_offset(), 0);
    ///
    /// // And this happens to be the first path visit to the node.
    /// let gbwt_node = support::encode_node(14, Orientation::Forward);
    /// assert_eq!(path_pos.handle(), gbwt_node);
    /// assert_eq!(path_pos.gbwt_offset(), 0);
    ///
    /// // And it is covered by the path fragment starting from offset 0.
    /// assert_eq!(path_name, FullPathName::generic("A"));
    ///
    /// // Now extract 1 bp context around interval 2..4.
    /// let result = subgraph.around_interval(path_pos, 2, 1, &mut |handle| {
    ///     GBZRecord::from_gbz(&graph, handle).ok_or(
    ///         format!("The graph does not contain handle {}", handle)
    ///     )
    /// });
    /// assert!(result.is_ok());
    ///
    /// // The interval corresponds to nodes 14 and 15.
    /// let true_nodes = vec![12, 13, 14, 15, 16, 17];
    /// assert_eq!(subgraph.nodes(), true_nodes.len());
    /// assert!(subgraph.node_iter().eq(true_nodes.iter().copied()));
    /// ```
    pub fn path_pos_from_gbz(
        &mut self,
        graph: &GBZ,
        path_index: &PathIndex,
        query_pos: &FullPathName
    ) -> Result<(PathPosition, FullPathName), String> {
        let path = GBZPath::with_name(graph, query_pos).ok_or(
            format!("Cannot find a path covering {}", query_pos)
        )?;
        // Transform the offset relative to the haplotype to the offset relative to the path.
        let query_offset = query_pos.fragment - path.name.fragment;

        // Path id to an indexed position.
        let index_offset = path_index.path_to_offset(path.handle).ok_or(
            format!("Path {} has not been indexed for random access", path.name())
        )?;
        let (path_offset, pos) = path_index.indexed_position(index_offset, query_offset).unwrap();

        self.find_path_pos(query_offset, path_offset, pos, path.name, &mut |handle| {
            GBZRecord::from_gbz(graph, handle).ok_or(
                format!("The graph does not contain handle {}", handle)
            )
        })
    }

    /// Returns the path position for the haplotype offset represented by the query position.
    ///
    /// `query_pos.fragment` is used as an offset in the haplotype.
    /// The return value consists of the position and metadata for the path covering the position.
    /// Updates the subgraph to include the path from the nearest indexed position to the query position.
    ///
    /// # Arguments
    ///
    /// * `graph`: A graph interface.
    /// * `query_pos`: Query position.
    ///
    /// # Errors
    ///
    /// Returns an error if database operations fail.
    /// Returns an error if there is no path covering the given position or the path has not been indexed for random access.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZBase, GraphInterface, Subgraph};
    /// use gbwt::{FullPathName, Orientation};
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
    /// // Query for path A offset 2.
    /// let mut query_pos = FullPathName::generic("A");
    /// query_pos.fragment = 2;
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.path_pos_from_db(&mut interface, &query_pos);
    /// assert!(result.is_ok());
    /// let (path_pos, path_name) = result.unwrap();
    /// assert_eq!(path_pos.seq_offset(), query_pos.fragment);
    ///
    /// // The query position is at the start of node 14 in forward orientation.
    /// assert_eq!(path_pos.node_id(), 14);
    /// assert_eq!(path_pos.orientation(), Orientation::Forward);
    /// assert_eq!(path_pos.node_offset(), 0);
    ///
    /// // And this happens to be the first path visit to the node.
    /// let gbwt_node = support::encode_node(14, Orientation::Forward);
    /// assert_eq!(path_pos.handle(), gbwt_node);
    /// assert_eq!(path_pos.gbwt_offset(), 0);
    ///
    /// // And it is covered by the path fragment starting from offset 0.
    /// assert_eq!(path_name, FullPathName::generic("A"));
    ///
    /// // Clean up.
    /// drop(interface);
    /// drop(database);
    /// fs::remove_file(&db_file).unwrap();
    /// ```
    pub fn path_pos_from_db(
        &mut self,
        graph: &mut GraphInterface,
        query_pos: &FullPathName
    ) -> Result<(PathPosition, FullPathName), String> {
        let path = graph.find_path(query_pos)?;
        let path = path.ok_or(format!("Cannot find a path covering {}", query_pos))?;
        if !path.is_indexed {
            return Err(format!("Path {} has not been indexed for random access", query_pos));
        }
        // Transform the offset relative to the haplotype to the offset relative to the path.
        let query_offset = query_pos.fragment - path.name.fragment;

        // Find an indexed position before the query position.
        let result = graph.indexed_position(path.handle, query_offset)?;
        let (path_offset, pos) = result.ok_or(
            format!("Path {} has not been indexed for random access", path.name())
        )?;

        self.find_path_pos(query_offset, path_offset, pos, path.name, &mut |handle| {
            let record = graph.get_record(handle)?;
            record.ok_or(format!("The graph does not contain handle {}", handle))
        })
    }

    // Shared functionality for `path_pos_from_gbz` and `path_pos_from_db`.
    fn find_path_pos(
        &mut self,
        query_offset: usize,
        path_offset: usize,
        pos: Pos,
        path_name: FullPathName,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(PathPosition, FullPathName), String> {
        // Iterate over the path until the query position, updating the subgraph.
        let mut path_offset = path_offset;
        let mut pos = pos;
        let mut node_offset: Option<usize> = None;
        let mut gbwt_pos: Option<Pos> = None;
        loop {
            let node_id = support::node_id(pos.node);
            if !self.has_node(node_id) {
                self.add_node_internal(node_id, get_record)?;
            }
            let record = self.record(pos.node).unwrap();
            if path_offset + record.sequence_len() > query_offset {
                node_offset = Some(query_offset - path_offset);
                gbwt_pos = Some(pos);
                break;
            }
            path_offset += record.sequence_len();
            let next = record.to_gbwt_record().lf(pos.offset);
            if next.is_none() {
                break;
            }
            pos = next.unwrap();
        }

        let node_offset = node_offset.ok_or(
            format!("Path {} does not contain offset {}", path_name, query_offset)
        )?;
        let gbwt_pos = gbwt_pos.unwrap();

        let path_position = PathPosition {
            seq_offset: query_offset,
            gbwt_node: gbwt_pos.node,
            node_offset,
            gbwt_offset: gbwt_pos.offset,
        };
        Ok((path_position, path_name))
    }

    /// Updates the subgraph to a context around the given graph position.
    ///
    /// Reuses existing records when possible.
    /// Removes node records outside the context as well as all path information.
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
    /// use gbz_base::{GBZRecord, Subgraph};
    /// use gbwt::{GBZ, GraphPosition, Orientation};
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
    /// let result = subgraph.around_position(pos, 1, &mut |handle| {
    ///     GBZRecord::from_gbz(&graph, handle).ok_or(
    ///         format!("The graph does not contain handle {}", handle)
    ///     )
    /// });
    /// assert!(result.is_ok());
    ///
    /// // The subgraph should be centered around 1 bp node 14 of degree 4.
    /// let true_nodes = [12, 13, 14, 15, 16];
    /// assert_eq!(subgraph.nodes(), true_nodes.len());
    /// assert!(subgraph.node_iter().eq(true_nodes.iter().copied()));
    ///
    /// // But there are no paths.
    /// assert_eq!(subgraph.paths(), 0);
    /// ```
    pub fn around_position(
        &mut self,
        pos: GraphPosition,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        // The initial node is always in the subgraph, so we might as well add it now to determine sequence length.
        if !self.has_node(pos.node) {
            self.add_node_internal(pos.node, get_record)?;
        }
        let handle = pos.to_gbwt();
        let record = self.record(handle).unwrap();

        // Start the graph traversal from both sides of the initial node.
        let mut active: BinaryHeap<Reverse<(usize, (usize, NodeSide))>> = BinaryHeap::new();
        let start = (pos.node, NodeSide::start(pos.orientation));
        active.push(Reverse((pos.offset, start)));
        let end_distance = record.sequence_len() - pos.offset - 1;
        let end = (pos.node, NodeSide::end(pos.orientation));
        active.push(Reverse((end_distance, end)));

        self.insert_context(active, context, get_record)
    }

    /// Updates the subgraph to a context around the given path interval.
    ///
    /// Reuses existing records when possible.
    /// Removes node records outside the context as well as all path information.
    /// See [`Self::path_pos_from_gbz`] for an example.
    ///
    /// # Arguments
    ///
    /// * `start_pos`: The starting position for the interval.
    /// * `len`: Length of the interval (in bp).
    /// * `context`: The context length around the reference interval (in bp).
    /// * `get_record`: A function that returns the record for the given node handle.
    ///
    /// # Errors
    ///
    /// Passes through any errors from `get_record`.
    /// Returns an error if the interval is empty, starts outside the initial node, or is longer than the remaining path.
    pub fn around_interval(
        &mut self,
        start_pos: PathPosition,
        len: usize,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        if len == 0 {
            return Err(String::from("Interval length must be greater than 0"));
        }
        let mut pos = start_pos.gbwt_pos();
        let mut offset = start_pos.node_offset();
        let mut len = len;

        // Insert all nodes in the interval to the subgraph and to active node sides.
        let mut active: BinaryHeap<Reverse<(usize, (usize, NodeSide))>> = BinaryHeap::new();
        loop {
            let node_id = support::node_id(pos.node);
            let orientation = support::node_orientation(pos.node);
            if !self.has_node(node_id) {
                self.add_node_internal(node_id, get_record)?;
            }
            let record = self.record(pos.node).unwrap();
            if offset >= record.sequence_len() {
                return Err(format!("Offset {} in node {} of length {}", offset, node_id, record.sequence_len()));
            }

            // Handle the current node.
            let start = (node_id, NodeSide::start(orientation));
            active.push(Reverse((offset, start)));
            let end = (node_id, NodeSide::end(orientation));
            let distance_to_next = record.sequence_len() - offset;
            if len <= distance_to_next {
                let end_distance = if len == distance_to_next { 0 } else { distance_to_next - len - 1 };
                active.push(Reverse((end_distance, end)));
                break;
            } else {
                active.push(Reverse((0, end)));
            }

            // Proceed to the next node.
            if let Some(next) = record.to_gbwt_record().lf(pos.offset) {
                pos = next;
            } else {
                return Err(format!("No successor for GBWT position ({}, {})", pos.node, pos.offset));
            }
            offset = 0;
            len -= distance_to_next;
        }

        self.insert_context(active, context, get_record)
    }

    /// Updates the subgraph to a context around the given nodes.
    ///
    /// Reuses existing records when possible.
    /// Removes node records outside the context as well as all path information.
    ///
    /// # Arguments
    ///
    /// * `nodes`: Set of reference nodes for the subgraph.
    /// * `context`: The context length around the reference node (in bp).
    /// * `get_record`: A function that returns the record for the given node handle.
    ///
    /// # Errors
    ///
    /// Passes through any errors from `get_record`.
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZRecord, Subgraph};
    /// use gbwt::GBZ;
    /// use gbwt::support;
    /// use simple_sds::serialize;
    /// use std::collections::BTreeSet;
    ///
    /// // Get the graph.
    /// let gbz_file = support::get_test_data("translation.gbz");
    /// let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    ///
    /// // Extract a subgraph that contains an 1 bp context around node 5.
    /// let mut subgraph = Subgraph::new();
    /// let mut nodes = BTreeSet::new();
    /// nodes.insert(5);
    /// let result = subgraph.around_nodes(&nodes, 1, &mut |handle| {
    ///     GBZRecord::from_gbz(&graph, handle).ok_or(
    ///         format!("The graph does not contain handle {}", handle)
    ///     )
    /// });
    /// assert!(result.is_ok());
    ///
    /// // The subgraph should be centered around 2 bp node 5 of degree 3.
    /// let true_nodes = [3, 4, 5, 6];
    /// assert_eq!(subgraph.nodes(), true_nodes.len());
    /// assert!(subgraph.node_iter().eq(true_nodes.iter().copied()));
    ///
    /// // But there are no paths.
    /// assert_eq!(subgraph.paths(), 0);
    /// ```
    pub fn around_nodes(
        &mut self,
        nodes: &BTreeSet<usize>,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        // Start the graph traversal from both sides of the initial nodes.
        let mut active: BinaryHeap<Reverse<(usize, (usize, NodeSide))>> = BinaryHeap::new();
        for &node_id in nodes {
            let start = (node_id, NodeSide::Left);
            active.push(Reverse((0, start)));
            let end = (node_id, NodeSide::Right);
            active.push(Reverse((0, end)));
        }

        self.insert_context(active, context, get_record)
    }

    // FIXME: check that if we call this twice, we do not add or remove any records.
    // Inserts all nodes within the context, starting from the active node sides.
    //
    // All nodes in `active` are assumed to be in the subgraph, even if the distances to the sides exceed the context.
    // Reuses existing records if possible.
    // Removes existing nodes that are not in the new context.
    // Clears all path information.
    fn insert_context(
        &mut self,
        active: BinaryHeap<Reverse<(usize, (usize, NodeSide))>>,
        context: usize,
        get_record: &mut dyn FnMut(usize) -> Result<GBZRecord, String>
    ) -> Result<(), String> {
        self.clear_paths();

        let mut active = active;
        let mut visited: BTreeSet<(usize, NodeSide)> = BTreeSet::new();
        let mut to_remove: BTreeSet<usize> = self.node_iter().collect();

        while !active.is_empty() {
            let (distance, node_side) = active.pop().unwrap().0;
            if visited.contains(&node_side) {
                continue;
            }
            visited.insert(node_side);
            to_remove.remove(&node_side.0);
            if !self.has_node(node_side.0) {
                self.add_node_internal(node_side.0, get_record)?;
            }

            // We can reach the other side by traversing the node.
            let other_side = (node_side.0, node_side.1.flip());
            if !visited.contains(&other_side) {
                let handle = support::encode_node(node_side.0, node_side.1.as_start());
                let record = self.record(handle).unwrap();
                let next_distance = distance + record.sequence_len() - 1;
                if next_distance <= context {
                    active.push(Reverse((next_distance, other_side)));
                }
            }

            // The predecessors of this node side are 1 bp away.
            let handle = support::encode_node(node_side.0, node_side.1.as_end());
            let record = self.record(handle).unwrap();
            let next_distance = distance + 1;
            if next_distance <= context {
                for successor in record.successors() {
                    let node_id = support::node_id(successor);
                    let side = NodeSide::start(support::node_orientation(successor));
                    if !visited.contains(&(node_id, side)) {
                        active.push(Reverse((next_distance, (node_id, side))));
                    }
                }
            }
        }

        for node_id in to_remove {
            self.remove_node_internal(node_id);
        }
        Ok(())
    }

    /// Extracts a subgraph around the given query position.
    ///
    /// Reuses existing records when possible.
    /// Removes node records not covered by the query.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ graph.
    /// * `path_index`: A path index for the graph, if the query is path-based.
    /// * `query`: Arguments for extracting the subgraph.
    ///
    /// # Errors
    ///
    /// Returns an error if the query or the graph is invalid.
    /// Returns an error if the graph does not contain the queried position.
    /// Returns an error if a path index is required but not provided, or if the query path has not been indexed.
    /// If an error occurs, the subgraph may contain arbitrary nodes but no paths.
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
    /// assert_eq!(subgraph.nodes(), 5);
    /// assert_eq!(subgraph.paths(), 3);
    ///
    /// // We get the same result using a node id.
    /// let query = SubgraphQuery::nodes([14], 1, HaplotypeOutput::All);
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.from_gbz(&graph, None, &query);
    /// assert!(result.is_ok());
    /// assert_eq!(subgraph.nodes(), 5);
    /// assert_eq!(subgraph.paths(), 3);
    /// ```
    pub fn from_gbz(&mut self, graph: &GBZ, path_index: Option<&PathIndex>, query: &SubgraphQuery) -> Result<(), String> {
        match query.query_type() {
            QueryType::PathOffset(query_pos) => {
                let path_index = path_index.ok_or(
                    String::from("Path index is required for path-based queries")
                )?;
                let reference_path = self.path_pos_from_gbz(graph, path_index, query_pos)?;
                self.around_position(reference_path.0.graph_pos(), query.context, &mut |handle| {
                    GBZRecord::from_gbz(graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(Some(reference_path), query.output())?;
            }
            QueryType::PathInterval((query_pos, len)) => {
                let path_index = path_index.ok_or(
                    String::from("Path index is required for path-based queries")
                )?;
                let reference_path = self.path_pos_from_gbz(graph, path_index, query_pos)?;
                self.around_interval(reference_path.0, *len, query.context, &mut |handle| {
                    GBZRecord::from_gbz(graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(Some(reference_path), query.output())?;
            }
            QueryType::Nodes(nodes) => {
                if query.output() == HaplotypeOutput::ReferenceOnly {
                    return Err(String::from("Cannot output a reference path in a node-based query"));
                }
                self.around_nodes(nodes, query.context, &mut |handle| {
                    GBZRecord::from_gbz(graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(None, query.output())?;
            }
        }

        Ok(())
    }

    /// Extracts a subgraph around the given query position.
    ///
    /// Reuses existing records when possible.
    /// Removes node records not covered by the query.
    ///
    /// # Errors
    ///
    /// Returns an error if the query or the graph is invalid or if there is a database error.
    /// Returns an error if the graph does not contain the queried position.
    /// Returns an error if the query path has not been indexed.
    /// If an error occurs, the subgraph may contain arbitrary nodes but no paths.
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
    /// assert_eq!(subgraph.nodes(), 5);
    /// assert_eq!(subgraph.paths(), 3);
    ///
    /// // Clean up.
    /// drop(interface);
    /// drop(database);
    /// fs::remove_file(&db_file).unwrap();
    /// ```
    pub fn from_db(&mut self, graph: &mut GraphInterface, query: &SubgraphQuery) -> Result<(), String> {
        match query.query_type() {
            QueryType::PathOffset(query_pos) => {
                let reference_path = self.path_pos_from_db(graph, query_pos)?;
                self.around_position(reference_path.0.graph_pos(), query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(Some(reference_path), query.output())?;
            }
            QueryType::PathInterval((query_pos, len)) => {
                let reference_path = self.path_pos_from_db(graph, query_pos)?;
                self.around_interval(reference_path.0, *len, query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(Some(reference_path), query.output())?;
            }
            QueryType::Nodes(nodes) => {
                if query.output() == HaplotypeOutput::ReferenceOnly {
                    return Err(String::from("Cannot output a reference path in a node-based query"));
                }
                self.around_nodes(nodes, query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                })?;
                self.extract_paths(None, query.output())?;
            }
        }

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
    ///
    /// # Examples
    ///
    /// ```
    /// use gbz_base::{GBZRecord, Subgraph, HaplotypeOutput};
    /// use gbwt::GBZ;
    /// use gbwt::support;
    /// use simple_sds::serialize;
    ///
    /// // Get the graph.
    /// let gbz_file = support::get_test_data("example.gbz");
    /// let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    /// let get_record = &mut |handle| {
    ///     GBZRecord::from_gbz(&graph, handle).ok_or(
    ///         format!("The graph does not contain handle {}", handle)
    ///     )
    /// };
    ///
    /// // Start with an empty subgraph and add a node.
    /// let mut subgraph = Subgraph::new();
    /// let result = subgraph.add_node(15, get_record);
    /// assert!(result.is_ok());
    ///
    /// // There are no paths until we extract them.
    /// assert_eq!(subgraph.paths(), 0);
    /// let result = subgraph.extract_paths(None, HaplotypeOutput::All);
    /// assert!(result.is_ok());
    /// assert_eq!(subgraph.paths(), 2);
    ///
    /// // When we change the subgraph, we need to extract the paths again.
    /// for node_id in [13, 14] {
    ///    let result = subgraph.add_node(node_id, get_record);
    ///    assert!(result.is_ok());
    /// }
    /// assert_eq!(subgraph.paths(), 0);
    /// let result = subgraph.extract_paths(None, HaplotypeOutput::All);
    /// assert!(result.is_ok());
    /// assert_eq!(subgraph.paths(), 3);
    ///
    /// // The same is true for removing nodes.
    /// for node_id in [14, 15] {
    ///     subgraph.remove_node(node_id);
    /// }
    /// assert_eq!(subgraph.paths(), 0);
    /// let result = subgraph.extract_paths(None, HaplotypeOutput::All);
    /// assert!(result.is_ok());
    /// assert_eq!(subgraph.paths(), 1);
    /// ```
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
                        if pos == position.gbwt_pos() {
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

    // Clears all path information from the subgraph.
    fn clear_paths(&mut self) {
        self.paths.clear();
        self.ref_id = None;
        self.ref_path = None;
        self.ref_interval = None;
    }
}

//-----------------------------------------------------------------------------

// FIXME Examples using these
/// Node/edge operations similar to [`GBZ`].
impl Subgraph {
    /// Returns the number of nodes in the subgraph.
    #[inline]
    pub fn nodes(&self) -> usize {
        self.records.len() / 2
    }

    /// Returns the smallest node identifier in the subgraph.
    #[inline]
    pub fn min_node(&self) -> Option<usize> {
        self.records.keys().next().map(|&handle| support::node_id(handle))
    }

    /// Returns the largest node identifier in the subgraph.
    #[inline]
    pub fn max_node(&self) -> Option<usize> {
        self.records.keys().next_back().map(|&handle| support::node_id(handle))
    }

    /// Returns `true` if the subgraph contains the given node.
    #[inline]
    pub fn has_node(&self, node_id: usize) -> bool {
        self.records.contains_key(&support::encode_node(node_id, Orientation::Forward))
    }

    // Returns `true` if the subgraph contains the given handle.
    #[inline]
    pub fn has_handle(&self, handle: usize) -> bool {
        self.records.contains_key(&handle)
    }

    /// Returns the sequence for the node in the subgraph, or [`None`] if there is no such node.
    #[inline]
    pub fn sequence(&self, node_id: usize) -> Option<&[u8]> {
        self.records.get(&support::encode_node(node_id, Orientation::Forward)).map(|record| record.sequence())
    }

    /// Returns the sequence for the node in the given orientation, or [`None`] if there is no such node.
    #[inline]
    pub fn oriented_sequence(&self, node_id: usize, orientation: Orientation) -> Option<&[u8]> {
        self.records.get(&support::encode_node(node_id, orientation)).map(|record| record.sequence())
    }

    /// Returns the length of the sequence for the node in the subgraph, or [`None`] if there is no such node.
    #[inline]
    pub fn sequence_len(&self, node_id: usize) -> Option<usize> {
        self.records.get(&support::encode_node(node_id, Orientation::Forward)).map(|record| record.sequence_len())
    }

    /// Returns an iterator over the node identifiers in the subgraph.
    ///
    /// The identifiers are listed in ascending order.
    pub fn node_iter(&'_ self) -> impl Iterator<Item = usize> + use<'_> {
        self.records.keys().step_by(2).map(|&handle| support::node_id(handle))
    }

    /// Returns an iterator over the successors of an oriented node, or [`None`] if there is no such node.
    pub fn successors(&self, node_id: usize, orientation: Orientation) -> Option<EdgeIter> {
        let handle = support::encode_node(node_id, orientation);
        let gbz_record = self.records.get(&handle)?;
        let record = gbz_record.to_gbwt_record();
        Some(EdgeIter::new(record, self, false))
    }

    /// Returns an iterator over the predecessors of an oriented node, or [`None`] if there is no such node.
    pub fn predecessors(&self, node_id: usize, orientation: Orientation) -> Option<EdgeIter> {
        let handle = support::encode_node(node_id, orientation.flip());
        let gbz_record = self.records.get(&handle)?;
        let record = gbz_record.to_gbwt_record();
        Some(EdgeIter::new(record, self, true))
    }

    // Returns the record for the given node handle, or [`None`] if the node is not in the subgraph.
    #[inline]
    fn record(&self, handle: usize) -> Option<&GBZRecord> {
        self.records.get(&handle)
    }

    /// Returns the number of paths in the subgraph.
    #[inline]
    pub fn paths(&self) -> usize {
        self.paths.len()
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
                if self.has_handle(successor) && support::edge_is_canonical(from, to) {
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
                if self.has_handle(successor) && support::edge_is_canonical(from, to) {
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

/// A path position represented in multiple ways.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PathPosition {
    // Sequence offset from the start of the path.
    seq_offset: usize,
    // GBWT node identifier encoding an oriented node.
    gbwt_node: usize,
    // Offset from the start of the node.
    node_offset: usize,
    // Offset within GBWT record.
    gbwt_offset: usize,
}

impl PathPosition {
    /// Returns the sequence offset in bp from the start of the path.
    #[inline]
    pub fn seq_offset(&self) -> usize {
        self.seq_offset
    }

    /// Returns the node identifier.
    #[inline]
    pub fn node_id(&self) -> usize {
        support::node_id(self.gbwt_node)
    }

    /// Returns the orientation of the node.
    #[inline]
    pub fn orientation(&self) -> Orientation {
        support::node_orientation(self.gbwt_node)
    }

    /// Returns the offset from the start of the node.
    #[inline]
    pub fn node_offset(&self) -> usize {
        self.node_offset
    }

    /// Returns the graph position.
    #[inline]
    pub fn graph_pos(&self) -> GraphPosition {
        GraphPosition {
            node: self.node_id(),
            orientation: self.orientation(),
            offset: self.node_offset(),
        }
    }

    /// Returns the GBWT node identifier / handle.
    #[inline]
    pub fn handle(&self) -> usize {
        self.gbwt_node
    }

    /// Returns the offset within the GBWT record.
    #[inline]
    pub fn gbwt_offset(&self) -> usize {
        self.gbwt_offset
    }

    /// Returns the GBWT position for the oriented node containing the position.
    #[inline]
    pub fn gbwt_pos(&self) -> Pos {
        Pos {
            node: self.handle(),
            offset: self.gbwt_offset(),
        }
    }
}

// TODO: This might belong to gbwt-rs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum NodeSide {
    Left,
    Right,
}

impl NodeSide {
    fn flip(&self) -> NodeSide {
        match self {
            NodeSide::Left => NodeSide::Right,
            NodeSide::Right => NodeSide::Left,
        }
    }

    fn as_start(&self) -> Orientation {
        match self {
            NodeSide::Left => Orientation::Forward,
            NodeSide::Right => Orientation::Reverse,
        }
    }

    fn as_end(&self) -> Orientation {
        match self {
            NodeSide::Left => Orientation::Reverse,
            NodeSide::Right => Orientation::Forward,
        }
    }

    fn start(orientation: Orientation) -> NodeSide {
        if orientation == Orientation::Forward {
            NodeSide::Left
        } else {
            NodeSide::Right
        }
    }

    fn end(orientation: Orientation) -> NodeSide {
        if orientation == Orientation::Forward {
            NodeSide::Right
        } else {
            NodeSide::Left
        }
    }
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
        let mut ref_pos = path_pos.node_offset();
        for handle in self.path.iter().take(path_offset) {
            ref_pos += subgraph.records.get(handle).unwrap().sequence_len();
        }

        let start = path_pos.seq_offset() - ref_pos;
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

// FIXME examples
/// An iterator over the predecessors or successors of a node.
///
/// The type of `Item` is `(`[`usize`]`, `[`Orientation`]`)`.
/// These values encode the identifier and the orientation of the node.
/// Successor nodes are always listed in sorted order.
/// Predecessor nodes are sorted by their identifiers, but the reverse orientation is listed before the forward orientation.
/// Like [`gbwt::gbz::EdgeIter`], but only for edges in the subgraph.
#[derive(Clone, Debug)]
pub struct EdgeIter<'a> {
    record: Record<'a>, // FIXME: We could simply have an edge slice.
    parent: &'a Subgraph,
    // The first outrank that has not been visited.
    next: usize,
    // The first outrank that we should not visit.
    limit: usize,
    // Flip the orientation in the iterated values.
    flip: bool,
}


impl<'a> EdgeIter<'a> {
    /// Creates a new iterator over the successors of the record.
    ///
    /// If `flip` is `true`, the iterator will flip the orientation of the successors.
    /// This is effectively the same as listing the predecessors of the reverse orientation of the record.
    pub fn new(record: Record<'a>, parent: &'a Subgraph, flip: bool) -> Self {
        let next = if record.outdegree() > 0 && record.successor(0) == ENDMARKER { 1 } else { 0 };
        let limit = record.outdegree();
        EdgeIter {
            record,
            parent,
            next,
            limit,
            flip,
        }
    }
}

impl<'a> Iterator for EdgeIter<'a> {
    type Item = (usize, Orientation);

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.limit {
            let handle = self.record.successor(self.next);
            self.next += 1;
            if self.parent.has_handle(handle) {
                let node_id = support::node_id(handle);
                let orientation = if self.flip {
                    support::node_orientation(handle).flip()
                } else {
                    support::node_orientation(handle)
                };
                return Some((node_id, orientation));
            }
        }
        None
    }
}

impl<'a> DoubleEndedIterator for EdgeIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        while self.next < self.limit {
            let handle = self.record.successor(self.limit - 1);
            self.limit -= 1;
            if self.parent.has_handle(handle) {
                let node_id = support::node_id(handle);
                let orientation = if self.flip {
                    support::node_orientation(handle).flip()
                } else {
                    support::node_orientation(handle)
                };
                return Some((node_id, orientation));
            }
        }
        None
    }
}

impl<'a> FusedIterator for EdgeIter<'a> {}

//-----------------------------------------------------------------------------
