//! Queries for extracting a subgraph from GBZ-base or a GBZ graph.

use gbz::{support, FullPathName};

use std::collections::BTreeSet;
use std::fmt::Display;
use std::ops::Range;

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
    /// No haplotypes in the output.
    None,
}

impl Display for HaplotypeOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            HaplotypeOutput::All => write!(f, "all"),
            HaplotypeOutput::Distinct => write!(f, "distinct"),
            HaplotypeOutput::ReferenceOnly => write!(f, "reference only"),
            HaplotypeOutput::None => write!(f, "none"),
        }
    }
}

/// Output options for extending the subgraph with snarls overlapping with it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SnarlOutput {
    /// Do not extract snarls.
    None,
    /// Extract snarls with both boundary nodes contained in the subgraph.
    Contained,
    /// Extract all snarls that overlap with the subgraph.
    Overlapping,
}

impl Display for SnarlOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SnarlOutput::None => write!(f, "none"),
            SnarlOutput::Contained => write!(f, "contained"),
            SnarlOutput::Overlapping => write!(f, "overlapping"),
        }
    }
}

/// Distance calculation heuristic for context extraction.
///
/// Distances are always measured to the sides of a node.
/// The modes differ in how a node with two reachable sides is handled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DistanceMode {
    /// Distances to the two sides of a node are determined independently (default).
    #[default]
    Side,
    /// Visiting a node fixes the distance to its other side to `distance + sequence_len - 1`,
    /// even if there is a shorter route to that side.
    ///
    /// This matches the heuristic used by `gbwtgraph::find_subgraph` in the C++ implementation.
    /// The resulting subgraph is always a subset of the one produced by [`DistanceMode::Side`].
    Node,
}

impl Display for DistanceMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DistanceMode::Side => write!(f, "side"),
            DistanceMode::Node => write!(f, "node"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum QueryType {
    // Path name and offset in bp stored in the fragment field.
    PathOffset(FullPathName),
    // Starting position as in `PathOffset` and length in bp.
    PathInterval(FullPathName, usize),
    // Set of node identifiers.
    Nodes(BTreeSet<usize>),
    // Subgraph between two handles in the same chain.
    Between(usize, usize),
}

//-----------------------------------------------------------------------------

/// Arguments for extracting a subgraph.
///
/// # Examples
///
/// ```
/// use gbz_base::{SubgraphQuery, DistanceMode, SnarlOutput};
/// use gbz::FullPathName;
///
/// let path_name = FullPathName::generic("path");
/// let query = SubgraphQuery::path_offset(&path_name, 123);
/// assert_eq!(query.context(), SubgraphQuery::DEFAULT_CONTEXT);
/// assert_eq!(query.snarls(), SubgraphQuery::DEFAULT_SNARLS);
/// assert_eq!(query.output(), SubgraphQuery::DEFAULT_OUTPUT);
/// assert_eq!(query.distance_mode(), SubgraphQuery::DEFAULT_DISTANCE_MODE);
///
/// let query = query.with_limit(Some(100));
/// assert_eq!(query.limit(), Some(100));
///
/// let query = query.with_context(1000);
/// assert_eq!(query.context(), 1000);
///
/// let query = query.with_snarls(SnarlOutput::Contained);
/// assert_eq!(query.snarls(), SnarlOutput::Contained);
///
/// let query = query.with_snarls(SnarlOutput::Overlapping);
/// assert_eq!(query.snarls(), SnarlOutput::Overlapping);
///
/// let query = query.with_distance_mode(DistanceMode::Node);
/// assert_eq!(query.distance_mode(), DistanceMode::Node);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubgraphQuery {
    query_type: QueryType,

    // Optional safety limit for the size of the subgraph in nodes.
    limit: Option<usize>,

    // Context size around the query position(s) in bp.
    context: usize,

    // What to do with top-level snarls overlapping with the subgraph.
    snarls: SnarlOutput,

    // How to output the haplotypes.
    output: HaplotypeOutput,

    // Distance calculation heuristic for context extraction.
    distance_mode: DistanceMode,
}

impl SubgraphQuery {
    /// Default value for context length (in bp).
    pub const DEFAULT_CONTEXT: usize = 100;

    /// Default value for the snarl extraction option.
    pub const DEFAULT_SNARLS: SnarlOutput = SnarlOutput::None;

    /// Default value for the haplotype output option.
    pub const DEFAULT_OUTPUT: HaplotypeOutput = HaplotypeOutput::All;

    /// Default value for the distance calculation mode.
    pub const DEFAULT_DISTANCE_MODE: DistanceMode = DistanceMode::Side;

    /// Creates a query that retrieves a subgraph around a path offset.
    ///
    /// The reference path should be specified by using a sample name, a contig name, and optionally a haplotype number.
    /// The fragment field should not be used.
    /// If the reference haplotype is fragmented, the query will try to find the right fragment.
    pub fn path_offset(path_name: &FullPathName, offset: usize) -> Self {
        let mut path_name = path_name.clone();
        path_name.fragment = offset;
        SubgraphQuery {
            query_type: QueryType::PathOffset(path_name),
            limit: None,
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
            distance_mode: Self::DEFAULT_DISTANCE_MODE,
        }
    }

    /// Cretes a query that retrieves a subgraph around a path interval.
    ///
    /// The reference path should be specified by using a sample name, a contig name, and optionally a haplotype number.
    /// The fragment field should not be used.
    /// If the reference haplotype is fragmented, the query will try to find the right fragment.
    pub fn path_interval(path_name: &FullPathName, interval: Range<usize>) -> Self {
        let mut path_name = path_name.clone();
        path_name.fragment = interval.start;
        SubgraphQuery {
            query_type: QueryType::PathInterval(path_name, interval.len()),
            limit: None,
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
            distance_mode: Self::DEFAULT_DISTANCE_MODE,
        }
    }

    /// Creates a query that retrieves a subgraph around a set of nodes.
    pub fn nodes(nodes: impl IntoIterator<Item = usize>) -> Self {
        SubgraphQuery {
            query_type: QueryType::Nodes(nodes.into_iter().collect()),
            limit: None,
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
            distance_mode: Self::DEFAULT_DISTANCE_MODE,
        }
    }

    /// Creates a query that extracts a subgraph between two handles in the same chain.
    ///
    /// This query ignores context length and the snarl extraction flag.
    /// If the nodes are not in the same chain in the given order, the subgraph can otherwise be arbitrarily large.
    pub fn between(start: usize, end: usize) -> Self {
        SubgraphQuery {
            query_type: QueryType::Between(start, end),
            limit: None,
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
            distance_mode: Self::DEFAULT_DISTANCE_MODE,
        }
    }

    /// Returns an updated query with the given limit for the size of the subgraph in nodes.
    pub fn with_limit(self, limit: Option<usize>) -> Self {
        SubgraphQuery { limit, ..self }
    }

    /// Returns an updated query with the given context length.
    ///
    /// See [`Self::DEFAULT_CONTEXT`] for the default value.
    pub fn with_context(self, context: usize) -> Self {
        SubgraphQuery { context, ..self }
    }

    /// Returns an updated query with the given snarl extraction option.
    ///
    /// See [`Self::DEFAULT_SNARLS`] for the default value.
    pub fn with_snarls(self, snarls: SnarlOutput) -> Self {
        SubgraphQuery { snarls, ..self }
    }

    /// Returns an updated query with the given distance calculation mode.
    ///
    /// See [`Self::DEFAULT_DISTANCE_MODE`] for the default value.
    pub fn with_distance_mode(self, distance_mode: DistanceMode) -> Self {
        SubgraphQuery { distance_mode, ..self }
    }

    #[deprecated(since = "0.6.0", note = "Use `with_haplotypes` instead")]
    pub fn with_output(self, output: HaplotypeOutput) -> Self {
        self.with_haplotypes(output)
    }

    /// Returns an updated query with the given haplotype output option.
    ///
    /// See [`Self::DEFAULT_OUTPUT`] for the default value.
    ///
    /// # Panics
    ///
    /// Panics if this is a node-based query and the output would be [`HaplotypeOutput::ReferenceOnly`].
    pub fn with_haplotypes(self, output: HaplotypeOutput) -> Self {
        if self.is_node_based() {
            assert!(output != HaplotypeOutput::ReferenceOnly, "Reference-only output is not supported for node-based queries");
        }
        SubgraphQuery { output, ..self }
    }

    pub(super) fn query_type(&self) -> &QueryType {
        &self.query_type
    }

    /// Returns `true` if this is a reference-based / path-based query.
    pub fn is_reference_based(&self) -> bool {
        matches!(self.query_type, QueryType::PathOffset(_) | QueryType::PathInterval(_, _))
    }

    /// Returns `true` if this is a node-based query.
    pub fn is_node_based(&self) -> bool {
        matches!(self.query_type, QueryType::Nodes(_) | QueryType::Between(_, _))
    }

    /// Returns `true` if this is a multi-node query.
    pub fn is_multi_node(&self) -> bool {
        matches!(self.query_type, QueryType::Nodes(ref nodes) if nodes.len() > 1)
    }

    /// Returns the safety limit for the size of the subgraph in nodes.
    pub fn limit(&self) -> Option<usize> {
        self.limit
    }

    /// Returns the context length (in bp) for the query.
    pub fn context(&self) -> usize {
        self.context
    }

    /// Returns the snarl extraction option for the query.
    pub fn snarls(&self) -> SnarlOutput {
        self.snarls
    }

    /// Returns the output format for the query.
    pub fn output(&self) -> HaplotypeOutput {
        self.output
    }

    /// Returns the distance calculation mode for the query.
    pub fn distance_mode(&self) -> DistanceMode {
        self.distance_mode
    }
}

impl Display for SubgraphQuery {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Query itself.
        match self.query_type() {
            QueryType::PathOffset(path_name) => write!(f, "(path {}", path_name)?,
            QueryType::PathInterval(path_name, len) => write!(f, "(path {}, len {}", path_name, len)?,
            QueryType::Nodes(nodes) => write!(f, "(nodes {:#?}", nodes)?,
            QueryType::Between(start, end) => {
                let (start_id, start_o) = support::decode_node(*start);
                let (end_id, end_o) = support::decode_node(*end);
                write!(f, "(between ({} {}) and ({} {})", start_id, start_o, end_id, end_o)?;
            }
        }

        // Safety limit.
        if let Some(limit) = self.limit() {
            write!(f, ", limit {}", limit)?;
        }

        // Context length and snarls.
        match self.snarls() {
            SnarlOutput::None => write!(f, ", context {}", self.context)?,
            SnarlOutput::Contained => write!(f, ", context {} with contained snarls", self.context)?,
            SnarlOutput::Overlapping => write!(f, ", context {} with overlapping snarls", self.context)?,
        }

        // Distance mode (only if not the default).
        if self.distance_mode != Self::DEFAULT_DISTANCE_MODE {
            write!(f, ", distance mode {}", self.distance_mode)?;
        }

        // Haplotype output.
        write!(f, ", {})", self.output)?;

        Ok(())
    }
}

//-----------------------------------------------------------------------------
