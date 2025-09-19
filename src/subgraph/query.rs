//! Queries for extracting a subgraph from GBZ-base or a GBZ graph.

use gbwt::{support, FullPathName};

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
pub(super) enum QueryType {
    // Path name and offset in bp stored in the fragment field.
    PathOffset(FullPathName),
    // Starting position as in `PathOffset` and length in bp.
    PathInterval(FullPathName, usize),
    // Set of node identifiers.
    Nodes(BTreeSet<usize>),
    // Subgraph between two handles in the same chain, with an optional safety limit for the number of nodes extracted.
    Between((usize, usize), Option<usize>),
}

//-----------------------------------------------------------------------------

/// Arguments for extracting a subgraph.
///
/// # Examples
///
/// ```
/// use gbz_base::SubgraphQuery;
/// use gbwt::FullPathName;
///
/// let path_name = FullPathName::generic("path");
/// let query = SubgraphQuery::path_offset(&path_name, 123);
/// assert_eq!(query.context(), SubgraphQuery::DEFAULT_CONTEXT);
/// assert_eq!(query.snarls(), SubgraphQuery::DEFAULT_SNARLS);
/// assert_eq!(query.output(), SubgraphQuery::DEFAULT_OUTPUT);
///
/// let query = query.with_context(1000);
/// assert_eq!(query.context(), 1000);
///
/// let query = query.with_snarls(true);
/// assert!(query.snarls());
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubgraphQuery {
    query_type: QueryType,

    // Context size around the reference position (in bp).
    context: usize,

    // Also extract nodes in covered top-level snarls.
    snarls: bool,

    // How to output the haplotypes.
    output: HaplotypeOutput,
}

impl SubgraphQuery {
    /// Default value for context length (in bp).
    pub const DEFAULT_CONTEXT: usize = 100;

    /// Default value for the snarl extraction flag.
    pub const DEFAULT_SNARLS: bool = false;

    /// Default value for the haplotype output option.
    pub const DEFAULT_OUTPUT: HaplotypeOutput = HaplotypeOutput::All;

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
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
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
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
        }
    }

    /// Creates a query that retrieves a subgraph around a set of nodes.
    pub fn nodes(nodes: impl IntoIterator<Item = usize>) -> Self {
        SubgraphQuery {
            query_type: QueryType::Nodes(nodes.into_iter().collect()),
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
        }
    }

    /// Creates a query that extracts a subgraph between two handles in the same chain.
    ///
    /// This query ignores context length and the snarl extraction flag.
    /// An optional safety limit for the size of the subgraph in nodes can be provided.
    /// If the nodes are not in the same chain in the given order, the subgraph can otherwise be arbitrarily large.
    pub fn between(start: usize, end: usize, limit: Option<usize>) -> Self {
        SubgraphQuery {
            query_type: QueryType::Between((start, end), limit),
            context: Self::DEFAULT_CONTEXT,
            snarls: Self::DEFAULT_SNARLS,
            output: Self::DEFAULT_OUTPUT,
        }
    }

    /// Returns an updated query with the given context length.
    ///
    /// See [`Self::DEFAULT_CONTEXT`] for the default value.
    pub fn with_context(self, context: usize) -> Self {
        SubgraphQuery { context, ..self }
    }

    /// Returns an updated query with the given snarl extraction flag.
    ///
    /// See [`Self::DEFAULT_SNARLS`] for the default value.
    pub fn with_snarls(self, snarls: bool) -> Self {
        SubgraphQuery { snarls, ..self }
    }

    /// Returns an updated query with the given haplotype output option.
    ///
    /// See [`Self::DEFAULT_OUTPUT`] for the default value.
    ///
    /// # Panics
    ///
    /// Panics if this is a node-based query and the output would be [`HaplotypeOutput::ReferenceOnly`].
    pub fn with_output(self, output: HaplotypeOutput) -> Self {
        if let QueryType::Nodes(_) = self.query_type {
            assert!(output != HaplotypeOutput::ReferenceOnly, "Reference-only output is not supported for node-based queries");
        }
        SubgraphQuery { output, ..self }
    }

    pub(super) fn query_type(&self) -> &QueryType {
        &self.query_type
    }

    /// Returns the context length (in bp) for the query.
    pub fn context(&self) -> usize {
        self.context
    }

    /// Returns `true` if the query also extracts nodes coverered top-level snarls.
    ///
    /// A snarl is covered, if both of its boundary nodes are contained in the query interval or in the context.
    pub fn snarls(&self) -> bool {
        self.snarls
    }

    /// Returns the output format for the query.
    pub fn output(&self) -> HaplotypeOutput {
        self.output
    }
}

impl Display for SubgraphQuery {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let context_str = if self.snarls() {
            format!("{} with snarls", self.context)
        } else {
            format!("{}", self.context)
        };
        match self.query_type() {
            QueryType::PathOffset(path_name) => write!(f, "(path {}, context {}, {})", path_name, context_str, self.output),
            QueryType::PathInterval(path_name, len) => write!(f, "(path {}, len {}, context {}, {})", path_name, len, context_str, self.output),
            QueryType::Nodes(nodes) => write!(f, "(nodes {:#?}, context {}, {})", nodes, context_str, self.output),
            QueryType::Between((start, end), limit) => {
                let (start_id, start_o) = support::decode_node(*start);
                let (end_id, end_o) = support::decode_node(*end);
                let limit_str = if let Some(limit) = limit {
                    format!(", limit {}", limit)
                } else {
                    String::new()
                };
                write!(f, "(between ({} {}) and ({} {}){}, {})", start_id, start_o, end_id, end_o, limit_str, self.output)
            },
        }
    }
}

//-----------------------------------------------------------------------------