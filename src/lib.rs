//! # GBZ-base: an immutable pangenome graph stored in a SQLite database
//!
//! This is a prototype for storing a GBZ graph in a SQLite database.
//! It is intended for interactive applications that need immediate access to the graph.
//! In such applications, the overhead from loading the GBZ graph into memory can be significant (e.g. 20 seconds for a human graph).
//! As long as the application needs only a fraction of the entire graph (e.g. 1 Mbp context in a human graph), using the database is faster than loading the graph.
//! This assumes that the database is stored on a local SSD.
//!
//! The prototype builds on the [`gbwt`] crate.
//! Once the implementation has stabilized, it will become an optional feature in the [`gbwt`] crate.
//!
//! ### Basic concepts
//!
//! Nodes are accessed by handles, which are [`gbwt::GBWT`] node identifiers.
//! A handle encodes both the identifier of the node in the underlying graph and its orientation.
//! Each node record corresponds to a row in table `Nodes`, with the handle as its primary key.
//!
//! Paths are accessed by handles, which are path identifiers in the original graph.
//! Each path record corresponds to a row in table `Paths`, with the handle as its primary key.
//! The record contains information for both orientations of the path.
//!
//! Paths can be indexed for random access, which can be useful for e.g. finding a graph region by its reference coordinates.
//! Indexing is based on storing the sequence offset and the GBWT position at the start of a node once every ~1000 bp.
//! The indexed positions are stored in table `ReferenceIndex`.
//! By default, only generic paths (sample name `_gbwt_ref`) and reference paths (sample name listed in GBWT tag `reference_samples`) are indexed.
//! The database can become excessively large if all paths are indexed.

pub mod db;
pub mod formats;
pub mod subgraph;

pub use db::{GBZBase, GBZPath, GBZRecord, GraphInterface};
pub use subgraph::{Subgraph, SubgraphQuery, HaplotypeOutput, PathIndex};
