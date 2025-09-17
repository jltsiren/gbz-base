//! # GBZ-base and GAF-base: pangenome file formats using SQLite databases.
//!
//! # GBZ-base
//!
//! This is a prototype for storing a GBZ graph in a SQLite database.
//! It is intended for interactive applications that need immediate access to the graph.
//! In such applications, the overhead from loading the GBZ graph into memory can be significant (e.g. 20 seconds for a human graph).
//! As long as the application needs only a fraction of the entire graph (e.g. 1 Mbp context in a human graph), using the database is faster than loading the graph.
//! This assumes that the database is stored on a local SSD.
//!
//! The prototype builds on the [`gbwt`] crate.
//!
//! See [`GBZBase`], [`GraphInterface`], and [`Subgraph`] for the database interface.
//! See [`GBZPath`] and [`GBZRecord`] for the related structures.
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
//!
//! # GAF-base
//!
//! This is a prototype for an SQLite-based file format for sequence alignments to a pangenome graph.
//! It is mostly compatible with the GAF format.
//! Target paths are stored as a GBWT index in table `Nodes`, which is similar to the table in GBZ-base.
//!
//! Alignment metadata is stored in table `Alignments`.
//! Each row in the table corresponds to a block of alignments, which are assumed to be close in the graph.
//! The metadata is stored space-efficiently using column-based compression.
//!
//! `Alignment` table is indexed by (minimum handle, maximum handle) in the target paths.
//! Given a query subgraph, we want to find blocks where the interval overlaps with the subgraph.
//! These blocks must then be decompressed, as they may contain alignments to the subgraph.
//!
//! See [`GAFBase`] and [`ReadSet`] for the database interface.
//! See [`alignment`], [`Alignment`], and [`AlignmentBlock`] for more details.

pub mod alignment;
pub mod db;
pub mod formats;
pub mod path_index;
pub mod subgraph;
pub mod utils;

pub use alignment::{Alignment, AlignmentBlock};
pub use alignment::mapping::{Difference, Mapping};
pub use db::{GBZBase, GBZPath, GBZRecord, GraphInterface, GraphReference};
pub use db::{GAFBase, GAFBaseParams, ReadSet};
pub use path_index::PathIndex;
pub use subgraph::Subgraph;
pub use subgraph::query::{SubgraphQuery, HaplotypeOutput};
pub use utils::Chains;
