# GBZ-base releases

## Current version

* New functionality:
  * `Subgraph::extract_haplotype_walks` for extracting the paths as `HaplotypeWalk` objects.
  * Tuning database behavior with `execute_pragma` and `set_pragma` in `GBZBase` and `GAFBase`.
* Compatibility with the old subgraph query algorithms still used in GBWTGraph:
  * `DistanceMode` parameter for extracting context based on shortest distances to node sides (default) or nodes (for compatibility).
  * `AlignmentMode` parameter for aligning other paths to the reference by LCS weighted by node length (default) or unweighted LCS (for compatibility).
  * These parameters are not exposed in `gbz-base query`, and they may be deprecated and removed in future releases.
* Database files can be opened using `file:` URIs.
* Support for GAF version 1.1.

## GBZ-base 0.6.2 (2026-09-26)

* Query benchmark can also use GBZ graphs in addition to GBZ-bases.

This is a patch version for the revised paper.

## GBZ-base 0.6.1 (2026-08-24)

* Support for GBZ version 3 and GBWT version 6 with Zstandard compressed BWT.

## GBZ-base 0.6.0 (2026-08-11)

* Database versions: GBZ-base version 4, GAF-base version 4
* Fallible operations return `Error` instead of `String`.
  * An `Error` combines a message with an `ErrorKind` that tells apart failures needing different responses.
  * `Result<T>` is a shorthand for `Result<T, Error>`.
* Consolidated binaries:
  * `gbz-base construct` replaces `gbz2db`.
  * `gbz-base query` replaces `query`.
  * `gaf-base sort` replaces `gafsort`.
  * `gaf-base construct` replaces `gaf2db`.
  * `gaf-base decompress` replaces `db2gaf`.
  * Old binaries are deprecated but still available.
* Parameter presets for short and long reads in `gaf-base sort` and `gaf-base construct`.
* Query improvements:
  * Haplotype output selection with `--haplotypes`, with an option for no haplotypes.
  * Option `--gaf-only` for writing the alignments instead of the subgraph to stdout.
  * Safety limit for subgraph size can be set for all query types.
* Bug fixes:
  * `gaf-base decompress` works correctly with a reference-free GAF-base.
  * Multithreaded stable GAF sorting works correctly.

## GBZ-base 0.5.1 (2026-06-02)

* GBZ-base construction prints a warning if it cannot find top-level chains for all graph components.
* Tool for benchmarking GBZ-base and GAF-base queries.

This is a patch version for the paper.

## GBZ-base 0.5.0 (2026-05-05)

* Database versions: GBZ-base version 4, GAF-base version 4
* GBZ-base version 4:
  * Same as version 0.4.0, but with the new version number scheme.
* GAF-base version 4:
  * Quality strings are compressed with rANS 4x16 (order-1, rle).
  * GAF-base is now substantially smaller than bgzip-compressed sorted GAF.
* Improved snarl-aware queries:
  * `SnarlOutput::Contained`: Extend the subgraph with snarls that have both boundary nodes in the subgraph (existing behavior).
  * `SnarlOutput::Overlapping`: Extend the subgraph with all overlapping snarls (requires a connected subgraph).

## GBZ-base 0.4.0 (2026-04-17)

* Database versions: GBZ-base v0.4.0, GAF-base version 3
* GBZ-base and GAF-base construction without vg:
  * GBZ-base construction will find top-level chains automatically if not provided.
  * GAF sorting with `gaf_sort()` and the `gafsort` tool.
  * GBWT construction in a background thread when building GAF-base.
* Bug fixes:
  * Exact alignments in a block with varying query lengths are now encoded correctly. GAF-bases must be rebuilt.
* Binaries:
  * Command line arguments can use suffixes (e.g. `k`, `MiB`) for sizes and counts that can plausibly be large.
  * Binaries report peak memory usage.

## GAF-base 0.3.0 (2026-02-18)

* Database versions: GBZ-base v0.4.0, GAF-base version 3
* Support for GBZ version 2 with Zstandard compressed sequences.
* GAF-base version 3:
  * More space-efficient representation of numerical values in table `Alignments`.
  * Database construction parameters stored in table `Tags`.
  * Optional reference-free GAF-base by storing node sequences in table `Nodes`.
  * Option to leave out base quality strings.
  * Also stores unknown optional fields, unless told otherwise.

## GBZ-base 0.2.0 (2025-12-26)

* Database versions: GBZ-base v0.4.0, GAF-base v0.2.0
* `db2gaf` tool for converting a GAF-base back to GAF format.
* Support for [stable graph names](https://github.com/jltsiren/pggname):
  * Uses `pggname::GraphName` for importing graph names and relationships between GBZ tags and GFA/GAF headers.
  * Stable graph names are stored in databases and included in GFA/GAF outputs.
  * `db2gaf` and `query` use the information for determining if the graph is a valid reference for the alignments.

## GBZ-base 0.1.0 (2025-11-06)

* Database versions: GBZ-base v0.4.0, GAF-base v0.2.0

This is the initial release of GBZ-base and GAF-base.

## Release process

* Clean up with `cargo clean`.
* Update database versions to non-dev versions in `db.rs`.
* Update version in `Cargo.toml`.
* Switch to crates.io versions of dependencies, if necessary.
* Update `RELEASES.md`.
* Run `cargo clippy --features benchmark`.
* Run tests with `cargo test`.
* Build documentation with `cargo doc`.
* Build the optimized version with `cargo build --release --features benchmark`.
* Commit the final changes for the release.
* Publish in crates.io with `cargo publish`.
* Push to GitHub.
* Draft a new release in GitHub.
