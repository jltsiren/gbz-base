# GBZ-base releases

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

* Run `cargo clippy`.
* Run tests with `cargo test`.
* Update database versions to non-dev versions in `db.rs`.
* Update version in `Cargo.toml`.
* Update `RELEASES.md`.
* Publish in crates.io with `cargo publish`.
* Push to GitHub.
* Draft a new release in GitHub.
