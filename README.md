# GBZ and GAF in SQLite

This is a prototype for SQLite-based file formats for:

* Pangenome graphs in [GBZ format](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).
* Sequence alignments to a pangenome graph in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).

The formats are intended for interactive applications, where you want to access parts of the graph immediately without loading the entire graph into memory.

Both file formats are under development and can change without warning.

## Building

To build the package, run:

```sh
cargo build --release
```

You will then have the following binaries in `target/release/`:

* `gbz2db`: Builds a GBZ-base.
* `gaf2db`: Builds a GAF-base.
* `db2gaf`: Converts GAF-base back to GAF.
* `query`: Queries in GBZ-base and GAF-base.
* `gafsort`: Sorts a GAF file for GAF-base construction.

## GBZ-base construction

### Basic construction

You can convert a GBZ graph `graph.gbz` into a database `graph.gbz.db` using:

```sh
gbz2db graph.gbz
```

The output file can be changed using `--output filename.db`.
If the output file already exists, the construction will fail.
An existing database can be overwritten with option `--overwrite`.

The database will be functionally equivalent to the GBZ graph, except that it will not contain a node-to-segment translation.
Generic paths (with sample name `_gbwt_ref`) and reference paths (samples specified in GBWT tag `reference_samples`) will be indexed for querying.

### Precomputed top-level chains

A GBZ-base stores links between the boundary nodes of snarls in top-level chains.
Such links enable better subgraph queries (see below).
By default, the `gbz2db` tries to find them using a simple algorithm that works with Minigraph–Cactus graphs.
It requires that each graph component contains exactly two tips with a directed path between them.

If the assumptions fail, the algorithm will not be able to find top-level chains for some components.
In such cases, prebuilt chains can be used with option `--chains graph.chains`.
[vg](https://github.com/vgteam/vg) version 1.69.0 or newer can extract a chains file from a distance index or a snarls file.

Example with chains extracted from a distance index:

```sh
vg chains graph.gbz graph.dist > graph.chains
gbz2db --chains graph.chains graph.gbz
```

## GAF-base construction

### Sorting the reads

Before building a GAF-base, you must first sort the alignments and build a GBWT index of the target paths.
This can be done using the bundled `gafsort` tool (or with `vg gamsort`):

```sh
gafsort --progress --threads 6 reads.gaf.gz | bgzip --threads 4 > sorted.gaf.gz
```

`sorted.gaf.gz` now contains the sorted reads.
Six worker threads and four compression threads should be sufficient due to sequential bottlenecks.

The default chunk size (1 million lines) is appropriate for short reads.
When sorting long reads, it can be changed with `--chunk-size N` (e.g. `--chunk-size 10k` for 20 kpb reads).
Chunk size can be specified using suffixes such as `k` or `M`.

### Building the GAF-base

GAF-base construction takes the sorted reads:

```sh
gaf2db sorted.gaf.gz
```

The reads can be uncompressed or compressed with gzip.
Default output file is `<input>.db`.
Options `--output` and `--overwrite` are the same as in GBZ-base construction.

A prebuilt GBWT index (e.g. from `vg gamsort`) can be provided with `--gbwt FILE`.
This will lower the memory requirements for GAF-base construction, but it will not make the construction faster.

The default block size (1000 alignments) is appropriate for short reads.
When building a database of long read alignments, it can be changed with `--block-size N` (e.g. `--block-size 10` for 20 kbp reads).

The default GAF-base is reference-based, like the GAF format itself.
The alignments can only be decoded by using the corresponding GBZ graph or GBZ-base.

A reference-free GAF-base can be built by providing a reference graph during construction with option `--ref-free FILE`.
The graph file should be a GBZ graph.
A GBZ-base can also be used, but the construction will be much slower.
Reference-free GAF-bases can be used without a reference graph.

If quality strings are not required, it is possible to drop them with option `--no-quality`.
This will make the database much smaller, particularly with long reads.
Similarly, optional fields not supported directly by GAF-base can be dropped with option `--no-optional`.

## Subgraph queries

The `query` tool supports extracting subgraphs from a database or a GBZ graph.
A database is available immediately, while loading a large graph into memory can take tens of seconds.
On the other hand, large and/or repeated queries can be tens of times faster with an in-memory graph.

Basic usage:

```sh
# Offset-based query.
query --sample GRCh38 --contig chr12 --offset 1234567 graph.db > out.gfa

# Interval-based query on a generic path.
query --contig chrM --interval 3000..4000 graph.db > out.gfa

# Node-based query.
query --node 12345 graph.db > out.gfa
```

Offsets are 0-based and intervals are half-open.
There will be proper metadata for the path used for an offset-based or interval-based query.
All other paths will be listed as unknown haplotypes.

### Other options

* `--handle N`: Node-based query using handles (GBWT node identifiers) encoded as `2 * node_id + is_reverse`.
* `--context N`: Extract `N` bp context around the query position (default: 100).
  Context length can be specified using suffixes such as `k` or `M`.
* `--distinct`: Collapse identical paths in the subgraph and report the number of copies using `WT:i` tags.
* `--reference-only`: Output the query path but no unknown haplotypes.
* `--cigar`: Output CIGAR strings relative to the query path as `CG:Z` tags.
* `--format json`: Extract the subgraph in JSON format instead of GFA.

### Snarl-based queries

Queries based on paths and nodes find the query position and then extract a greedy context around it.
This does not always result in a meaningful subgraph.
Nodes outside the region of interest will be included, if the shortest path to them is within the context.
Nodes within the region of interest may be excluded, if the region contains large enough variants.

Snarl-based queries avoid these issues:

```sh
# Snarl-based query.
query --between 12345:12401 graph.db > out.gfa
```

The query specifies two boundary nodes, which are assumed to be in the same chain in the snarl decomposition.
If no orientation is provided  (e.g. `12345+` or `12401-`), the boundary nodes are assumed to be in the forward orientation.
The query extracts all nodes and snarls in the chain between (and including) the boundary nodes but no greedy context.

If the boundary nodes are not in the same chain or they are given in the wrong order, the outcome is unpredictable.
To avoid extracting most of the chromosome, a safety limit for the number of nodes may be given with `--limit N`.
The limit can be specified using suffixes such as `k` or `M`.
If the limit is exceeded, the query will fail.

Other queries can also be made snarl-aware:

```sh
query --contig chrM --interval 3000..4000 --snarls graph.db > out.gfa
```

After extracting a subgraph, the query determines all top-level snarls that have their boundary nodes in the extracted subgraph.
Then it ensures that those snarls are fully in the subgraph.
These snarls are based on the top-level chains provided or computed during GBZ-base construction.

`--extend-snarls` handles two additional cases that `--snarls` does not:

* **Partial overlap**: the query interval covers only one boundary node of a top-level snarl.
  The subgraph is extended to fully include such snarls.
* **Interval inside snarls**: the query interval starts and/or ends inside a top-level snarl, so one or both outer chain boundary nodes do not appear in the initial subgraph.
  The path is walked backward from the interval start and forward from the interval end to find the nearest chain boundary on each side, and the subgraph is extended to include them. This covers both the case where the interval lies entirely within a single snarl and the case where its endpoints lie in different snarls.

```sh
query --contig chrM --interval 3000..4000 --extend-snarls graph.db > out.gfa
```

`--extend-snarls` implies `--snarls` and should be used with care, as some snarls can be very large.

### Extracting alignments

If you have a GBZ-base for a graph and a GAF-base for reads aligned to the graph, you can extract the reads aligned to the subgraph:

```sh
query --sample GRCh38 --contig chr12 --offset 1234567 \
    --gaf-base reads.db --gaf-output out.gaf \
    graph.db > out.gfa
```

By default, this extracts all alignments overlapping with the subgraph and clips them to the subgraph.
Use option `--alignments overlapping` to avoid clipping or `--alignments contained` to select only alignments fully within the subgraph.

The GBZ-base can be for the graph the reads were aligned to, or for any supergraph.
For example, a GBZ-base for a clipped (default) Minigraph–Cactus graph can be used with reads aligned to a corresponding frequency-filtered or personalized (haplotype-sampled) graph.

## Interface

See `cargo doc --open`.

## Building for WebAssembly

To build for a WASI WebAssembly runtime, run:

```sh
./build-wasm.sh
```

This will make sure you have the `wasm32-wasi` target installed via `rustup`, and produce `gbz2db.wasm` and `query.wasm` files in `target/wasm32-wasi/release/`.

#### Publishing to NPM

After building the `.wasm` files, and increasing the version in `package.json`, you can publish to [the `gbz-base` NPM package](https://www.npmjs.com/package/gbz-base).

First do a dry run and make sure that you are publishing the files you want to publish:

```sh
npm publish --dry-run
```

Then, publish for real:

```sh
npm publish
```

#### Using from NPM

To use the published NPM module in a project, first install it:

```sh
npm install --save gbz-base
```

The `gbz2db.wasm` and `query.wasm` files are available for `import` as `gbz-db/gbz2db.wasm` and `gbz2db/query.wasm`. If you are targeting the browser using Webpack, you can `await import()` them and then get a fetchable URL in the `default` field of the result:

```
let blobImport = await import("gbz-base/query.wasm");
let blob = await fetch(blobImport.default);
```

You can then use a [library like `browser-wasi-shim`](https://github.com/bjorn3/browser_wasi_shim#readme) to execute the command-line tool with particular arguments in a particular filesystem.

If you are not targeting the browser, you will have to open and read the files at their actual filesystem paths, `node_modules/gbz-base/target/wasm32-wasi/release/gbz2db.wasm` and `node_modules/gbz-base/target/wasm32-wasi/release/query.wasm`. Node itself cannot `import` or `require()` a WASI file in a useful way, and this module doesn't itself include any JavaScript code to help find and load the WASM files.
