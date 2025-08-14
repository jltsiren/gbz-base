# GBZ and GAF in SQLite

This is a prototype for SQLite-based file formats for:

* Pangenome graphs in [GBZ graph](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).
* Sequence alignments to a pangenome graph in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).

The formats are intended for interactive applications, where you want to access parts of the graph immediately without loading the entire graph into memory.

Both file formats are under development and can change without warning.

## Building

To build the package, run:

```sh
cargo build --release
```

You will then have the `gbz2db`, `gaf2db`, and `query` tool binaries in `target/release/`.

## GBZ-base construction

You can convert a GBZ graph `graph.gbz` into a database `graph.gbz.db` using:

```sh
gbz2db graph.gbz
```

The output file can be changed using `--output filename.db`.
If the output file already exists, the construction will fail.
An existing database can be overwritten with option `--overwrite`.

The database will be functionally equivalent to the GBZ graph, except that it will not contain a node-to-segment translation.
Generic paths (with sample name `_gbwt_ref`) and reference paths (samples specified in GBWT tag `reference_samples`) will be indexed for querying.

## GAF-base construction

### Sorting the reads

Before building a GAF-base, you must first sort the alignments and build a GBWT index of the target paths.
This can be done using [vg](https://github.com/vgteam/vg):

```sh
vg gamsort --progress --threads 6 \
    --gbwt-output sorted.gbwt \
    --gaf-input reads.gaf.gz | bgzip --threads 4 > sorted.gaf.gz
```

`sorted.gaf.gz` now contains the sorted reads and `sorted.gbwt` the target paths in the same order.
Six worker threads are generally enough when reading a gzip-compressed GAF file.
Four compression threads should be sufficient when building a GBWT index.

The default chunk size (1 million lines) is appropriate for short reads.
When sorting long reads, it can be changed with `--chunk-size N` (e.g. `--chunk-size 10000` for 20 kpb reads).

### Building the GAF-base

GAF-base construction takes both the sorted reads and the GBWT index:

```sh
gaf2db --gbwt sorted.gbwt sorted.gaf.gz
```

The reads can be uncompressed or compressed with gzip.
Default output file is `<input>.db`.
Options `--output` and `--overwrite` are the same as in GBZ-base construction.

The default block size (1000 alignments) is appropriate for short reads.
When building a database of long read alignments, it can be changed with `--block-size N` (e.g. `--block-size 10` for 20 kbp reads).

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

* `--context N`: Extract `N` bp context around the query position (default: 100).
* `--distinct`: Collapse identical paths in the subgraph and report the number of copies using `WT:i` tags.
* `--reference-only`: Output the query path but no unknown haplotypes.
* `--cigar`: Output CIGAR strings relative to the query path as `CG:Z` tags.
* `--format json`: Extract the subgraph in JSON format instead of GFA.

### Extracting alignments

If you have a GBZ-base for a graph and a GAF-base for reads aligned to the graph, you can extract the reads aligned to the subgraph:

```sh
query --sample GRCh38 --contig chr12 --offset 1234567 \
    --gaf-base reads.db --gaf-output out.gaf \
    graph.db > out.gfa
```

By default, this extracts all alignments overlapping with the subgraph.
Use option `--contained` to extract only alignments fully within the subgraph.

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
