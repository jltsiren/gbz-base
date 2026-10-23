# GBZ and GAF in SQLite

SQLite-based file formats for:

* Pangenome graphs in [GBZ format](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).
* Sequence alignments to a pangenome graph in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).

The formats are intended for interactive applications, where you want to access parts of the graph immediately without loading the entire graph into memory.
For best performance, the database files should be on a fast local SSD.
The GBZ-base file format has stabilized, while GAF-base is still under active development.

## Building

To build the package, run:

```sh
cargo build --release
```

You will then have the following binaries in `target/release/`:

* `gbz-base`:
  * `gbz-base construct` builds a GBZ-base (`gbz2db` until version 0.5.1).
  * `gbz-base query` for querying GBZ-bases and GAf-bases (`query` until version 0.5.1).
* `gaf-base`:
  * `gaf-base sort` sorts a GAF file for GAF-base construction (`gafsort` until version 0.5.1).
  * `gaf-base construct` builds a GAF-base (`gaf2db` until version 0.5.1).
  * `gaf-base decompress` converts GAF-base back to GAF (`db2gaf` until version 0.5.1).

## GBZ-base construction

### Basic construction

You can convert a GBZ graph `graph.gbz` into a database `graph.gbz.db` using:

```sh
gbz-base construct graph.gbz
```

The output file can be changed using `--output filename.db`.
If the output file already exists, the construction will fail.
An existing database can be overwritten with option `--overwrite`.

The database will be functionally equivalent to the GBZ graph, except that it will not contain a node-to-segment translation.
Generic paths (with sample name `_gbwt_ref`) and reference paths (samples specified in GBWT tag `reference_samples`) will be indexed for querying.

### Precomputed top-level chains

A GBZ-base stores links between the boundary nodes of snarls in top-level chains.
Such links enable better subgraph queries (see below).
By default, the `gbz-base construct` tries to find them using a simple algorithm that works with Minigraph–Cactus graphs.
It requires that each graph component contains exactly two tips with a directed path between them.

If the assumptions fail, the algorithm will not be able to find top-level chains for some components.
It will then print a warning message.
In such cases, prebuilt chains can be used with option `--chains graph.chains`.
[vg](https://github.com/vgteam/vg) version 1.69.0 or newer can extract a chains file from a distance index or a snarls file.
A top-level distance index (built with `vg index` option `--no-nested-distance`) is sufficient for this purpose.

Example with chains extracted from a distance index:

```sh
vg chains graph.gbz graph.dist > graph.chains
gbz-base construct --chains graph.chains graph.gbz
```

## GAF-base construction

### Sorting the reads

Before building a GAF-base, you must first sort the alignments and build a GBWT index of the target paths.
This can be done using `gaf-base sort` (or with `vg gamsort`):

```sh
gaf-base sort --progress --threads 6 reads.gaf.gz | bgzip --threads 4 > sorted.gaf.gz
```

`sorted.gaf.gz` now contains the sorted reads.
Six worker threads and four compression threads should be sufficient due to sequential bottlenecks.

The default temporary file size (1 million lines) is appropriate for short reads.
You can get a more appropriate size (10 thousand lines) for long reads with `--preset long`.
File size can be changed manually with `--records-per-file N` (e.g. `--records-per-file 10k` for 20 kpb reads).
The size can be specified using suffixes such as `k` or `M`.

### Building the GAF-base

GAF-base construction takes the sorted reads:

```sh
gaf-base construct sorted.gaf.gz
```

The reads can be uncompressed or compressed with gzip.
Default output file is `<input>.db`.
Options `--output` and `--overwrite` are the same as in GBZ-base construction.

A prebuilt GBWT index (e.g. from `vg gamsort`) can be provided with `--gbwt FILE`.
This will lower the memory requirements for GAF-base construction, but it will not make the construction faster.

The default block size (1000 alignments) is appropriate for short reads.
With long reads, you can use smaller blocks (10 alignments) with `--preset long`.
Block size can be set manually with `--block-size N` (e.g. `--block-size 10` for 20 kbp reads).

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

`gbz-base query` supports extracting subgraphs from a database or a GBZ graph.
A database is available immediately, while loading a large graph into memory can take tens of seconds.
On the other hand, large and/or repeated queries can be tens of times faster with an in-memory graph.

Basic usage:

```sh
# Offset-based query.
gbz-base query --sample GRCh38 --contig chr12 --offset 1234567 graph.db > out.gfa

# Interval-based query on a generic path.
gbz-base query --contig chrM --interval 3000..4000 graph.db > out.gfa

# Node-based query.
gbz-base query --node 12345 graph.db > out.gfa
```

Offsets are 0-based and intervals are half-open.
There will be proper metadata for the path used for an offset-based or interval-based query.
All other paths will be listed as unknown haplotypes.

### Other options

* `--handle N`: Node-based query using handles (GBWT node identifiers) encoded as `2 * node_id + is_reverse`.
* `--limit N`: Abort the query if subgraph size exceeds `N` nodes.
  The limit can be specified using suffixes such as `k` or `M`.
* `--context N`: Extract `N` bp context around the query position (default: 100).
  Context length can be specified using suffixes such as `k` or `M`.
* `--haplotypes X`: Haplotype output selection:
  * `all`: Output all haplotypes paths (default).
  * `distinct`: Collapse identical paths in the subgraph and report the number of copies using `WT:i` tags.
  * `reference-only`: Output the query path but no unknown haplotypes.
  * `none`: Do not include haplotype paths in the output.
* `--cigar`: Output CIGAR strings relative to the query path as `CG:Z` tags.
* `--format json`: Extract the subgraph in JSON format instead of GFA.

### Extracting a sequence of snarls

Queries based on paths and nodes find the query position and then extract a greedy context around it.
This does not always result in a meaningful subgraph.
Nodes outside the region of interest will be included, if the shortest path to them is within the context.
Nodes within the region of interest may be excluded, if the region contains large enough variants.

Snarl-based queries avoid these issues:

```sh
# Snarl-based query.
gbz-base query --between 12345:12401 graph.db > out.gfa
```

The query specifies two boundary nodes, which are assumed to be in the same chain in the snarl decomposition.
If no orientation is provided  (e.g. `12345+` or `12401-`), the boundary nodes are assumed to be in the forward orientation.
The query extracts all nodes and snarls in the chain between (and including) the boundary nodes but no greedy context.

If the boundary nodes are not in the same chain or they are given in the wrong order, the outcome is unpredictable.
A safety limit (`--limit`) can be useful for dealing with such situations.

### Extracting snarls covered by a query

Other queries can also be made snarl-aware:

```sh
gbz-base query --contig chrM --interval 3000..4000 --snarls graph.db > out.gfa
```

After extracting a subgraph, the query extends it to cover all top-level snarls contained in the subgraph.
A top-level snarl is contained in a subgraph if both of its boundary nodes are in the subgraph.
These snarls are based on the top-level chains provided or computed during GBZ-base construction.

The extended subgraph can be much larger than the original one, if:

* The subgraph is not a single component.
* The reference sample contains a large deletion in the subgraph.
* Context length (`--context`) is non-zero and there is a large deletion in any sample within the context.

A safety limit (`--limit`) can be used to deal with such situations.

### Extracting snarls overlapping with a query

We can also extend the subgraph with all overlapping snarls:

```sh
gbz-base query --contig chrM --interval 3000..4000 --extend-snarls graph.db > out.gfa
```

In this case, we also include partially overlapping snarls.
A snarl is partially overlapping with a subgraph if a boundary node and one of its successors in the snarl are in the subgraph.
If there are no chain links in the subgraph, the query instead tries to find a snarl that contains the entire subgraph.

Option `--extend-snarls` requires that the subgraph is weakly connected.
It does not work with node-based queries with multiple nodes.
If the graph contains large snarls, the extended subgraph can be much larger than the original one.
A safety limit (`--limit`) can be useful here as well.

### Extracting alignments

If you have a GBZ-base for a graph and a GAF-base for reads aligned to the graph, you can extract the reads aligned to the subgraph:

```sh
gbz-base query --sample GRCh38 --contig chr12 --offset 1234567 \
    --gaf-base reads.db --gaf-output out.gaf \
    graph.db > out.gfa
```

By default, this extracts all alignments overlapping with the subgraph and clips them to the subgraph.
Use option `--alignments overlapping` to avoid clipping or `--alignments contained` to select only alignments fully within the subgraph.

If the subgraph is not needed, GAF-only output (to stdout) can be selected with option `--gaf-only`:

```sh
gbz-base query --sample GRCh38 --contig chr12 --offset 1234567 \
    --gaf-base reads.db --gaf-only \
    graph.db > out.gaf
```

The GBZ-base can be for the graph the reads were aligned to, or for any supergraph.
For example, a GBZ-base for a clipped (default) Minigraph–Cactus graph can be used with reads aligned to a corresponding frequency-filtered or personalized (haplotype-sampled) graph.

## Test data

Directory `test-data` contains test data.
In particular, there is a subset of a HPRC release 1 graph also used as the [haplotype sampling test case](https://github.com/vgteam/vg/blob/master/test/haplotype-sampling/README.md) in vg.
The test graph contains the following regions:

* `GRCh38#chr6:31498145-31511124` (micb)
* `GRCh38#chr19:54816468-54830778` (kir3dl1)

The relevant files are:

* `micb-kir3dl1.gbz` / `micb-kir3dl1.gfa`: The graph
* `micb-kir3dl1.chains`: Top-level chains for the graph
* `micb-kir3dl1_HG003.gaf` / `micb-kir3dl1_HG003.gaf.gz`: HG003 reads aligned to the graph
* `micb-kir3dl1_HG003.gbwt`: Pre-built GBWT index for the reads

## Batched GAF-base queries

`gaf-base batch` can be used for querying GAF-base without the overhead of opening a database connection and extracting a subgraph.
The subcommand reads queries from standard input and writes GAF output to standard output.
It can be used interactively (typically by launching it as a child process) or for batched queries.
In interactive use, the parent process is expected to issue one query at a time and read the result afterwards.

Each query is given as two integers that define a closed range of node identifiers.
An empty line or end-of-file tells `gaf-base batch` to stop.
The results from a single query are similarly terminated by an empty line.

Example queries:
```
179551168	179552548
179552548	179554243
179554243	179558906
179558906	179559573
179559573	179560180
```

Basic usage:
```sh
gaf-base batch reads.db < queries.tsv > out.gaf
```

If the GAF-base is reference-based, the output will be written as GAF version 1.1.
As a reference-based database does not store node sequences, target path lengths will be undefined in the output.
Difference strings will use `?` and `!` operations for mismatches and deletions, respectively.
The `*` and `-` operations used in GAF 1.0 would require target bases.

### Options

* `-p` / `--progress`: Print progress information to standard error.
* `--ids`: Include 1-based database row ids and 0-based row offsets as GAF tags `ri:i` and `ro:i`, respectively.
  With block size `b`, the 0-based row number in the GAF file should be `(ri - 1) * b + ro`, excluding headers.
* `--aln-cache N`: Cache size in bytes for the database connection used for alignments.
  This can be specified using suffixes such as `50M` or `128MiB`.
* `--node-cache N`: Cache size in bytes for the database connection used for node records.
  This can be specified using suffixes such as `50M` or `128MiB`.
* `--record-cache N`: Cache size threshold in handles for reusing the node record cache in the next query.
  This can be specified using suffixes such as `100k`.
  When subsequent queries cover adjacent graph regions, this can speed them up significantly.

## Interface

See `cargo doc --open`.

## Building for WebAssembly

**This section is probably outdated. We will update it once we have the time to continue working on a version of Sequence Tube Map backed with GBZ-base.**

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

## Citing GBZ-base and GAF-base

Jouni Sirén, Benedict Paten, and the Human Pangenome Reference Consortium: **GBZ-base and GAF-base: Indexed pangenome file formats**.
bioRxiv, 2026.
DOI: [10.64898/2026.07.10.737775](https://doi.org/10.64898/2026.07.10.737775)

## Notes

* The included `.cargo/config.toml` sets the target CPU to `native`.
