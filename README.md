# GBZ in SQLite

This is a prototype for storing a [GBZ graph](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md) in an SQLite database.
It is intended for interactive applications, where you want to access parts of the graph immediately without loading the entire graph into memory.

Once the implementation has stabilized, it will be included as an optional feature in [GBWT-rs](https://github.com/jltsiren/gbwt-rs).

## Building

To build the package, run:

```
cargo build --release
```

You will then have the `gbz2db` and `query` tool binaries in `target/release/`.

## Building for WebAssembly

To build for a WASI WebAssembly runtime, run:

```
./build-wasm.sh
```

This will make sure you have the `wasm32-wasi` target installed via `rustup`, and produce `gbz2db.wasm` and `query.wasm` files in `target/wasm32-wasi/release/`.

#### Publishing to NPM

After building the `.wasm` files, and increasing the version in `package.json`, you can publish to [the `gbz-base` NPM package](https://www.npmjs.com/package/gbz-base).

First do a dry run and make sure that you are publishing the files you want to publish:
```
npm publish --dry-run
```

Then, publish for real:
```
npm publish
```

#### Using from NPM

To use the published NPM module in a project, first install it:

```
npm install --save gbz-base
```

The `gbz2db.wasm` and `query.wasm` files are available for `import` as `gbz-db/gbz2db.wasm` and `gbz2db/query.wasm`. If you are targeting the browser using Webpack, you can `await import()` them and then get a fetchable URL in the `default` field of the result:

```
let blobImport = await import("gbz-base/query.wasm");
let blob = await fetch(blobImport.default);
```

You can then use a [library like `browser-wasi-shim`](https://github.com/bjorn3/browser_wasi_shim#readme) to execute the command-line tool with particular arguments in a particular filesystem.

If you are not targeting the browser, you will have to open and read the files at their actual filesystem paths, `node_modules/gbz-base/target/wasm32-wasi/release/gbz2db.wasm` and `node_modules/gbz-base/target/wasm32-wasi/release/query.wasm`. Node itself cannot `import` or `require()` a WASI file in a useful way, and this module doesn't itself include any JavaScript code to help find and load the WASM files.
