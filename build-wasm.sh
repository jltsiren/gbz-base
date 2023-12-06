#!/usr/bin/env bash

# build-wasm.sh: Build for WebAssembly. Must run from the root of the project.

# Get the toolchain
rustup target add wasm32-wasi

if [[ -z "${CC_wasm32_wasi}" ]] ; then
    if [[ ! -e wasi-sdk-20.0/bin/clang ]] ; then
        echo >&2 "Installing WASI C SDK..."
        if [[ "$(uname)" == "Darwin" ]] ; then
            # Get the toolchain for Mac
            curl -O https://github.com/WebAssembly/wasi-sdk/releases/download/wasi-sdk-20/wasi-sdk-20.0-macos.tar.gz
        else
            # Get the toolchain for Linux
            wget https://github.com/WebAssembly/wasi-sdk/releases/download/wasi-sdk-20/wasi-sdk-20.0-linux.tar.gz
        fi
        tar -xvf wasi-sdk-20.0-*.tar.gz
        rm wasi-sdk-20.0-*.tar.gz
    fi
    export CC_wasm32_wasi="$(pwd)/wasi-sdk-20.0/bin/clang"
fi

# Make the sqlite build not use long double (due to missing intrinsic implementations in the SDK libc.
# Also don't try and use pthreads for mutexes since we don't have a pthreads implementation.
export LIBSQLITE3_FLAGS="-DLONGDOUBLE_TYPE=double -DSQLITE_THREADSAFE=0"

# There is a wasm32-wasi-preview1-threads toolchain in Rust but that doesn't work with the SDK because the SDK calls it wasm32-wasi-threads. See <https://github.com/WebAssembly/wasi-libc/pull/434#issuecomment-1843395319>
cargo build --release --target=wasm32-wasi
