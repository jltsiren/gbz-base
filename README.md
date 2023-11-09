# GBZ in SQLite

This is a prototype for storing a [GBZ graph](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md) in an SQLite database.
It is intended for interactive applications, where you want to access parts of the graph immediately without loading the entire graph into memory.

Once the implementation has stabilized, it will be included as an optional feature in [GBWT-rs](https://github.com/jltsiren/gbwt-rs).
