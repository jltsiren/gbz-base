use super::*;

use crate::{internal, utils};

use gbz::GENERIC_SAMPLE;
use gbz::{Metadata, GBWT};

use std::collections::HashSet;
use std::path::PathBuf;

//-----------------------------------------------------------------------------

fn path_by_handle(interface: &mut GraphInterface, path_handle: usize) -> Option<GBZPath> {
    let path = interface.get_path(path_handle);
    assert!(path.is_ok(), "Failed to get path {}: {}", path_handle, path.unwrap_err());
    path.unwrap()
}

fn existing_path_by_handle(interface: &mut GraphInterface, path_handle: usize) -> GBZPath {
    let path = path_by_handle(interface, path_handle);
    assert!(path.is_some(), "Missing path {}", path_handle);
    path.unwrap()
}

fn path_by_name(interface: &mut GraphInterface, name: &FullPathName) -> Option<GBZPath> {
    let path = interface.find_path(name);
    assert!(path.is_ok(), "Failed to get path by name: {}", path.unwrap_err());
    path.unwrap()
}

fn paths_by_sample(interface: &mut GraphInterface, sample: &str) -> Vec<GBZPath> {
    let paths = interface.paths_for_sample(sample);
    assert!(paths.is_ok(), "Failed to get paths for sample {}: {}", sample, paths.unwrap_err());
    paths.unwrap()
}

fn indexed_pos(interface: &mut GraphInterface, path_handle: usize, offset: usize) -> Option<(usize, Pos)> {
    let result = interface.indexed_position(path_handle, offset);
    assert!(result.is_ok(), "Failed to get indexed position for path {} at offset {}: {}", path_handle, offset, result.unwrap_err());
    result.unwrap()
}

fn existing_indexed_pos(interface: &mut GraphInterface, path_handle: usize, offset: usize) -> (usize, Pos) {
    let result = indexed_pos(interface, path_handle, offset);
    assert!(result.is_some(), "Missing indexed position for path {} at offset {}", path_handle, offset);
    result.unwrap()
}

//-----------------------------------------------------------------------------

fn check_header(database: &GBZBase, graph: &GBZ, chains: &Chains) {
    let metadata = graph.metadata().unwrap();
    assert_eq!(database.nodes(), graph.nodes(), "Wrong number of nodes");
    assert_eq!(database.chains(), chains.len(), "Wrong number of chains");
    assert_eq!(database.chain_links(), chains.links(), "Wrong number of chain links");
    assert_eq!(database.paths(), metadata.paths(), "Wrong number of paths");
    assert_eq!(database.samples(), metadata.samples(), "Wrong number of samples");
    assert_eq!(database.haplotypes(), metadata.haplotypes(), "Wrong number of haplotypes");
    assert_eq!(database.contigs(), metadata.contigs(), "Wrong number of contigs");
}

fn check_tags(interface: &mut GraphInterface, graph: &GBZ) {
    let gbwt: &GBWT = graph.as_ref();

    // GBWT tags.
    for (key, value) in gbwt.tags().iter() {
        let tag = interface.get_gbwt_tag(key);
        assert!(tag.is_ok(), "Failed to get GBWT tag: {}", tag.unwrap_err());
        let tag = tag.unwrap();
        assert!(tag.is_some(), "Missing GBWT tag {}", key);
        assert_eq!(tag.unwrap(), *value, "Wrong GBWT tag value for {}", key);
    }
    let extracted_tags = interface.get_gbwt_tags();
    assert!(extracted_tags.is_ok(), "Failed to get GBWT tags: {}", extracted_tags.unwrap_err());
    let extracted_tags = extracted_tags.unwrap();
    assert_eq!(extracted_tags, *gbwt.tags(), "Wrong GBWT tags from the database");

    // GBZ tags.
    for (key, value) in graph.tags().iter() {
        let tag = interface.get_gbz_tag(key);
        assert!(tag.is_ok(), "Failed to get GBZ tag: {}", tag.unwrap_err());
        let tag = tag.unwrap();
        assert!(tag.is_some(), "Missing GBZ tag {}", key);
        assert_eq!(tag.unwrap(), *value, "Wrong GBZ tag value for {}", key);
    }
    let extracted_tags = interface.get_gbz_tags();
    assert!(extracted_tags.is_ok(), "Failed to get GBZ tags: {}", extracted_tags.unwrap_err());
    let extracted_tags = extracted_tags.unwrap();
    assert_eq!(extracted_tags, *graph.tags(), "Wrong GBZ tags from the database");
}

fn check_nodes(interface: &mut GraphInterface, graph: &GBZ, chains: &Chains) {
    let gbwt: &GBWT = graph.as_ref();
    let bwt: &BWT = gbwt.as_ref();

    for node_id in graph.node_iter() {
        let sequence = graph.sequence(node_id).unwrap();
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let handle = support::encode_node(node_id, orientation);
            let record = interface.get_record(handle);
            assert!(record.is_ok(), "Failed to get node record for handle {}: {}", handle, record.unwrap_err());
            let record = record.unwrap();
            assert!(record.is_some(), "Missing node record for handle {}", handle);
            let record = record.unwrap();

            // Handle, sequence, and link.
            assert_eq!(record.handle(), handle, "Wrong node record handle");
            if orientation == Orientation::Forward {
                assert_eq!(record.sequence(), sequence, "Wrong sequence for handle {}", handle);
            } else {
                let rc = support::reverse_complement(sequence);
                assert_eq!(record.sequence(), rc, "Wrong sequence for handle {}", handle);
            }
            assert_eq!(record.next(), chains.next(handle), "Wrong chain link for handle {}", handle);

            // Successors using the graph / iterator.
            let truth: Vec<(usize, Orientation)> = graph.successors(node_id, orientation).unwrap().collect();
            let successors: Vec<usize> = record.successors().collect();
            assert_eq!(successors.len(), truth.len(), "Wrong number of successors for handle {} using successors()", handle);
            let edges: Vec<Pos> = record.edges().collect();
            assert_eq!(edges.len(), truth.len(), "Wrong number of successors for handle {} using edges()", handle);
            for i in 0..successors.len() {
                let (truth_id, truth_orientation) = truth[i];
                let truth_handle = support::encode_node(truth_id, truth_orientation);
                assert_eq!(successors[i], truth_handle, "Wrong {}-th successor for handle {} using successors()", i, handle);
                assert_eq!(edges[i].node, truth_handle, "Wrong {}-th successor for handle {} using edges()", i, handle);
            }

            // Records: edges and BWT.
            let db_record = record.to_gbwt_record();
            let record_index = gbwt.node_to_record(handle);
            let gbwt_record = bwt.record(record_index).unwrap();
            assert_eq!(db_record.outdegree(), gbwt_record.outdegree(), "Wrong outdegree for handle {}", handle);
            for i in 0..gbwt_record.outdegree() {
                assert_eq!(db_record.successor(i), gbwt_record.successor(i), "Wrong successor {} for handle {}", i, handle);
                assert_eq!(db_record.offset(i), gbwt_record.offset(i), "Wrong offset {} for handle {}", i, handle);
            }
            assert_eq!(db_record.decompress(), gbwt_record.decompress(), "Wrong BWT for handle {}", handle);
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn gbz_record_from_graph() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let gbwt: &GBWT = graph.as_ref();
    let bwt: &BWT = gbwt.as_ref();

    // Check records for all handles, including ones that do not exist.
    // Only check the raw data and rely on other tests for GBZRecord methods.
    for handle in 0..=gbwt.alphabet_size() {
        let record = GBZRecord::from_gbz(&graph, handle);
        let node_id = support::node_id(handle);
        if !graph.has_node(node_id) {
            assert!(record.is_none(), "Found a record for nonexistent handle {}", handle);
            continue;
        } else {
            assert!(record.is_some(), "Missing record for handle {}", handle);
        }
        let record = record.unwrap();

        // Data from the GBWT record.
        assert_eq!(record.handle, handle, "Wrong handle for handle {}", handle);
        let record_id = gbwt.node_to_record(handle);
        let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
        let (edges, _) = Record::decompress_edges(edge_bytes).unwrap();
        assert_eq!(record.edges, edges, "Wrong edges for handle {}", handle);
        assert_eq!(record.bwt, bwt_bytes, "Wrong BWT for handle {}", handle);

        // Sequence.
        let sequence = graph.sequence(node_id).unwrap();
        if support::node_orientation(handle) == Orientation::Forward {
            assert_eq!(record.sequence, sequence, "Wrong sequence for handle {}", handle);
        } else {
            let rc = support::reverse_complement(sequence);
            assert_eq!(record.sequence, rc, "Wrong sequence for handle {}", handle);
        }
    }
}

//-----------------------------------------------------------------------------

#[test]
fn create_from_graph() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains = Chains::new();

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_graph(&graph, &chains);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph, &chains);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph, &chains);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn create_from_file() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains = Chains::new();

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_files(&gbz_file, None);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph, &chains);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph, &chains);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn large_test_case() {
    // Load the graph.
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains_file = utils::get_test_data("micb-kir3dl1.chains");
    let chains = Chains::load_from(&chains_file).unwrap();

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_graph(&graph, &chains);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph, &chains);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph, &chains);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

// TODO: handle empty paths
#[test]
fn get_path() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains = Chains::new();

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_graph(&graph, &chains);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check all fields except `is_indexed` for all paths. Also check that the
    // path can be found by its metadata.
    let gbwt: &GBWT = graph.as_ref();
    let metadata: &Metadata = gbwt.metadata().unwrap();
    for path_handle in 0..database.paths() {
        let path = existing_path_by_handle(&mut interface, path_handle);

        assert_eq!(path.handle, path_handle, "Wrong path handle");
        let fw_start = gbwt.start(support::encode_path(path_handle, Orientation::Forward)).unwrap();
        assert_eq!(path.fw_start, fw_start, "Wrong forward start position for path {}", path_handle);
        let rev_start = gbwt.start(support::encode_path(path_handle, Orientation::Reverse)).unwrap();
        assert_eq!(path.rev_start, rev_start, "Wrong reverse start position for path {}", path_handle);

        let path_name = metadata.path(path_handle).unwrap();
        let sample_name = metadata.sample(path_name.sample()).unwrap();
        assert_eq!(path.name.sample, sample_name, "Wrong sample name for path {}", path_handle);
        let contig_name = metadata.contig(path_name.contig()).unwrap();
        assert_eq!(path.name.contig, contig_name, "Wrong contig name for path {}", path_handle);
        assert_eq!(path.name.haplotype, path_name.phase(), "Wrong haplotype number for path {}", path_handle);
        assert_eq!(path.name.fragment, path_name.fragment(), "Wrong fragment number for path {}", path_handle);

        let found_path = path_by_name(&mut interface, &path.name);
        assert!(found_path.is_some(), "Could not find path {} by name", path_handle);
        let found_path = found_path.unwrap();
        assert_eq!(found_path, path, "Wrong path found by name for path {}", path_handle);

        // Now we know that `path` is correct, so we can use it to check `GBZPath` creation from a GBZ graph.
        let with_id = GBZPath::with_id(&graph, path_handle);
        assert!(with_id.is_some(), "Could not create GBZPath {} from GBZ using id", path_handle);
        let mut with_id = with_id.unwrap();
        with_id.is_indexed = path.is_indexed;
        assert_eq!(with_id, path, "Wrong GBZPath {} created from GBZ using id", path_handle);

        // And again, but using the path name.
        let with_name = GBZPath::with_name(&graph, &path.name);
        assert!(with_name.is_some(), "Could not create GBZPath {} from GBZ using name", path.name);
        let mut with_name = with_name.unwrap();
        with_name.is_indexed = path.is_indexed;
        assert_eq!(with_name, path, "Wrong GBZPath {} created from GBZ using name", path.name);
    }

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn paths_for_sample() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains = Chains::new();

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_graph(&graph, &chains);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check that we get correct paths for each existing sample name.
    let metadata = graph.metadata().unwrap();
    for sample_id in 0..metadata.samples() {
        let path_handles: Vec<usize> = metadata.path_iter().enumerate()
            .filter(|(_, path_name)| path_name.sample() == sample_id)
            .map(|(path_id, _)| path_id)
            .collect();
        let sample_name = metadata.sample(sample_id).unwrap();
        let paths = paths_by_sample(&mut interface, sample_name);
        assert_eq!(paths.len(), path_handles.len(), "Wrong number of paths for sample {}", sample_name);
        for (index, path_handle) in path_handles.iter().enumerate() {
            assert_eq!(*path_handle, paths[index].handle, "Wrong {}-th path handle for sample {}", index, sample_name);
        }
    }

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn nonexistent_paths() {
    // Create and open the database and create a graph interface.
    let gbz_file = support::get_test_data("example.gbz");
    let db_file = internal::create_gbz_base_from_files(&gbz_file, None);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Paths by handle.
    let by_handle = path_by_handle(&mut interface, database.paths());
    assert!(by_handle.is_none(), "Found a nonexistent path by handle");

    // Paths by name.
    let path_name = FullPathName::reference("fake", "path");
    let by_metadata = path_by_name(&mut interface, &path_name);
    assert!(by_metadata.is_none(), "Found a nonexistent path by name");

    // Paths by sample name.
    let by_sample = paths_by_sample(&mut interface, "fake");
    assert!(by_sample.is_empty(), "Found nonexistent paths by sample name");

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

fn visited_positions(gbwt: &GBWT, path_handle: usize) -> HashSet<Pos> {
    let encoded = support::encode_path(path_handle, Orientation::Forward);
    let mut curr = gbwt.start(encoded);
    let mut result = HashSet::new();
    while let Some(p) = curr {
        result.insert(p);
        curr = gbwt.forward(p);
    }
    result
}

fn check_indexed_positions(gbz_file: &PathBuf, chains_file: Option<&PathBuf>, ref_samples: Vec<String>) {
    let graph: GBZ = serialize::load_from(gbz_file).unwrap();
    let chains = if let Some(chains_file) = chains_file {
        Chains::load_from(chains_file).unwrap()
    } else {
        Chains::new()
    };

    // Create and open the database and create a graph interface.
    let db_file = internal::create_gbz_base_from_graph(&graph, &chains);
    let database = internal::open_gbz_base(&db_file);
    let mut interface = internal::create_graph_interface(&database);

    // Check that paths corresponding to the reference samples have been indexed.
    let mut offsets: Vec<usize> = (0..20).collect();
    let past_the_end = 1000000;
    offsets.push(past_the_end); // This is past the end in all test cases.
    let gbwt: &GBWT = graph.as_ref();
    for path_handle in 0..database.paths() {
        let path = existing_path_by_handle(&mut interface, path_handle);
        if path.is_indexed {
            assert!(ref_samples.contains(&path.name.sample), "Path {} with sample name {} is not indexed", path_handle, path.name.sample);
            let mut previous = 0;
            let visited = visited_positions(gbwt, path_handle);
            for &offset in offsets.iter() {
                let (indexed, pos) = existing_indexed_pos(&mut interface, path_handle, offset);
                assert!(indexed <= offset, "Indexed position {} is after query position {} for path {}", indexed, offset, path_handle);
                assert!(indexed >= previous, "Indexed position {} is before the previous position {} for path {}", indexed, previous, path_handle);
                previous = indexed;
                assert!(visited.contains(&pos), "Path {} does not visit indexed position ({}, {}) at offset {}", path_handle, pos.node, pos.offset, indexed);
            }
        } else {
            assert!(!ref_samples.contains(&path.name.sample), "Path {} with sample name {} is indexed", path_handle, path.name.sample);
            let result = indexed_pos(&mut interface, path_handle, past_the_end);
            assert!(result.is_none(), "Found indexed position for non-indeded path {}", path_handle);
            continue;
        }
    }

    let nonexistent = indexed_pos(&mut interface, database.paths(), past_the_end);
    assert!(nonexistent.is_none(), "Found indexed position for nonexistent path");

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn index_generic_paths() {
    let gbz_file = support::get_test_data("example.gbz");
    let chains_file = None;
    let ref_samples = vec![String::from(GENERIC_SAMPLE)];
    check_indexed_positions(&gbz_file, chains_file, ref_samples);
}

#[test]
fn index_reference_paths() {
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let chains_file = Some(utils::get_test_data("micb-kir3dl1.chains"));
    let ref_samples = vec![String::from("GRCh38"), String::from("CHM13")];
    check_indexed_positions(&gbz_file, chains_file.as_ref(), ref_samples);
}

//-----------------------------------------------------------------------------

// Helper functions for GAF-base tests.

fn check_gaf_base(db: &GAFBase, nodes: usize, alignments: usize, rows: usize, bidirectional_gbwt: bool) {
    assert_eq!(db.nodes(), nodes, "Wrong number of nodes in GAF-base database");
    assert_eq!(db.alignments(), alignments, "Wrong number of alignments in GAF-base database");
    assert_eq!(db.blocks(), rows, "Wrong number of blocks in GAF-base database");
    assert_eq!(db.bidirectional_gbwt(), bidirectional_gbwt, "Wrong GBWT bidirectionality in GAF-base database");
}

fn check_gaf_base_tags(db: &GAFBase, name_tags: usize) {
    let tags = db.tags();
    let expected_tags = 4 + name_tags; // version, nodes, alignments, bidirectional
    assert_eq!(tags.len(), expected_tags, "Wrong number of tags in GAF-base database");
    if name_tags > 0 {
        let graph_name = GraphName::from_tags(&tags);
        assert!(graph_name.is_ok(), "Failed to parse GraphName from GAF-base tags: {}", graph_name.unwrap_err());
        let graph_name = graph_name.unwrap();
        assert!(graph_name.name().is_some(), "Missing graph name in GAF-base tags");
    }
}

//-----------------------------------------------------------------------------

// Tests for GAF-base.

const UNIDIRECTIONAL_NODES: usize = 2291;
const BIDIRECTIONAL_NODES: usize = 2324;
const ALIGNMENTS: usize = 12439;
const BLOCKS: usize = 14; // 1000 reads/block, with unaligned reads as a separate block.

#[test]
fn gaf_base_empty() {
    let db_file = internal::create_gaf_base("empty.gaf", "empty.gbwt");
    let db = internal::open_gaf_base(&db_file);
    check_gaf_base(&db, 0, 0, 0, true);
    check_gaf_base_tags(&db, 0);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_plain_unidirectional() {
    let db_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt");
    let db = internal::open_gaf_base(&db_file);
    check_gaf_base(&db, UNIDIRECTIONAL_NODES, ALIGNMENTS, BLOCKS, false);
    check_gaf_base_tags(&db, 1);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_plain_bidirectional() {
    let db_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf", "bidirectional.gbwt");
    let db = internal::open_gaf_base(&db_file);
    check_gaf_base(&db, BIDIRECTIONAL_NODES, ALIGNMENTS, BLOCKS, true);
    check_gaf_base_tags(&db, 1);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_gz_unidirectional() {
    let db_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf.gz", "micb-kir3dl1_HG003.gbwt");
    let db = internal::open_gaf_base(&db_file);
    check_gaf_base(&db, UNIDIRECTIONAL_NODES, ALIGNMENTS, BLOCKS, false);
    check_gaf_base_tags(&db, 1);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_gz_bidirectional() {
    let db_file = internal::create_gaf_base("micb-kir3dl1_HG003.gaf.gz", "bidirectional.gbwt");
    let db = internal::open_gaf_base(&db_file);
    check_gaf_base(&db, BIDIRECTIONAL_NODES, ALIGNMENTS, BLOCKS, true);
    check_gaf_base_tags(&db, 1);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

