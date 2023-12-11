use super::*;

use gbwt::{GBWT, Metadata};
use std::path::PathBuf;

//-----------------------------------------------------------------------------

fn create_database_from_graph(graph: &GBZ) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!GBZBase::exists(&db_file), "Database {} already exists", db_file.display());
    let result = GBZBase::create(&graph, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

fn create_database_from_file(filename: &PathBuf) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!GBZBase::exists(&db_file), "Database {} already exists", db_file.display());
    let result = GBZBase::create_from_file(&filename, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

fn open_database(filename: &PathBuf) -> GBZBase {
    let database = GBZBase::open(filename);
    assert!(database.is_ok(), "Failed to open database: {}", database.unwrap_err());
    database.unwrap()
}

fn create_interface(database: &GBZBase) -> GraphInterface {
    let interface = GraphInterface::new(database);
    assert!(interface.is_ok(), "Failed to create graph interface: {}", interface.unwrap_err());
    interface.unwrap()
}

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

fn path_by_metadata(interface: &mut GraphInterface, sample: &str, contig: &str, haplotype: usize, fragment: usize) -> Option<GBZPath> {
    let path = interface.find_path(sample, contig, haplotype, fragment);
    assert!(path.is_ok(), "Failed to get path by metadata: {}", path.unwrap_err());
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

fn check_header(database: &GBZBase, graph: &GBZ) {
    let metadata = graph.metadata().unwrap();
    assert_eq!(database.nodes(), graph.nodes(), "Wrong number of nodes");
    assert_eq!(database.samples(), metadata.samples(), "Wrong number of samples");
    assert_eq!(database.haplotypes(), metadata.haplotypes(), "Wrong number of haplotypes");
    assert_eq!(database.contigs(), metadata.contigs(), "Wrong number of contigs");
    assert_eq!(database.paths(), metadata.paths(), "Wrong number of paths");
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

    // GBZ tags.
    for (key, value) in graph.tags().iter() {
        let tag = interface.get_gbz_tag(key);
        assert!(tag.is_ok(), "Failed to get GBZ tag: {}", tag.unwrap_err());
        let tag = tag.unwrap();
        assert!(tag.is_some(), "Missing GBZ tag {}", key);
        assert_eq!(tag.unwrap(), *value, "Wrong GBZ tag value for {}", key);
    }
}

fn check_nodes(interface: &mut GraphInterface, graph: &GBZ) {
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

            // Handle and sequence.
            assert_eq!(record.handle(), handle, "Wrong node record handle");
            if orientation == Orientation::Forward {
                assert_eq!(record.sequence(), sequence, "Wrong sequence for handle {}", handle);
            } else {
                let rc = support::reverse_complement(sequence);
                assert_eq!(record.sequence(), rc, "Wrong sequence for handle {}", handle);
            }

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
            assert!(db_record.is_some(), "Failed to convert node record for handle {} to GBWT record", handle);
            let db_record = db_record.unwrap();
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
fn create_from_graph() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_graph(&graph);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn create_from_file() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_file(&gbz_file);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

#[test]
fn get_path() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_graph(&graph);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
        assert_eq!(path.sample, sample_name, "Wrong sample name for path {}", path_handle);
        let contig_name = metadata.contig(path_name.contig()).unwrap();
        assert_eq!(path.contig, contig_name, "Wrong contig name for path {}", path_handle);
        assert_eq!(path.haplotype, path_name.phase(), "Wrong haplotype number for path {}", path_handle);
        assert_eq!(path.fragment, path_name.fragment(), "Wrong fragment number for path {}", path_handle);

        let found_path = path_by_metadata(&mut interface, &path.sample, &path.contig, path.haplotype, path.fragment);
        assert!(found_path.is_some(), "Could not find path {} by metadata", path_handle);
        let found_path = found_path.unwrap();
        assert_eq!(found_path, path, "Wrong path found by metadata for path {}", path_handle);
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

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_graph(&graph);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
    let db_file = create_database_from_file(&gbz_file);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

    // Paths by handle.
    let by_handle = path_by_handle(&mut interface, database.paths());
    assert!(by_handle.is_none(), "Found a nonexistent path by handle");

    // Paths by metadata.
    let by_metadata = path_by_metadata(&mut interface, "fake", "path", 0, 0);
    assert!(by_metadata.is_none(), "Found a nonexistent path by metadata");

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

#[test]
fn indexed_position() {
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_graph(&graph);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

    // Check that paths corresponding to `gbwt::REF_SAMPLE` have been indexed.
    // TODO: We should also check reference paths specified in the GFA file, but the example
    // does not have any.
    const MAX_LEN: usize = 20; // This is enough for the example.
    let gbwt: &GBWT = graph.as_ref();
    for path_handle in 0..database.paths() {
        let path = existing_path_by_handle(&mut interface, path_handle);
        if path.is_indexed {
            assert_eq!(path.sample, REF_SAMPLE, "Path {} with sample name {} is indexed", path_handle, path.sample);
            let mut previous = 0;
            let visited = visited_positions(gbwt, path_handle);
            for offset in 0..MAX_LEN {
                let (indexed, pos) = existing_indexed_pos(&mut interface, path_handle, offset);
                assert!(indexed <= offset, "Indexed position {} is after query position {} for path {}", indexed, offset, path_handle);
                assert!(indexed >= previous, "Indexed position {} is before the previous position {} for path {}", indexed, previous, path_handle);
                previous = indexed;
                assert!(visited.contains(&pos), "Path {} does not visit indexed position ({}, {}) at offset {}", path_handle, pos.node, pos.offset, indexed);
            }
        } else {
            assert_ne!(path.sample, REF_SAMPLE, "Path {} with sample name {} is not indexed", path_handle, path.sample);
            let result = indexed_pos(&mut interface, path_handle, MAX_LEN);
            assert!(result.is_none(), "Found indexed position for non-indeded path {}", path_handle);
            continue;
        }
    }

    let nonexistent = indexed_pos(&mut interface, database.paths(), MAX_LEN);
    assert!(nonexistent.is_none(), "Found indexed position for nonexistent path");

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------
