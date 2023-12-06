use super::*;

use gbwt::GBWT;
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
            assert_eq!(record.handle, handle, "Wrong node record handle");
            if orientation == Orientation::Forward {
                assert_eq!(record.sequence, sequence, "Wrong sequence for handle {}", handle);
            } else {
                let rc = support::reverse_complement(sequence);
                assert_eq!(record.sequence, rc, "Wrong sequence for handle {}", handle);
            }

            // Edges and BWT.
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

//-----------------------------------------------------------------------------

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

// TODO: get_path, find_path, paths_for_sample

//-----------------------------------------------------------------------------

// TODO: get_indexed_position

//-----------------------------------------------------------------------------
