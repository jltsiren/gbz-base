use crate::{HaplotypeOutput, SubgraphQuery, Alignment};
use crate::utils;

use super::*;

use gbwt::REF_SAMPLE;
use gbwt::{Metadata, GBWT};

use std::collections::HashSet;
use std::path::PathBuf;

//-----------------------------------------------------------------------------

fn create_database_from_graph(graph: &GBZ, chains: &Chains) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!utils::file_exists(&db_file), "Database {} already exists", db_file.display());
    let result = GBZBase::create(&graph, chains, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

fn create_database_from_files(gbz_file: &PathBuf, chains_file: Option<&PathBuf>) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!utils::file_exists(&db_file), "Database {} already exists", db_file.display());
    let chains_file = chains_file.map(|x| x.as_ref());
    let result = GBZBase::create_from_files(&gbz_file, chains_file, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

fn open_database(filename: &PathBuf) -> GBZBase {
    let database = GBZBase::open(filename);
    assert!(database.is_ok(), "Failed to open database: {}", database.unwrap_err());
    database.unwrap()
}

fn create_interface(database: &GBZBase) -> GraphInterface<'_> {
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

    // GBZ tags.
    for (key, value) in graph.tags().iter() {
        let tag = interface.get_gbz_tag(key);
        assert!(tag.is_ok(), "Failed to get GBZ tag: {}", tag.unwrap_err());
        let tag = tag.unwrap();
        assert!(tag.is_some(), "Missing GBZ tag {}", key);
        assert_eq!(tag.unwrap(), *value, "Wrong GBZ tag value for {}", key);
    }
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
    let db_file = create_database_from_graph(&graph, &chains);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
    let db_file = create_database_from_files(&gbz_file, None);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
    let db_file = create_database_from_graph(&graph, &chains);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

    // Check header, tags, and nodes.
    check_header(&database, &graph, &chains);
    check_tags(&mut interface, &graph);
    check_nodes(&mut interface, &graph, &chains);

    drop(interface);
    drop(database);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

// FIXME handle empty paths
#[test]
fn get_path() {
    // Load the graph.
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let chains = Chains::new();

    // Create and open the database and create a graph interface.
    let db_file = create_database_from_graph(&graph, &chains);
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
    let db_file = create_database_from_graph(&graph, &chains);
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
    let db_file = create_database_from_files(&gbz_file, None);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
    let db_file = create_database_from_graph(&graph, &chains);
    let database = open_database(&db_file);
    let mut interface = create_interface(&database);

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
    let ref_samples = vec![String::from(REF_SAMPLE)];
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

fn create_gaf_base(gaf_part: &'static str, gbwt_part: &'static str) -> PathBuf {
    let gaf_file = utils::get_test_data(gaf_part);
    let gbwt_file = utils::get_test_data(gbwt_part);
    let db_file = serialize::temp_file_name("gaf-base");
    let params = GAFBaseParams::default();
    let result = GAFBase::create_from_files(&gaf_file, &gbwt_file, &db_file, &params);
    assert!(result.is_ok(), "Failed to create GAF-base database: {}", result.unwrap_err());
    db_file
}

fn open_gaf_base(db_file: &Path) -> GAFBase {
    let db = GAFBase::open(db_file);
    assert!(db.is_ok(), "Failed to open GAF-base database: {}", db.unwrap_err());
    // We do not check the filename, as the one reported by SQLite may differ.
    db.unwrap()
}

fn check_gaf_base(db: &GAFBase, nodes: usize, alignments: usize, bidirectional_gbwt: bool) {
    assert_eq!(db.nodes(), nodes, "Wrong number of nodes in GAF-base database");
    assert_eq!(db.alignments(), alignments, "Wrong number of alignments in GAF-base database");
    assert_eq!(db.bidirectional_gbwt(), bidirectional_gbwt, "Wrong GBWT bidirectionality in GAF-base database");
}

fn gaf_base_gbz() -> GBZ {
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let result = serialize::load_from(&gbz_file);
    assert!(result.is_ok(), "Failed to load GBZ graph: {}", result.unwrap_err());
    result.unwrap()
}

fn gaf_base_db() -> (GBZBase, PathBuf) {
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let db_file = create_database_from_files(&gbz_file, None);
    let db = open_database(&db_file);
    (db, db_file)
}

fn gaf_base_reads(gaf_part: &'static str) -> Vec<Alignment> {
    let gaf_file = utils::get_test_data(gaf_part);

    let reader = utils::open_file(gaf_file);
    assert!(reader.is_ok(), "Failed to open GAF file: {}", reader.err().unwrap());
    let mut reader = reader.unwrap();

    let mut result = Vec::new();
    let mut line_num = 1;
    loop {
        let mut buf: Vec<u8> = Vec::new();
        let len = reader.read_until(b'\n', &mut buf);
        assert!(len.is_ok(), "Failed to read line {}: {}", line_num, len.err().unwrap());
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.is_empty() {
            break;
        }
        let aln = Alignment::from_gaf(&buf);
        assert!(aln.is_ok(), "Failed to parse GAF line {}: {}", line_num, aln.unwrap_err());
        result.push(aln.unwrap());
        line_num += 1;
    };

    result

}

fn gaf_base_subgraph_queries() -> Vec<(SubgraphQuery, String)> {
    vec![
        (SubgraphQuery::nodes(vec![]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("empty")),
        (SubgraphQuery::nodes(vec![1000]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("single component")),
        (SubgraphQuery::nodes(vec![500, 1500]).with_context(100).with_snarls(false).with_output(HaplotypeOutput::All), String::from("two components")),
    ]
}

fn validate_read_set(read_set: &ReadSet, subgraph: &Subgraph, all_reads: &[Alignment], name: &str, contained: bool) {
    // Select the alignments that should be included in the read set.
    let mut expected_reads = Vec::new();
    for aln in all_reads {
        let target_path = aln.target_path().unwrap();
        if target_path.is_empty() {
            continue;
        }
        let mut take = false;
        for handle in target_path {
            if subgraph.has_handle(*handle) {
                take = true;
            } else if contained {
                take = false;
                break;
            }
        }
        if take {
            let mut aln = aln.clone();
            aln.optional.clear(); // We do not store unknown optional fields for the moment.
            expected_reads.push(aln);
        }
    }
    assert_eq!(read_set.len(), expected_reads.len(), "Wrong number of reads in read set for {}", name);

    // GAF-base construction and read set extraction maintain order.
    for (i, (aln, truth)) in read_set.iter().zip(expected_reads.iter()).enumerate() {
        assert_eq!(aln, truth, "Wrong read {} for {}", i, name);
    }

    // Generate the expected GAF output.
    let mut expected_gaf = Vec::new();
    for aln in expected_reads {
        let mut target_sequence = Vec::new();
        let target_path = aln.target_path().unwrap();
        for handle in target_path {
            let record = read_set.nodes.get(handle);
            assert!(record.is_some(), "Missing GBZ record for handle {} in {}", handle, name);
            let record = record.unwrap();
            target_sequence.extend_from_slice(record.sequence());
        }
        target_sequence = target_sequence[aln.path_interval.clone()].to_vec();
        let gaf_line = aln.to_gaf(&target_sequence);
        expected_gaf.extend_from_slice(&gaf_line);
        expected_gaf.push(b'\n');
    }

    // Check the GAF output.
    let mut gaf_output = Vec::new();
    let result = read_set.to_gaf(&mut gaf_output);
    assert!(result.is_ok(), "Failed to output GAF for {}: {}", name, result.unwrap_err());
    assert_eq!(gaf_output.len(), expected_gaf.len(), "Wrong GAF output length for {}", name);
    assert_eq!(gaf_output, expected_gaf, "Wrong GAF output for {}", name);
}

//-----------------------------------------------------------------------------

// Tests for GAF-base.

#[test]
fn gaf_base_empty() {
    let db_file = create_gaf_base("empty.gaf", "empty.gbwt");
    let db = open_gaf_base(&db_file);
    check_gaf_base(&db, 0, 0, true);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_plain_unidirectional() {
    let db_file = create_gaf_base("micb-kir3dl1_HG003.gaf", "micb-kir3dl1_HG003.gbwt");
    let db = open_gaf_base(&db_file);
    check_gaf_base(&db, 2291, 12439, false);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_plain_bidirectional() {
    let db_file = create_gaf_base("micb-kir3dl1_HG003.gaf", "bidirectional.gbwt");
    let db = open_gaf_base(&db_file);
    check_gaf_base(&db, 2324, 12439, true);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_gz_unidirectional() {
    let db_file = create_gaf_base("micb-kir3dl1_HG003.gaf.gz", "micb-kir3dl1_HG003.gbwt");
    let db = open_gaf_base(&db_file);
    check_gaf_base(&db, 2291, 12439, false);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

#[test]
fn gaf_base_gz_bidirectional() {
    let db_file = create_gaf_base("micb-kir3dl1_HG003.gaf.gz", "bidirectional.gbwt");
    let db = open_gaf_base(&db_file);
    check_gaf_base(&db, 2324, 12439, true);

    drop(db);
    let _ = std::fs::remove_file(&db_file);
}

//-----------------------------------------------------------------------------

// Tests for ReadSet.

fn test_read_set_gbz(gbwt_part: &'static str) {
    // Load GBZ.
    let graph = gaf_base_gbz();

    // Build and open GAF-base.
    let gaf_base_file = create_gaf_base("micb-kir3dl1_HG003.gaf", gbwt_part);
    let gaf_base = open_gaf_base(&gaf_base_file);

    // Parse the reads as a source of truth.
    let all_reads = gaf_base_reads("micb-kir3dl1_HG003.gaf");

    let queries = gaf_base_subgraph_queries();
    for (query, name) in queries {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_gbz(&graph, None, None, &query);
        assert!(result.is_ok(), "Failed to extract subgraph for {}: {}", name, result.unwrap_err());

        let overlapping = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &gaf_base, false);
        assert!(overlapping.is_ok(), "Failed to extract overlapping reads for {}: {}", name, overlapping.unwrap_err());
        let overlapping = overlapping.unwrap();
        validate_read_set(&overlapping, &subgraph, &all_reads, &name, false);

        let contained = ReadSet::new(GraphReference::Gbz(&graph), &subgraph, &gaf_base, true);
        assert!(contained.is_ok(), "Failed to extract fully contained reads for {}: {}", name, contained.unwrap_err());
        let contained = contained.unwrap();
        validate_read_set(&contained, &subgraph, &all_reads, &name, true);
    }

    // Cleanup.
    drop(gaf_base);
    let _ = std::fs::remove_file(&gaf_base_file);
}

#[test]
fn read_set_gbz_unidirectional() {
    test_read_set_gbz("micb-kir3dl1_HG003.gbwt");
}

#[test]
fn read_set_gbz_bidirectional() {
    test_read_set_gbz("bidirectional.gbwt");
}

fn test_read_set_db(gbwt_part: &'static str) {
    // Build and open GAF-base.
    let gaf_base_file = create_gaf_base("micb-kir3dl1_HG003.gaf", gbwt_part);
    let gaf_base = open_gaf_base(&gaf_base_file);

    // Open GBZ-base.
    let (gbz_base, gbz_base_file) = gaf_base_db();
    let mut graph = create_interface(&gbz_base);

    // Parse the reads as a source of truth.
    let all_reads = gaf_base_reads("micb-kir3dl1_HG003.gaf");

    let queries = gaf_base_subgraph_queries();
    for (query, name) in queries {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_db(&mut graph, &query);
        assert!(result.is_ok(), "Failed to extract subgraph for {}: {}", name, result.unwrap_err());

        let overlapping = ReadSet::new(GraphReference::Db(&mut graph), &subgraph, &gaf_base, false);
        assert!(overlapping.is_ok(), "Failed to extract overlapping reads for {}: {}", name, overlapping.unwrap_err());
        let overlapping = overlapping.unwrap();
        validate_read_set(&overlapping, &subgraph, &all_reads, &name, false);

        let contained = ReadSet::new(GraphReference::Db(&mut graph), &subgraph, &gaf_base, true);
        assert!(contained.is_ok(), "Failed to extract fully contained reads for {}: {}", name, contained.unwrap_err());
        let contained = contained.unwrap();
        validate_read_set(&contained, &subgraph, &all_reads, &name, true);
    }

    // Cleanup.
    drop(graph);
    drop(gbz_base);
    let _ = std::fs::remove_file(&gbz_base_file);
    drop(gaf_base);
    let _ = std::fs::remove_file(&gaf_base_file);
}

#[test]
fn read_set_db_unidirectional() {
    test_read_set_db("micb-kir3dl1_HG003.gbwt");
}

#[test]
fn read_set_db_bidirectional() {
    test_read_set_db("bidirectional.gbwt");
}

//-----------------------------------------------------------------------------
