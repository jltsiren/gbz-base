use crate::{GBZBase, GraphInterface, formats};
use crate::{GAFBase, GAFBaseParams, GraphReference};
use crate::{Alignment, PathIndex, Chains};
use crate::utils;

use gbz::GBZ;
use gbz::support;

use simple_sds::serialize;

use std::path::{Path, PathBuf};

//-----------------------------------------------------------------------------

// GBZ-base utilities.

pub(crate) fn create_gbz_base_from_graph(graph: &GBZ, chains: &Chains) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!utils::file_exists(&db_file), "Database {} already exists", db_file.display());
    let result = GBZBase::create(&graph, chains, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

pub(crate) fn create_gbz_base_from_files(gbz_file: &Path, chains_file: Option<&Path>) -> PathBuf {
    let db_file = serialize::temp_file_name("gbz-base");
    assert!(!utils::file_exists(&db_file), "Database {} already exists", db_file.display());
    let chains_file = chains_file.map(|x| x.as_ref());
    let result = GBZBase::create_from_files(&gbz_file, chains_file, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    db_file
}

pub(crate) fn open_gbz_base(filename: &Path) -> GBZBase {
    let database = GBZBase::open(filename);
    assert!(database.is_ok(), "Failed to open database: {}", database.unwrap_err());
    database.unwrap()
}

pub(crate) fn create_graph_interface(database: &GBZBase) -> GraphInterface<'_> {
    let interface = GraphInterface::new(database);
    assert!(interface.is_ok(), "Failed to create graph interface: {}", interface.unwrap_err());
    interface.unwrap()
}

//-----------------------------------------------------------------------------

// GAF-base utilities.

pub(crate) fn create_gaf_base(gaf_part: &'static str, gbwt_part: &'static str) -> PathBuf {
    create_gaf_base_with_params(gaf_part, gbwt_part, GraphReference::None, &GAFBaseParams::default())
}

pub(crate) fn create_gaf_base_with_params(gaf_part: &'static str, gbwt_part: &'static str, graph: GraphReference<'_, '_>, params: &GAFBaseParams) -> PathBuf {
    let gaf_file = utils::get_test_data(gaf_part);
    let gbwt_file = utils::get_test_data(gbwt_part);
    let db_file = serialize::temp_file_name("gaf-base");
    let result = GAFBase::create_from_files(&gaf_file, &gbwt_file, &db_file, graph, params);
    assert!(result.is_ok(), "Failed to create GAF-base database: {}", result.unwrap_err());
    db_file
}

pub(crate) fn open_gaf_base(db_file: &Path) -> GAFBase {
    let db = GAFBase::open(db_file);
    assert!(db.is_ok(), "Failed to open GAF-base database: {}", db.unwrap_err());
    // We do not check the filename, as the one reported by SQLite may differ.
    db.unwrap()
}

//-----------------------------------------------------------------------------

// Files for the GAF-base test case (from vg haplotype sampling).

pub(crate) fn load_gaf_base_gbz() -> GBZ {
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let result = serialize::load_from(&gbz_file);
    assert!(result.is_ok(), "Failed to load GBZ graph: {}", result.unwrap_err());
    result.unwrap()
}

pub(crate) fn create_gaf_base_db() -> (GBZBase, PathBuf) {
    let gbz_file = utils::get_test_data("micb-kir3dl1.gbz");
    let db_file = create_gbz_base_from_files(&gbz_file, None);
    let db = open_gbz_base(&db_file);
    (db, db_file)
}

pub(crate) fn load_gaf_base_reads(gzip_compressed: bool) -> Vec<Alignment> {
    let gaf_part = if gzip_compressed {
        "micb-kir3dl1_HG003.gaf.gz"
    } else {
        "micb-kir3dl1_HG003.gaf"
    };
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
        if formats::is_gaf_header_line(&buf) {
            line_num += 1;
            continue;
        }
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

//-----------------------------------------------------------------------------

// Files for the GBWT-rs test case.

pub(crate) fn load_gbz(filename: &'static str) -> GBZ {
    let gbz_file = support::get_test_data(filename);
    let graph = serialize::load_from(&gbz_file);
    if let Err(err) = graph {
        panic!("Failed to load GBZ graph from {}: {}", gbz_file.display(), err);
    }
    graph.unwrap()
}

pub(crate) fn load_gbz_and_create_path_index(filename: &'static str, interval: usize) -> (GBZ, PathIndex) {
    let graph = load_gbz(filename);
    let path_index = PathIndex::new(&graph, interval, false);
    if let Err(err) = path_index {
        panic!("Failed to create path index with interval {}: {}", interval, err);
    }
    (graph, path_index.unwrap())
}

pub(crate) fn load_chains(filename: &'static str) -> Chains {
    let chains_file = utils::get_test_data(filename);
    let chains = Chains::load_from(&chains_file);
    if let Err(err) = chains {
        panic!("Failed to load chains from {}: {}", chains_file.display(), err);
    }
    chains.unwrap()
}

//-----------------------------------------------------------------------------
