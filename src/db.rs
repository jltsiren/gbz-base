//! GBZ-base and GAF-base: SQLite databases storing a GBZ graph and sequence alignments to the graph.

use crate::{Alignment, AlignmentBlock, Chains};
use crate::{formats, utils};

use std::io::BufRead;
use std::path::Path;
use std::sync::{mpsc, Arc};
use std::{fs, thread};

use rusqlite::{Connection, OpenFlags, OptionalExtension, Row, Statement};

use gbz::{FullPathName, Orientation, Pos, GBWT, GBZ};
use gbz::bwt::{BWT, Record};
use gbz::support::{self, Tags};

use pggname::GraphName;

use simple_sds::serialize;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// A database connection to a GBZ-base database.
///
/// This structure stores a database connection and some header information.
/// In multi-threaded applications, each thread should have its own connection.
/// Graph operations are supported through the [`GraphInterface`] structure.
///
/// # Examples
///
/// ```
/// use gbz_base::{utils, GBZBase};
/// use gbz::support;
/// use simple_sds::serialize;
/// use std::fs;
///
/// // Create the database.
/// let gbz_file = support::get_test_data("example.gbz");
/// let db_file = serialize::temp_file_name("gbz-base");
/// assert!(!utils::file_exists(&db_file));
/// // Here we build a database without chains.
/// let result = GBZBase::create_from_files(&gbz_file, None, &db_file);
/// assert!(result.is_ok());
///
/// // Open the database and check some header information.
/// let database = GBZBase::open(&db_file).unwrap();
/// assert_eq!(database.nodes(), 12);
/// assert_eq!(database.samples(), 2);
/// assert_eq!(database.haplotypes(), 3);
/// assert_eq!(database.contigs(), 2);
/// assert_eq!(database.paths(), 6);
///
/// // Clean up.
/// drop(database);
/// fs::remove_file(&db_file).unwrap();
/// ```
#[derive(Debug)]
pub struct GBZBase {
    connection: Connection,
    version: String,
    nodes: usize,
    chains: usize,
    chain_links: usize,
    paths: usize,
    samples: usize,
    haplotypes: usize,
    contigs: usize,
}

/// Using the database.
impl GBZBase {
    /// Index positions at the start of a node on reference paths approximately every this many base pairs.
    pub const INDEX_INTERVAL: usize = 1000;

    // Key for database version.
    const KEY_VERSION: &'static str = "version";

    /// Current database version.
    pub const VERSION: &'static str = "GBZ-base v0.4.0";

    // Key for node count.
    const KEY_NODES: &'static str = "nodes";

    // Key for chain count.
    const KEY_CHAINS: &'static str = "chains";

    // Key for the total number of links in the chains.
    const KEY_CHAIN_LINKS: &'static str = "chain_links";

    // Key for path count.
    const KEY_PATHS: &'static str = "paths";

    // Key for sample count.
    const KEY_SAMPLES: &'static str = "samples";

    // Key for sample count.
    const KEY_HAPLOTYPES: &'static str = "haplotypes";

    // Key for sample count.
    const KEY_CONTIGS: &'static str = "contigs";

    // Prefix for GBWT tag keys.
    const KEY_GBWT: &'static str = "gbwt_";

    // Prefix for GBZ tag keys.
    const KEY_GBZ: &'static str = "gbz_";

    /// Opens a connection to the database in the given file.
    ///
    /// Reads the header information and passes through any database errors.
    pub fn open<P: AsRef<Path>>(filename: P) -> Result<Self, String> {
        let flags = OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX;
        let connection = Connection::open_with_flags(filename, flags).map_err(|x| x.to_string())?;

        // Get some header information.
        let mut get_tag = connection.prepare(
            "SELECT value FROM Tags WHERE key = ?1"
        ).map_err(|x| x.to_string())?;
        let version = get_string_value(&mut get_tag, Self::KEY_VERSION)?;
        if version != Self::VERSION {
            return Err(format!("Unsupported database version: {} (expected {})", version, Self::VERSION));
        }
        let nodes = get_numeric_value(&mut get_tag, Self::KEY_NODES)?;
        let chains = get_numeric_value(&mut get_tag, Self::KEY_CHAINS)?;
        let chain_links = get_numeric_value(&mut get_tag, Self::KEY_CHAIN_LINKS)?;
        let paths = get_numeric_value(&mut get_tag, Self::KEY_PATHS)?;
        let samples = get_numeric_value(&mut get_tag, Self::KEY_SAMPLES)?;
        let haplotypes = get_numeric_value(&mut get_tag, Self::KEY_HAPLOTYPES)?;
        let contigs = get_numeric_value(&mut get_tag, Self::KEY_CONTIGS)?;
        drop(get_tag);

        Ok(GBZBase {
            connection,
            version,
            nodes, chains, chain_links,
            paths, samples, haplotypes, contigs,
        })
    }

    /// Returns the filename of the database or an error if there is no filename.
    pub fn filename(&self) -> Option<&str> {
        self.connection.path()
    }

    /// Returns the size of the database file in a human-readable format.
    pub fn file_size(&self) -> Option<String> {
        let filename = self.filename()?;
        utils::file_size(filename)
    }

    /// Returns the version of the database.
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Returns the number of nodes in the graph.
    pub fn nodes(&self) -> usize {
        self.nodes
    }

    /// Returns the number of top-level chains with stored links between boundary nodes.
    pub fn chains(&self) -> usize {
        self.chains
    }

    /// Returns the total number of links in the top-level chains.
    pub fn chain_links(&self) -> usize {
        self.chain_links
    }

    /// Returns the number of paths in the graph.
    pub fn paths(&self) -> usize {
        self.paths
    }

    /// Returns the number of samples in path metadata.
    pub fn samples(&self) -> usize {
        self.samples
    }

    /// Returns the number of haplotypes in path metadata.
    pub fn haplotypes(&self) -> usize {
        self.haplotypes
    }

    /// Returns the number of contigs in path metadata.
    pub fn contigs(&self) -> usize {
        self.contigs
    }
}

//-----------------------------------------------------------------------------

/// Creating the database.
impl GBZBase {
    /// Creates a new database from the input files.
    ///
    /// # Arguments
    ///
    /// * `gbz_file`: Name of the file containing the GBZ graph.
    /// * `chains_file`: Name of the file containing top-level chains.
    /// * `db_file`: Name of the database file to be created.
    ///
    /// # Errors
    ///
    /// Returns an error if the database already exists or if the GBZ graph does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create_from_files(gbz_file: &Path, chains_file: Option<&Path>, db_file: &Path) -> Result<(), String> {
        eprintln!("Loading GBZ graph {}", gbz_file.display());
        let graph: GBZ = serialize::load_from(gbz_file).map_err(|x| x.to_string())?;
        let chains = if let Some(filename) = chains_file {
            eprintln!("Loading top-level chain file {}", filename.display());
            Chains::load_from(filename)?
        } else {
            Chains::new()
        };
        Self::create(&graph, &chains, db_file)
    }

    // Sanity checks for the GBZ graph. We do not want to handle graphs without sufficient metadata.
    fn sanity_checks(graph: &GBZ) -> Result<(), String> {
        let metadata = graph.metadata().ok_or(
            String::from("The graph does not contain metadata")
        )?;

        if !metadata.has_path_names() {
            return Err("The metadata does not contain path names".to_string());
        }
        if !metadata.has_sample_names() {
            return Err("The metadata does not contain sample names".to_string());
        }
        if !metadata.has_contig_names() {
            return Err("The metadata does not contain contig names".to_string());
        }

        Ok(())
    }

    /// Creates a new database from the given inputs.
    ///
    /// # Arguments
    ///
    /// * `graph`: GBZ graph.
    /// * `chains`: Top-level chains.
    /// * `filename`: Name of the database file to be created.
    ///
    /// # Errors
    ///
    /// Returns an error if the database already exists or if the GBZ graph does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create<P: AsRef<Path>>(graph: &GBZ, chains: &Chains, filename: P) -> Result<(), String> {
        eprintln!("Creating database {}", filename.as_ref().display());
        if utils::file_exists(&filename) {
            return Err(format!("Database {} already exists", filename.as_ref().display()));
        }
        Self::sanity_checks(graph)?;

        let mut connection = Connection::open(filename).map_err(|x| x.to_string())?;
        Self::insert_tags(graph, chains, &mut connection).map_err(|x| x.to_string())?;
        Self::insert_nodes(graph, chains, &mut connection).map_err(|x| x.to_string())?;
        Self::insert_paths(graph, &mut connection).map_err(|x| x.to_string())?;
        Self::index_reference_paths(graph, &mut connection).map_err(|x| x.to_string())?;
        Ok(())
    }

    fn insert_tags(graph: &GBZ, chains: &Chains, connection: &mut Connection) -> rusqlite::Result<()> {
        eprintln!("Inserting header and tags");

        // Create the tags table.
        connection.execute(
            "CREATE TABLE Tags (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            ) STRICT",
            (),
        )?;

        // Insert header and tags.
        let mut inserted = 0;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Tags(key, value) VALUES (?1, ?2)"
            )?;

            // Header.
            let metadata = graph.metadata().unwrap();
            insert.execute((Self::KEY_VERSION, Self::VERSION))?;
            insert.execute((Self::KEY_NODES, graph.nodes()))?;
            insert.execute((Self::KEY_CHAINS, chains.len()))?;
            insert.execute((Self::KEY_CHAIN_LINKS, chains.links()))?;
            insert.execute((Self::KEY_PATHS, metadata.paths()))?;
            insert.execute((Self::KEY_SAMPLES, metadata.samples()))?;
            insert.execute((Self::KEY_HAPLOTYPES, metadata.haplotypes()))?;
            insert.execute((Self::KEY_CONTIGS, metadata.contigs()))?;
            inserted += 6;

            // GBWT tags.
            let index: &GBWT = graph.as_ref();
            for (key, value) in index.tags().iter() {
                let key = format!("{}{}", Self::KEY_GBWT, key);
                insert.execute((key, value))?;
                inserted += 1;
            }

            // GBZ tags.
            for (key, value) in graph.tags().iter() {
                let key = format!("{}{}", Self::KEY_GBZ, key);
                insert.execute((key, value))?;
                inserted += 1;
            }
        }
        transaction.commit()?;

        eprintln!("Inserted {} key-value pairs", inserted);
        Ok(())
    }

    fn insert_nodes(graph: &GBZ, chains: &Chains, connection: &mut Connection) -> rusqlite::Result<()> {
        eprintln!("Inserting nodes");

        // Create the nodes table.
        connection.execute(
            "CREATE TABLE Nodes (
                handle INTEGER PRIMARY KEY,
                edges BLOB NOT NULL,
                bwt BLOB NOT NULL,
                sequence BLOB NOT NULL,
                next INTEGER
            ) STRICT",
            (),
        )?;

        // Insert the nodes.
        let mut inserted = 0;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Nodes(handle, edges, bwt, sequence, next) VALUES (?1, ?2, ?3, ?4, ?5)"
            )?;
            let index: &GBWT = graph.as_ref();
            let bwt: &BWT = index.as_ref();
            for node_id in graph.node_iter() {
                // Forward orientation.
                let forward_id = support::encode_node(node_id, Orientation::Forward);
                let record_id = index.node_to_record(forward_id);
                let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
                let sequence = graph.sequence(node_id).unwrap();
                let encoded_sequence = utils::encode_sequence(sequence);
                let next: Option<usize> = chains.next(forward_id);
                insert.execute((forward_id, edge_bytes, bwt_bytes, encoded_sequence, next))?;
                inserted += 1;
        
                // Reverse orientation.
                let reverse_id = support::encode_node(node_id, Orientation::Reverse);
                let record_id = index.node_to_record(reverse_id);
                let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
                let sequence = support::reverse_complement(sequence);
                let encoded_sequence = utils::encode_sequence(&sequence);
                let next: Option<usize> = chains.next(reverse_id);
                insert.execute((reverse_id, edge_bytes, bwt_bytes, encoded_sequence, next))?;
                inserted += 1;
            }
        }
        transaction.commit()?;

        eprintln!("Inserted {} node records", inserted);
        Ok(())
    }

    fn insert_paths(graph: &GBZ, connection: &mut Connection) -> rusqlite::Result<()> {
        eprintln!("Inserting paths");

        // Create the paths table.
        connection.execute(
            "CREATE TABLE Paths (
                handle INTEGER PRIMARY KEY,
                fw_node INTEGER NOT NULL,
                fw_offset INTEGER NOT NULL,
                rev_node INTEGER NOT NULL,
                rev_offset INTEGER NOT NULL,
                sample TEXT NOT NULL,
                contig TEXT NOT NULL,
                haplotype INTEGER NOT NULL,
                fragment INTEGER NOT NULL,
                is_indexed INTEGER NOT NULL
            ) STRICT",
            (),
        )?;

        // Insert path starts and metadata.
        let mut inserted = 0;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO
                    Paths(handle, fw_node, fw_offset, rev_node, rev_offset, sample, contig, haplotype, fragment, is_indexed)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, FALSE)"
            )?;
            let index: &GBWT = graph.as_ref();
            let metadata = graph.metadata().unwrap();
            for handle in 0..metadata.paths() {
                let fw_start = path_start(index, handle, Orientation::Forward);
                let rev_start = path_start(index, handle, Orientation::Reverse);
                let name = FullPathName::from_metadata(metadata, handle).unwrap();
                insert.execute((
                    handle,
                    fw_start.node, fw_start.offset,
                    rev_start.node, rev_start.offset,
                    name.sample, name.contig, name.haplotype, name.fragment
                ))?;
                inserted += 1;
            }
        }
        transaction.commit()?;

        eprintln!("Inserted information on {} paths", inserted);
        Ok(())
    }

    fn index_reference_paths(graph: &GBZ, connection: &mut Connection) -> rusqlite::Result<()> {
        eprintln!("Indexing reference paths");

        // Create the reference path table.
        connection.execute(
            "CREATE TABLE ReferenceIndex (
                path_handle INTEGER NOT NULL,
                path_offset INTEGER NOT NULL,
                node_handle INTEGER NOT NULL,
                node_offset INTEGER NOT NULL,
                PRIMARY KEY (path_handle, path_offset)
            ) STRICT",
            (),
        )?;

        // NOTE: If a reference haplotype is fragmented, the indexed positions will be relative
        // to the start of each fragment.
        let reference_paths = graph.reference_positions(Self::INDEX_INTERVAL, true);
        if reference_paths.is_empty() {
            eprintln!("No reference paths to index");
            return Ok(());
        }

        // Insert the reference positions into the database.
        eprintln!("Inserting indexed positions");
        let transaction = connection.transaction()?;
        {
            let mut set_as_indexed = transaction.prepare(
                "UPDATE Paths SET is_indexed = TRUE WHERE handle = ?1"
            )?;
            let mut insert_position = transaction.prepare(
                "INSERT INTO ReferenceIndex(path_handle, path_offset, node_handle, node_offset)
                VALUES (?1, ?2, ?3, ?4)"
            )?;
            for ref_path in reference_paths.iter() {
                set_as_indexed.execute((ref_path.id,))?;
                for (path_offset, gbwt_position) in ref_path.positions.iter() {
                    insert_position.execute((ref_path.id, path_offset, gbwt_position.node, gbwt_position.offset))?;
                }
            }
        }
        transaction.commit()?;

        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// A database connection to a GAF-base database.
///
/// This structure stores a database connection and some header information.
/// In multi-threaded applications, each thread should have its own connection.
/// A set of alignments overlapping with a subgraph can be extracted using the [`crate::ReadSet`] structure.
///
/// # Examples
///
/// ```
/// use gbz_base::{GAFBase, GAFBaseParams};
/// use gbz_base::utils;
/// use simple_sds::serialize;
///
/// let gaf_file = utils::get_test_data("micb-kir3dl1_HG003.gaf.gz");
/// let gbwt_file = utils::get_test_data("micb-kir3dl1_HG003.gbwt");
/// let db_file = serialize::temp_file_name("gaf-base");
///
/// // Create the database.
/// let params = GAFBaseParams::default();
/// let db = GAFBase::create_from_files(&gaf_file, &gbwt_file, &db_file, &params);
/// assert!(db.is_ok());
///
/// // Now open it and check some statistics.
/// let db = GAFBase::open(&db_file);
/// assert!(db.is_ok());
/// let db = db.unwrap();
/// assert_eq!(db.nodes(), 2291);
/// assert_eq!(db.alignments(), 12439);
/// assert!(!db.bidirectional_gbwt());
///
/// drop(db);
/// let _ = std::fs::remove_file(&db_file);
/// ```
#[derive(Debug)]
pub struct GAFBase {
    pub(crate) connection: Connection,
    version: String,
    nodes: usize,
    alignments: usize,
    blocks: usize,
    bidirectional_gbwt: bool,
    tags: Tags,
}

/// Using the database.
impl GAFBase {
    // Key for database version.
    const KEY_VERSION: &'static str = "version";

    // FIXME: "GAF-base version 3" for release
    /// Current database version.
    pub const VERSION: &'static str = "GAF-base version 3-dev-2";

    // Key for node count.
    const KEY_NODES: &'static str = "nodes";

    // Key for alignment count.
    const KEY_ALIGNMENTS: &'static str = "alignments";

    // Key for bidirectional GBWT flag.
    const KEY_BIDIRECTIONAL_GBWT: &'static str = "bidirectional_gbwt";

    /// Default block size in alignments.
    pub const BLOCK_SIZE: usize = 1000;

    fn get_string_value(tags: &Tags, key: &str) -> String {
        tags.get(key).cloned().unwrap_or_default()
    }

    fn get_numeric_value(tags: &Tags, key: &str) -> Result<usize, String> {
        let value = Self::get_string_value(tags, key);
        value.parse::<usize>().map_err(|x| format!("Invalid numeric value for key {}: {}", key, x))
    }

    fn get_boolean_value(tags: &Tags, key: &str) -> Result<bool, String> {
        let value = Self::get_numeric_value(tags, key)?;
        match value {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(format!("Invalid boolean value for key {}: {}", key, value)),
        }
    }

    /// Opens a connection to the database in the given file.
    ///
    /// Reads the header information and passes through any database errors.
    pub fn open<P: AsRef<Path>>(filename: P) -> Result<Self, String> {
        let flags = OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX;
        let connection = Connection::open_with_flags(filename, flags).map_err(|x| x.to_string())?;

        // Read all tags from the database.
        let mut get_tags = connection.prepare(
            "SELECT key, value FROM Tags"
        ).map_err(|x| x.to_string())?;
        let mut tags = Tags::new();
        let mut rows = get_tags.query(()).map_err(|x| x.to_string())?;
        while let Some(row) = rows.next().map_err(|x| x.to_string())? {
            let key: String = row.get(0).map_err(|x| x.to_string())?;
            let value: String = row.get(1).map_err(|x| x.to_string())?;
            tags.insert(&key, &value);
        }
        drop(rows);
        drop(get_tags);

        let version = Self::get_string_value(&tags, Self::KEY_VERSION);
        if version != Self::VERSION {
            return Err(format!("Unsupported database version: {} (expected {})", version, Self::VERSION));
        }
        let nodes = Self::get_numeric_value(&tags, Self::KEY_NODES)?;
        let alignments = Self::get_numeric_value(&tags, Self::KEY_ALIGNMENTS)?;
        let bidirectional_gbwt = Self::get_boolean_value(&tags, Self::KEY_BIDIRECTIONAL_GBWT)?;

        // Also determine the number of rows in the Alignments table.
        let mut count_rows = connection.prepare(
            "SELECT COUNT(*) FROM Alignments"
        ).map_err(|x| x.to_string())?;
        let blocks = count_rows.query_row((), |row|
            row.get::<_, usize>(0)
        ).map_err(|x| x.to_string())?;
        drop(count_rows);

        Ok(GAFBase {
            connection,
            version,
            nodes, alignments, blocks, bidirectional_gbwt,
            tags,
        })
    }

    /// Returns the filename of the database or an error if there is no filename.
    pub fn filename(&self) -> Option<&str> {
        self.connection.path()
    }

    /// Returns the size of the database file in a human-readable format.
    pub fn file_size(&self) -> Option<String> {
        let filename = self.filename()?;
        utils::file_size(filename)
    }

    /// Returns the version of the database.
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Returns the number of nodes in the graph.
    pub fn nodes(&self) -> usize {
        self.nodes
    }

    /// Returns the number of alignments in the database.
    pub fn alignments(&self) -> usize {
        self.alignments
    }

    /// Returns the number of database rows storing the alignments.
    ///
    /// Each row corresponds to an [`AlignmentBlock`].
    pub fn blocks(&self) -> usize {
        self.blocks
    }

    /// Returns `true` if the paths are stored in a bidirectional GBWT.
    pub fn bidirectional_gbwt(&self) -> bool {
        self.bidirectional_gbwt
    }

    /// Returns all tags stored in the database.
    pub fn tags(&self) -> &Tags {
        &self.tags
    }

    /// Returns the stable graph name (pggname) for the graph used as the reference for the alignments.
    ///
    /// Returns an error if the tags cannot be parsed.
    pub fn graph_name(&self) -> Result<GraphName, String> {
        GraphName::from_tags(self.tags())
    }
}

//-----------------------------------------------------------------------------

// Statistics for the encoded alignments in the database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
 struct AlignmentStats {
    pub alignments: usize,
    pub start_bytes: usize,
    pub name_bytes: usize,
    pub quality_bytes: usize,
    pub difference_bytes: usize,
    pub flag_bytes: usize,
    pub number_bytes: usize,
}

impl AlignmentStats {
    pub fn new() -> Self {
        Self {
            alignments: 0,
            start_bytes: 0,
            name_bytes: 0,
            quality_bytes: 0,
            difference_bytes: 0,
            flag_bytes: 0,
            number_bytes: 0,
        }
    }

    pub fn update(&mut self, block: &AlignmentBlock) {
        self.alignments += block.len();
        self.start_bytes += block.gbwt_starts.len();
        self.name_bytes += block.names.len();
        self.quality_bytes += block.quality_strings.len();
        self.difference_bytes += block.difference_strings.len();
        self.flag_bytes += block.flags.bytes();
        self.number_bytes += block.numbers.len();
    }
}

// FIXME: store sequences, store quality strings, store optional fields
// FIXME: to tag
/// GAF-base construction parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GAFBaseParams {
    /// Number of alignments in a block (database row).
    pub block_size: usize,
}

impl Default for GAFBaseParams {
    fn default() -> Self {
        Self {
            block_size: GAFBase::BLOCK_SIZE,
        }
    }
}

//-----------------------------------------------------------------------------

// TODO: Once we have GBWT construction in the gbz crate, we can build the GBWT in the background
// and insert it into the database after alignments. We can determine GBWT starting positions for
// the alignments in advance. We know the node from the path, and we know the GBWT sequence id
// from alignment id. GBWT path starts come from the endmarker, and they are therefore before all
// other path visits to that node. If a path is the (i+1)-th path starting from GBWT node v,
// the GBWT starting position is (v, i).

// FIXME: Option to create with sequences; without quality strings
// FIXME: Store options as a tag
/// Creating the database.
impl GAFBase {
    /// Creates a new database from the [`GBWT`] index in file `gbwt_file` and stores the database in file `db_file`.
    ///
    /// The GBWT index can be forward-only or bidirectional.
    /// Path `i` in the GBWT index corresponds to line `i` in the GAF file.
    ///
    /// # Arguments
    ///
    /// * `gaf_file`: GAF file storing the alignments. Can be gzip-compressed.
    /// * `gbwt_file`: GBWT file storing the target paths.
    /// * `db_file`: Output database file.
    /// * `params`: Construction parameters.
    ///
    /// # Errors
    ///
    /// Returns an error if the input files do not exist or the database already exists.
    /// Passes through any database errors.
    pub fn create_from_files(
        gaf_file: &Path, gbwt_file: &Path, db_file: &Path,
        params: &GAFBaseParams
    ) -> Result<(), String> {
        eprintln!("Loading GBWT index {}", gbwt_file.display());
        let index: Arc<GBWT> = Arc::new(serialize::load_from(gbwt_file).map_err(|x| x.to_string())?);
        Self::create(gaf_file, index, db_file, params)
    }

    /// Creates a new database in file `filename` from the given [`GBWT`] index.
    ///
    /// The GBWT index can be forward-only or bidirectional.
    /// Path `i` in the GBWT index corresponds to line `i` in the GAF file.
    ///
    /// # Arguments
    ///
    /// * `gaf_file`: GAF file storing the alignments. Can be gzip-compressed.
    /// * `index`: GBWT index storing the target paths.
    /// * `db_file`: Output database file.
    /// * `params`: Construction parameters.
    ///
    /// # Errors
    ///
    /// Returns an error if the GAF file does not exist or if the database already exists.
    /// Passes through any database errors.
    pub fn create<P: AsRef<Path>, Q: AsRef<Path>>(
        gaf_file: P, index: Arc<GBWT>, db_file: Q,
        params: &GAFBaseParams
    ) -> Result<(), String> {
        eprintln!("Creating database {}", db_file.as_ref().display());
        if utils::file_exists(&db_file) {
            return Err(format!("Database {} already exists", db_file.as_ref().display()));
        }

        let mut connection = Connection::open(&db_file).map_err(|x| x.to_string())?;
        let nodes = Self::insert_nodes(&index, &mut connection).map_err(|x| x.to_string())?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));

        // We parse additional tags from GAF headers.
        let mut gaf_file = utils::open_file(gaf_file)?;
        Self::insert_tags(&index, nodes, &mut gaf_file, &mut connection)?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));

        // `insert_alignments` consumes the connection, as it is moved to another thread.
        Self::insert_alignments(index, gaf_file, connection, params)?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));

        Ok(())
    }

    // Returns the number of paths / alignments in the GBWT index.
    fn gbwt_paths(index: &GBWT) -> usize {
        if index.is_bidirectional() { index.sequences() / 2 } else { index.sequences() }
    }

    // We also include tags parsed from GAF headers.
    fn insert_tags(index: &GBWT, nodes: usize, gaf_file: &mut impl BufRead, connection: &mut Connection) -> Result<(), String> {
        eprintln!("Inserting header and tags");

        // Create the tags table.
        connection.execute(
            "CREATE TABLE Tags (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            ) STRICT",
            (),
        ).map_err(|x| x.to_string())?;

        // We currently only care about tags related to stable graph names.
        let header_lines = formats::read_gaf_header_lines(gaf_file).map_err(|x| x.to_string())?;
        let graph_name = GraphName::from_header_lines(&header_lines).unwrap_or_default();
        let mut additional_tags = Tags::new();
        graph_name.set_tags(&mut additional_tags);

        // Insert header and selected tags.
        let mut inserted = 0;
        let transaction = connection.transaction().map_err(|x| x.to_string())?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Tags(key, value) VALUES (?1, ?2)"
            ).map_err(|x| x.to_string())?;

            let alignments = Self::gbwt_paths(index);
            insert.execute((Self::KEY_VERSION, Self::VERSION)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_NODES, nodes)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_ALIGNMENTS, alignments)).map_err(|x| x.to_string())?;
            let bidirectional: usize = if index.is_bidirectional() { 1 } else { 0 };
            insert.execute((Self::KEY_BIDIRECTIONAL_GBWT, bidirectional)).map_err(|x| x.to_string())?;
            inserted += 4;

            for (key, value) in additional_tags.iter() {
                insert.execute((key, value)).map_err(|x| x.to_string())?;
                inserted += 1;
            }
        }
        transaction.commit().map_err(|x| x.to_string())?;

        eprintln!("Inserted {} key-value pairs", inserted);

        Ok(())
    }

    // Returns the number of nodes in the graph.
    fn insert_nodes(index: &GBWT, connection: &mut Connection) -> rusqlite::Result<usize> {
        eprintln!("Inserting nodes");

        // FIXME: sequence NOT NULL; keep empty if no sequence
        // Create the nodes table.
        connection.execute(
            "CREATE TABLE Nodes (
                handle INTEGER PRIMARY KEY,
                edges BLOB NOT NULL,
                bwt BLOB NOT NULL
            ) STRICT",
            (),
        )?;

        // Insert the nodes.
        let mut inserted = 0;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Nodes(handle, edges, bwt) VALUES (?1, ?2, ?3)"
            )?;
            let bwt: &BWT = index.as_ref();
            for record_id in bwt.id_iter() {
                if record_id == gbz::ENDMARKER {
                    continue;
                }
                let handle = index.record_to_node(record_id);
                let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
                insert.execute((handle, edge_bytes, bwt_bytes))?;
                inserted += 1;
            }
        }
        transaction.commit()?;

        eprintln!("Inserted {} node records", inserted);
        Ok(inserted / 2)
    }

    fn insert_alignments(index: Arc<GBWT>, gaf_file: Box<dyn BufRead>, connection: Connection, params: &GAFBaseParams) -> Result<(), String> {
        eprintln!("Inserting alignments");
        let mut gaf_file = gaf_file;
        let mut connection = connection;

        // FIXME: set everything except min_handle, max_handle, read_length NOT NULL
        connection.execute(
            "CREATE TABLE Alignments (
                id INTEGER PRIMARY KEY,
                min_handle INTEGER,
                max_handle INTEGER CHECK (min_handle <= max_handle),
                alignments INTEGER NOT NULL,
                read_length INTEGER,
                gbwt_starts BLOB,
                names BLOB,
                quality_strings BLOB,
                difference_strings BLOB,
                flags BLOB,
                numbers BLOB
            ) STRICT",
            (),
        ).map_err(|x| x.to_string())?;

        // TODO: Use rtree?
        // Create indexes for min/max nodes.
        connection.execute(
            "CREATE INDEX AlignmentNodeInterval
                ON Alignments(min_handle, max_handle)",
            (),
        ).map_err(|x| x.to_string())?;

        // The main thread parses the GAF file and sends blocks of alignments to an encoder thread.
        // That thread sends the encoded blocks to a third thread that inserts them into the database.
        // If something fails within a thread, it sends an error message to the next thread and stops.
        // If a thread receives an error message, it passes it through and stops.
        let (to_encoder, from_parser) = mpsc::sync_channel(4);
        let (to_insert, from_encoder) = mpsc::sync_channel(4);
        let (to_report, from_insert) = mpsc::sync_channel(1);

        // Encoder thread.
        let gbwt_index = index.clone();
        let encoder_thread = thread::spawn(move || {
            let mut alignment_id = 0;
            loop {
                // This can only fail if the sender is disconnected.
                let block: Result<Vec<Alignment>, String> = from_parser.recv().unwrap();
                match block {
                    Ok(block) => {
                        let encoded = AlignmentBlock::new(&block, &gbwt_index, alignment_id);
                        match encoded {
                            Ok(data) => {
                                let _ = to_insert.send(Ok(data));
                            },
                            Err(message) => {
                                let _ = to_insert.send(Err(message));
                                return;
                            }
                        }
                        alignment_id += block.len();
                        if block.is_empty() {
                            // An empty block indicates that we are done.
                            return;
                        }
                    },
                    Err(message) => {
                        let _ = to_insert.send(Err(message));
                        return;
                    },
                }
            }
        });

        // Insertion thread.
        let insert_thread = thread::spawn(move || {
            let transaction = connection.transaction().map_err(|x| x.to_string());
            if let Err(message) = transaction {
                let _ = to_report.send(Err(message));
                return;
            }
            let transaction = transaction.unwrap();

            let insert = transaction.prepare(
                "INSERT INTO
                    Alignments(min_handle, max_handle, alignments, read_length, gbwt_starts, names, quality_strings, difference_strings, flags, numbers)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)"
            ).map_err(|x| x.to_string());
            if let Err(message) = insert {
                let _ = to_report.send(Err(message));
                return;
            }
            let mut insert = insert.unwrap();

            let mut statistics = AlignmentStats::new();
            loop {
                // This can only fail if the sender is disconnected.
                let block: Result<AlignmentBlock, String> = from_encoder.recv().unwrap();
                match block {
                    Ok(block) => {
                        if block.is_empty() {
                            // An empty block indicates that we are done.
                            break;
                        }
                        statistics.update(&block);
                        let result = insert.execute((
                            block.min_handle, block.max_handle, block.alignments, block.read_length,
                            block.gbwt_starts, block.names,
                            block.quality_strings, block.difference_strings,
                            block.flags.as_ref(), block.numbers
                        )).map_err(|x| x.to_string());
                        if let Err(message) = result {
                            let _ = to_report.send(Err(message));
                            return;
                        }
                    },
                    Err(message) => {
                        let _ = to_report.send(Err(message));
                        return;
                    },
                }
            }

            drop(insert);
            let result = transaction.commit().map_err(|x| x.to_string());
            if let Err(message) = result {
                let _ = to_report.send(Err(message));
                return;
            }
            let _ = to_report.send(Ok(statistics));
        });

        // TODO: Maybe we should start a new block when the min node changes and the block is large enough.
        // Main thread that parses the alignments.
        let mut line_num: usize = 0;
        let mut unaligned_block = false;
        let mut block: Vec<Alignment> = Vec::new();
        let mut failed = false;
        loop {
            let mut buf: Vec<u8> = Vec::new();
            let len = gaf_file.read_until(b'\n', &mut buf).map_err(|x| x.to_string());
            match len {
                Ok(len) => {
                    if len == 0 {
                        // End of file.
                        break;
                    }
                },
                Err(message) => {
                    let _ = to_encoder.send(Err(message));
                    failed = true;
                    break;
                },
            };
            if formats::is_gaf_header_line(&buf) {
                // We should not encounter header lines between alignment lines.
                line_num += 1;
                continue;
            }
            let aln = Alignment::from_gaf(&buf).map_err(
                |x| format!("Failed to parse the alignment on line {}: {}", line_num, x)
            );
            match aln {
                Ok(aln) => {
                    if aln.is_unaligned() != unaligned_block {
                        // We have a new block.
                        if !block.is_empty() {
                            let _ = to_encoder.send(Ok(block));
                            block = Vec::new();
                        }
                        unaligned_block = aln.is_unaligned();
                    }
                    block.push(aln);
                    if block.len() >= params.block_size {
                        let _ = to_encoder.send(Ok(block));
                        block = Vec::new();
                    }
                },
                Err(message) => {
                    let _ = to_encoder.send(Err(message));
                    failed = true;
                    break;
                },
            }
            line_num += 1;
        }

        // Send the last block and a termination signal.
        // Then wait for the threads to finish and get the number of inserted alignments.
        if !failed {
            if !block.is_empty() {
                let _ = to_encoder.send(Ok(block));
            }
            let _ = to_encoder.send(Ok(Vec::new()));
        }
        let _ = encoder_thread.join();
        let _ = insert_thread.join();
        let statistics = from_insert.recv().unwrap()?;

        eprintln!("Inserted information on {} alignments", statistics.alignments);
        let expected_alignments = Self::gbwt_paths(&index);
        if statistics.alignments != expected_alignments {
            eprintln!("Warning: Expected {} alignments", expected_alignments);
        }
        eprintln!(
            "Field sizes: gbwt_starts {}, names {}, quality_strings {}, difference_strings {}, flags {}, numbers {}",
            utils::human_readable_size(statistics.start_bytes),
            utils::human_readable_size(statistics.name_bytes),
            utils::human_readable_size(statistics.quality_bytes),
            utils::human_readable_size(statistics.difference_bytes),
            utils::human_readable_size(statistics.flag_bytes),
            utils::human_readable_size(statistics.number_bytes)
        );

        Ok(())
    }
}

//-----------------------------------------------------------------------------

/// A record for an oriented node.
///
/// The record corresponds to one row in table `Nodes`.
/// It stores the information in a GBWT node record ([`Record`]) and the sequence of the node in the correct orientation.
/// The edges and the sequence are decompressed, while the BWT fragment remains in compressed form.
///
/// If the node is a boundary node in a top-level chain, there may also be a link to the next node in the chain.
/// That requires that the link was also present in the source of the record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBZRecord {
    handle: usize,
    edges: Vec<Pos>,
    bwt: Vec<u8>,
    sequence: Vec<u8>,
    next: Option<usize>,
}

impl GBZRecord {
    // TODO: add next from chains?
    /// Creates a new GBZ record from the given GBZ graph and node handle.
    ///
    /// Returns [`None`] if the node does not exist in the graph.
    pub fn from_gbz(graph: &GBZ, handle: usize) -> Option<Self> {
        let node_id = support::node_id(handle);
        if !graph.has_node(node_id) {
            return None;
        }

        let index: &GBWT = graph.as_ref();
        let bwt_records: &BWT = index.as_ref();
        let record_id = index.node_to_record(handle);
        let (edge_bytes, bwt_bytes) = bwt_records.compressed_record(record_id)?;
        let (edges, _) = Record::decompress_edges(edge_bytes)?;
        let bwt = bwt_bytes.to_vec();
        let sequence = graph.sequence(node_id)?;
        let sequence = if support::node_orientation(handle) == Orientation::Forward {
            sequence.to_vec()
        } else {
            support::reverse_complement(sequence)
        };
        let next = None;

        Some(GBZRecord { handle, edges, bwt, sequence, next })
    }

    /// Creates a new GBZ record from the raw parts.
    ///
    /// This is primarily for testing.
    ///
    /// # Safety
    ///
    /// This is probably safe, even if the record breaks some invariants.
    /// However, I do not have the time and the energy to determine the consequences.
    #[doc(hidden)]
    pub unsafe fn from_raw_parts(handle: usize, edges: Vec<Pos>, bwt: Vec<u8>, sequence: Vec<u8>, next: Option<usize>) -> Self {
        GBZRecord { handle, edges, bwt, sequence, next }
    }

    /// Returns a GBWT record based on this record.
    ///
    /// The lifetime of the returned record is tied to this record.
    ///
    /// # Panics
    ///
    /// Will panic if the record would be empty.
    /// This should never happen with a valid database.
    pub fn to_gbwt_record(&self) -> Record<'_> {
        if self.edges.is_empty() {
            panic!("GBZRecord::to_gbwt_record: Empty record");
        }
        unsafe { Record::from_raw_parts(self.handle, self.edges.clone(), &self.bwt) }
    }

    /// Returns the handle of the record.
    ///
    /// The handle is a [`GBWT`] node identifier.
    #[inline]
    pub fn handle(&self) -> usize {
        self.handle
    }

    /// Returns the node identifier of the record.
    #[inline]
    pub fn id(&self) -> usize {
        support::node_id(self.handle)
    }

    /// Returns the orientation of the record.
    #[inline]
    pub fn orientation(&self) -> Orientation {
        support::node_orientation(self.handle)
    }

    /// Returns an iterator over the handles of successor nodes.
    #[inline]
    pub fn successors(&self) -> impl Iterator<Item = usize> + '_ {
        self.edges.iter().filter(|x| x.node != gbz::ENDMARKER).map(|x| x.node)
    }

    /// Returns an iterator over the outgoing edges from the record.
    ///
    /// Each edge is a pair consisting of a destination node handle and a offset in the corresponding GBWT record.
    /// The edges are sorted by destination node.
    /// This iterator does not list the possible edge to [`gbz::ENDMARKER`], as it only exists for technical purposes.
    #[inline]
    pub fn edges(&self) -> impl Iterator<Item = Pos> + '_ {
        self.edges.iter().filter(|x| x.node != gbz::ENDMARKER).copied()
    }

    /// Returns the slice of edges stored in the record.
    ///
    /// Each edge is a pair consisting of a destination node handle and a offset in the corresponding GBWT record.
    /// The edges are sorted by destination node.
    /// This slice does not list the possible edge to [`gbz::ENDMARKER`], as it only exists for technical purposes.
    pub(crate) fn edges_slice(&self) -> &[Pos] {
        if !self.edges.is_empty() && self.edges[0].node == gbz::ENDMARKER {
            return &self.edges[1..];
        }
        &self.edges
    }

    /// Returns the sequence of the record.
    ///
    /// This is the sequence of the node in the correct orientation.
    /// The sequence is a valid [`String`] over the alphabet `ACGTN`.
    #[inline]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    /// Returns the length of the sequence stored in the record.
    #[inline]
    pub fn sequence_len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns the next handle for the next boundary node in the chain, or [`None`] if there is none.
    #[inline]
    pub fn next(&self) -> Option<usize> {
        self.next
    }
}

//-----------------------------------------------------------------------------

/// A record for a path in the graph.
///
/// The record corresponds to one row in table `Paths`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBZPath {
    /// Path handle / identifier of the path in the original graph.
    pub handle: usize,

    /// Starting position of the path in the forward orientation.
    pub fw_start: Pos,

    /// Starting position of the path in the reverse orientation.
    pub rev_start: Pos,

    /// Path name.
    pub name: FullPathName,

    /// Has this path been indexed for random access with [`GraphInterface::indexed_position`]?
    pub is_indexed: bool,
}

impl GBZPath {
    /// Returns a string representation of the path name.
    pub fn name(&self) -> String {
        self.name.to_string()
    }

    /// Creates a new path record from the given GBZ graph and a path identifier.
    pub fn with_id(graph: &GBZ, path_id: usize) -> Option<Self> {
        let metadata = graph.metadata()?;
        let index: &GBWT = graph.as_ref();
        let fw_start = path_start(index, path_id, Orientation::Forward);
        let rev_start = path_start(index, path_id, Orientation::Reverse);
        let name = FullPathName::from_metadata(metadata, path_id)?;
        Some(GBZPath {
            handle: path_id,
            fw_start, rev_start,
            name,
            is_indexed: false,
        })
    }

    /// Creates a new path record from the given GBZ graph and path name.
    ///
    /// The fragment field is assumed to be an offset in the haplotype.
    /// If the haplotype is fragmented, this returns the last fragment starting at or before the given offset.
    pub fn with_name(graph: &GBZ, name: &FullPathName) -> Option<Self> {
        let metadata = graph.metadata()?;
        let path_id = metadata.find_fragment(name)?;
        let index: &GBWT = graph.as_ref();
        let fw_start = path_start(index, path_id, Orientation::Forward);
        let rev_start = path_start(index, path_id, Orientation::Reverse);
        Some(GBZPath {
            handle: path_id,
            fw_start, rev_start,
            name: FullPathName::from_metadata(metadata, path_id)?,
            is_indexed: false,
        })
    }
}

impl AsRef<FullPathName> for GBZPath {
    fn as_ref(&self) -> &FullPathName {
        &self.name
    }
}

//-----------------------------------------------------------------------------

/// Database query interface.
///
/// This structure stores prepared statements for accessing the graph.
///
/// # Examples
///
/// ```
/// use gbz_base::{GBZBase, GraphInterface};
/// use gbz::{Orientation, FullPathName};
/// use gbz::support;
/// use simple_sds::serialize;
/// use std::fs;
///
/// // Create the database.
/// let gbz_file = support::get_test_data("example.gbz");
/// let db_file = serialize::temp_file_name("graph-interface");
/// let result = GBZBase::create_from_files(&gbz_file, None, &db_file);
/// assert!(result.is_ok());
///
/// // Open the database and create a graph interface.
/// let database = GBZBase::open(&db_file).unwrap();
/// let mut interface = GraphInterface::new(&database).unwrap();
///
/// // The example graph does not have a reference samples tag.
/// assert!(interface.get_gbwt_tag("reference_samples").unwrap().is_none());
///
/// // Node 21 with edges to 22 and 23, all in forward orientation.
/// let id = 21;
/// let orientation = Orientation::Forward;
/// let handle = support::encode_node(id, orientation);
/// let record = interface.get_record(handle).unwrap().unwrap();
/// assert_eq!(record.id(), id);
/// assert_eq!(record.orientation(), orientation);
/// let successors: Vec<usize> = record.successors().collect();
/// assert_eq!(
///     successors,
///     vec![
///         support::encode_node(22, Orientation::Forward),
///         support::encode_node(23, Orientation::Forward)
///     ]
/// );
///
/// // Reference path for contig B goes from 21 to 22.
/// let path_name = FullPathName::generic("B");
/// let path = interface.find_path(&path_name).unwrap().unwrap();
/// assert_eq!(path.fw_start.node, handle);
/// let next = record.to_gbwt_record().lf(path.fw_start.offset).unwrap();
/// assert_eq!(next.node, support::encode_node(22, Orientation::Forward));
///
/// // The first indexed position is at the start of the path.
/// assert!(path.is_indexed);
/// let indexed_pos = interface.indexed_position(path.handle, 3).unwrap().unwrap();
/// assert_eq!(indexed_pos, (0, path.fw_start));
///
/// // Clean up.
/// drop(interface);
/// drop(database);
/// fs::remove_file(&db_file).unwrap();
/// ```
#[derive(Debug)]
pub struct GraphInterface<'a> {
    get_tag: Statement<'a>,
    get_tags: Statement<'a>,
    get_record: Statement<'a>,
    get_path: Statement<'a>,
    find_path: Statement<'a>,
    paths_for_sample: Statement<'a>,
    indexed_position: Statement<'a>,
}

impl<'a> GraphInterface<'a> {
    /// Returns a new interface to the given database.
    ///
    /// Passes through any database errors.
    pub fn new(database: &'a GBZBase) -> Result<Self, String> {
        let get_tag = database.connection.prepare(
            "SELECT value FROM Tags WHERE key = ?1"
        ).map_err(|x| x.to_string())?;

        let get_tags = database.connection.prepare(
            "SELECT key, value FROM Tags WHERE key LIKE ?1"
        ).map_err(|x| x.to_string())?;

        let get_record = database.connection.prepare(
            "SELECT edges, bwt, sequence, next FROM Nodes WHERE handle = ?1"
        ).map_err(|x| x.to_string())?;

        let get_path = database.connection.prepare(
            "SELECT * FROM Paths WHERE handle = ?1"
        ).map_err(|x| x.to_string())?;

        let find_path = database.connection.prepare(
            "SELECT * FROM Paths
            WHERE sample = ?1 AND contig = ?2 AND haplotype = ?3 AND fragment <= ?4
            ORDER BY fragment DESC
            LIMIT 1"
        ).map_err(|x| x.to_string())?;

        let paths_for_sample = database.connection.prepare(
            "SELECT * FROM Paths WHERE sample = ?1"
        ).map_err(|x| x.to_string())?;

        let indexed_position = database.connection.prepare(
            "SELECT path_offset, node_handle, node_offset FROM ReferenceIndex
            WHERE path_handle = ?1 AND path_offset <= ?2
            ORDER BY path_offset DESC
            LIMIT 1"
        ).map_err(|x| x.to_string())?;

        Ok(GraphInterface {
            get_tag, get_tags,
            get_record,
            get_path, find_path, paths_for_sample,
            indexed_position,
        })
    }

    /// Returns the value of the [`GBWT`] tag with the given key, or [`None`] if the tag does not exist.
    pub fn get_gbwt_tag(&mut self, key: &str) -> Result<Option<String>, String> {
        let key = format!("{}{}", GBZBase::KEY_GBWT, key);
        self.get_tag.query_row(
            (key,),
            |row| row.get(0)
        ).optional().map_err(|x| x.to_string())
    }

    /// Returns the value of the [`GBZ`] tag with the given key, or [`None`] if the tag does not exist.
    pub fn get_gbz_tag(&mut self, key: &str) -> Result<Option<String>, String> {
        let key = format!("{}{}", GBZBase::KEY_GBZ, key);
        self.get_tag.query_row(
            (key,),
            |row| row.get(0)
        ).optional().map_err(|x| x.to_string())
    }

    // Returns all tags with the given prefix.
    fn get_tags_with_prefix(&mut self, prefix: &str) -> Result<Tags, String> {
        let mut tags = Tags::new();
        let pattern = format!("{}%", prefix);
        let mut rows = self.get_tags.query((pattern,)).map_err(|x| x.to_string())?;
        while let Some(row) = rows.next().map_err(|x| x.to_string())? {
            let key: String = row.get(0).map_err(|x| x.to_string())?;
            let value: String = row.get(1).map_err(|x| x.to_string())?;
            let key = key.trim_start_matches(prefix).to_string();
            tags.insert(&key, &value);
        }
        Ok(tags)
    }

    /// Returns all [`GBWT`] tags.
    pub fn get_gbwt_tags(&mut self) -> Result<Tags, String> {
        self.get_tags_with_prefix(GBZBase::KEY_GBWT)
    }

    /// Returns all [`GBZ`] tags.
    pub fn get_gbz_tags(&mut self) -> Result<Tags, String> {
        self.get_tags_with_prefix(GBZBase::KEY_GBZ)
    }

    /// Returns the stable graph name (pggname) for the graph.
    ///
    /// Passes through any database errors.
    /// Returns an empty name if the corresponding GBZ tags cannot be parsed.
    pub fn graph_name(&mut self) -> Result<GraphName, String> {
        let tags = self.get_gbz_tags()?;
        Ok(GraphName::from_tags(&tags).unwrap_or_default())
    }

    /// Returns the node record for the given handle, or [`None`] if the node does not exist.
    pub fn get_record(&mut self, handle: usize) -> Result<Option<GBZRecord>, String> {
        self.get_record.query_row(
            (handle,),
            |row| {
                let edge_bytes: Vec<u8> = row.get(0)?;
                let (edges, _) = Record::decompress_edges(&edge_bytes).unwrap();
                let bwt: Vec<u8> = row.get(1)?;
                let encoded_sequence: Vec<u8> = row.get(2)?;
                let sequence: Vec<u8> = utils::decode_sequence(&encoded_sequence);
                let next: Option<usize> = row.get(3)?;
                Ok(GBZRecord { handle, edges, bwt, sequence, next })
            }
        ).optional().map_err(|x| x.to_string())
    }

    fn row_to_gbz_path(row: &Row) -> rusqlite::Result<GBZPath> {
        let handle = row.get(0)?;
        let fw_start = Pos::new(row.get(1)?, row.get(2)?);
        let rev_start = Pos::new(row.get(3)?, row.get(4)?);
        let name = FullPathName {
            sample: row.get(5)?,
            contig: row.get(6)?,
            haplotype: row.get(7)?,
            fragment: row.get(8)?,
        };
        let is_indexed = row.get(9)?;
        Ok(GBZPath {
            handle,
            fw_start, rev_start,
            name,
            is_indexed,
        })
    }

    /// Returns the path with the given handle, or [`None`] if the path does not exist.
    pub fn get_path(&mut self, handle: usize) -> Result<Option<GBZPath>, String> {
        self.get_path.query_row((handle,), Self::row_to_gbz_path).optional().map_err(|x| x.to_string())
    }

    /// Returns the path with the given name, or [`None`] if the path does not exist.
    ///
    /// The fragment field is assumed to be an offset in the haplotype.
    /// If the haplotype is fragmented, this returns the last fragment starting at or before the given offset.
    pub fn find_path(&mut self, name: &FullPathName) -> Result<Option<GBZPath>, String> {
        self.find_path.query_row(
            (&name.sample, &name.contig, name.haplotype, name.fragment),
            Self::row_to_gbz_path
        ).optional().map_err(|x| x.to_string())
    }

    /// Returns all paths with the given sample name.
    pub fn paths_for_sample(&mut self, sample_name: &str) -> Result<Vec<GBZPath>, String> {
        let mut result: Vec<GBZPath> = Vec::new();
        let mut rows = self.paths_for_sample.query((sample_name,)).map_err(|x| x.to_string())?;
        while let Some(row) = rows.next().map_err(|x| x.to_string())? {
            let path = Self::row_to_gbz_path(row).map_err(|x| x.to_string())?;
            result.push(path);
        }
        Ok(result)
    }

    // TODO: List of all paths that are indexed. Handle fragmented haplotypes properly.

    /// Returns the last indexed position at or before offset `path_offset` on path `path_handle`.
    ///
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    pub fn indexed_position(&mut self, path_handle: usize, path_offset: usize) -> Result<Option<(usize, Pos)>, String> {
        self.indexed_position.query_row(
            (path_handle, path_offset),
            |row| {
                let path_offset = row.get(0)?;
                let node_handle = row.get(1)?;
                let node_offset = row.get(2)?;
                Ok((path_offset, Pos::new(node_handle, node_offset)))
            }
        ).optional().map_err(|x| x.to_string())
    }
}

//-----------------------------------------------------------------------------

// Helper types and functions for using the databases.

/// A reference to a GBZ-compatible graph.
///
/// Graph operations with a [`GBZ`] graph take an immutable reference to the graph.
/// The corresponding operations with GBZ-base using [`GraphInterface`] take a mutable reference instead.
/// This wrapper can be created on demand to encapsulate a reference to either type of graph.
pub enum GraphReference<'reference, 'graph> {
    /// A [`GBZ`] graph.
    Gbz(&'reference GBZ),
    /// A [`GraphInterface`].
    Db(&'reference mut GraphInterface<'graph>),
}

impl<'reference, 'graph> GraphReference<'reference, 'graph> {
    /// Returns the record for the oriented node corresponding to the given handle.
    ///
    /// # Errors
    ///
    /// Returns an error if the handle does not exist in the graph.
    /// Passes through any errors from the graph implementation.
    pub fn gbz_record(&mut self, handle: usize) -> Result<GBZRecord, String> {
        match self {
            GraphReference::Gbz(gbz) => {
                GBZRecord::from_gbz(gbz, handle).ok_or_else(|| format!("The graph does not contain handle {}", handle))
            },
            GraphReference::Db(db) => {
                db.get_record(handle)?.ok_or_else(|| format!("The graph does not contain handle {}", handle))
            },
        }
    }

    /// Returns the stable graph name (pggname) for the graph.
    ///
    /// # Errors
    ///
    /// Returns an empty object if the corresponding GBZ tags cannot be parsed.
    /// Passes through any errors from the graph implementation.
    pub fn graph_name(&mut self) -> Result<GraphName, String> {
        match self {
            GraphReference::Gbz(gbz) => {
                Ok(GraphName::from_gbz(gbz))
            },
            GraphReference::Db(db) => {
                db.graph_name()
            },
        }
    }
}

/// Returns the starting position of the given path in the given orientation.
///
/// Returns `(gbz::ENDMARKER, 0)` if the path is empty or does not exist.
pub fn path_start(index: &GBWT, path_id: usize, orientation: Orientation) -> Pos {
    index.start(support::encode_path(path_id, orientation)).unwrap_or(Pos::new(gbz::ENDMARKER, 0))
}

/// Type of a potential database file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DatabaseFileType {
    /// The file does not exist.
    Missing,
    /// The file is not a valid SQLite database.
    NotDatabase,
    /// The file is an unknown SQLite database.
    UnknownDatabase,
    /// The file is a known SQLite database with the given version string.
    Version(String),
}

/// Determines the type of the given file, which may be a SQLite database.
pub fn identify_database<P: AsRef<Path>>(filename: P) -> DatabaseFileType {
    let metadata = fs::metadata(&filename);
    if metadata.is_err() {
        return DatabaseFileType::Missing;
    }
    let metadata = metadata.unwrap();
    if !metadata.is_file() {
        return DatabaseFileType::NotDatabase;
    }

    let flags = OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX;
    let connection = Connection::open_with_flags(filename, flags);
    if connection.is_err() {
        return DatabaseFileType::NotDatabase;
    }
    let connection = connection.unwrap();

    let statement = connection.prepare("SELECT value FROM Tags WHERE key = 'version'");
    if statement.is_err() {
        return DatabaseFileType::UnknownDatabase;
    }
    let mut statement = statement.unwrap();

    let version: rusqlite::Result<String> = statement.query_row([], |row| row.get(0));
    if version.is_err() {
        return DatabaseFileType::UnknownDatabase;
    }
    DatabaseFileType::Version(version.unwrap())
}

// Executes the statement, which is expected to return a single string value.
// Then returns the value.
fn get_string_value(statement: &mut Statement, key: &str) -> Result<String, String> {
    let result: rusqlite::Result<String> = statement.query_row(
        (key,),
        |row| row.get(0)
    );
    match result {
        Ok(value) => Ok(value),
        Err(x) => Err(format!("Key not found: {} ({})", key, x)),
    }
}

// Executes the statement, which is expected to return a single string value.
// Then returns the value as an integer.
fn get_numeric_value(statement: &mut Statement, key: &str) -> Result<usize, String> {
    let value = get_string_value(statement, key)?;
    value.parse::<usize>().map_err(|x| x.to_string())
}

//-----------------------------------------------------------------------------

