//! GBZ-base: A SQLite database storing a GBZ graph.
// FIXME: also GAF-base

use crate::Alignment;
use crate::alignment::{QualityEncoder, TargetPath};
use crate::utils;

use std::path::Path;
use std::io::BufRead;
use std::sync::{mpsc, Arc};
use std::{fs, thread};

use rusqlite::{Connection, OpenFlags, OptionalExtension, Row, Statement};

use gbwt::{GBWT, GBZ, Orientation, Pos, FullPathName};
use gbwt::bwt::{BWT, Record};
use gbwt::support;

use simple_sds::raw_vector::RawVector;
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
/// use gbwt::support;
/// use simple_sds::serialize;
/// use std::fs;
///
/// // Create the database.
/// let gbz_file = support::get_test_data("example.gbz");
/// let db_file = serialize::temp_file_name("gbz-base");
/// assert!(!utils::file_exists(&db_file));
/// let result = GBZBase::create_from_file(&gbz_file, &db_file);
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
    samples: usize,
    haplotypes: usize,
    contigs: usize,
    paths: usize,
}

/// Using the database.
impl GBZBase {
    /// Index positions at the start of a node on reference paths approximately every this many base pairs.
    pub const INDEX_INTERVAL: usize = 1000;

    // Key for database version.
    const KEY_VERSION: &'static str = "version";

    /// Current database version.
    pub const VERSION: &'static str = "GBZ-base v0.3.0";

    // Key for node count.
    const KEY_NODES: &'static str = "nodes";

    // Key for sample count.
    const KEY_SAMPLES: &'static str = "samples";

    // Key for sample count.
    const KEY_HAPLOTYPES: &'static str = "haplotypes";

    // Key for sample count.
    const KEY_CONTIGS: &'static str = "contigs";

    // Key for path count.
    const KEY_PATHS: &'static str = "paths";

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
        let samples = get_numeric_value(&mut get_tag, Self::KEY_SAMPLES)?;
        let haplotypes = get_numeric_value(&mut get_tag, Self::KEY_HAPLOTYPES)?;
        let contigs = get_numeric_value(&mut get_tag, Self::KEY_CONTIGS)?;
        let paths = get_numeric_value(&mut get_tag, Self::KEY_PATHS)?;
        drop(get_tag);

        Ok(GBZBase {
            connection,
            version,
            nodes, samples, haplotypes, contigs, paths,
        })
    }

    /// Returns the filename of the database.
    pub fn filename(&self) -> Result<&str, String> {
        self.connection.path().ok_or("No filename for the database".to_string())
    }

    /// Returns the size of the database file in a human-readable format.
    pub fn file_size(&self) -> Option<String> {
        let filename = self.filename();
        if filename.is_err() {
            return None;
        }
        utils::file_size(filename.unwrap())
    }

    /// Returns the version of the database.
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Returns the number of nodes in the graph.
    pub fn nodes(&self) -> usize {
        self.nodes
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

    /// Returns the number of paths in the graph.
    pub fn paths(&self) -> usize {
        self.paths
    }
}

//-----------------------------------------------------------------------------

/// Creating the database.
impl GBZBase {
    /// Creates a new database from the GBZ graph in file `gbz_file` and stores the database in file `db_file`.
    ///
    /// # Errors
    ///
    /// Returns an error if the database already exists or if the GBZ graph does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create_from_file<P: AsRef<Path>, Q: AsRef<Path>>(gbz_file: P, db_file: Q) -> Result<(), String> {
        eprintln!("Loading GBZ graph {}", gbz_file.as_ref().display());
        let graph: GBZ = serialize::load_from(&gbz_file).map_err(|x| x.to_string())?;
        Self::create(&graph, db_file)
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

    /// Creates a new database in file `filename` from the given GBZ graph.
    ///
    /// # Errors
    ///
    /// Returns an error if the database already exists or if the GBZ graph does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create<P: AsRef<Path>>(graph: &GBZ, filename: P) -> Result<(), String> {
        eprintln!("Creating database {}", filename.as_ref().display());
        if utils::file_exists(&filename) {
            return Err(format!("Database {} already exists", filename.as_ref().display()));
        }
        Self::sanity_checks(graph)?;

        let mut connection = Connection::open(filename).map_err(|x| x.to_string())?;
        Self::insert_tags(graph, &mut connection).map_err(|x| x.to_string())?;
        Self::insert_nodes(graph, &mut connection).map_err(|x| x.to_string())?;
        Self::insert_paths(graph, &mut connection).map_err(|x| x.to_string())?;
        Self::index_reference_paths(graph, &mut connection).map_err(|x| x.to_string())?;
        Ok(())
    }

    fn insert_tags(graph: &GBZ, connection: &mut Connection) -> rusqlite::Result<()> {
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
            insert.execute((Self::KEY_SAMPLES, metadata.samples()))?;
            insert.execute((Self::KEY_HAPLOTYPES, metadata.haplotypes()))?;
            insert.execute((Self::KEY_CONTIGS, metadata.contigs()))?;
            insert.execute((Self::KEY_PATHS, metadata.paths()))?;
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

    fn insert_nodes(graph: &GBZ, connection: &mut Connection) -> rusqlite::Result<()> {
        eprintln!("Inserting nodes");

        // Create the nodes table.
        connection.execute(
            "CREATE TABLE Nodes (
                handle INTEGER PRIMARY KEY,
                edges BLOB NOT NULL,
                bwt BLOB NOT NULL,
                sequence BLOB NOT NULL
            ) STRICT",
            (),
        )?;

        // Insert the nodes.
        let mut inserted = 0;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Nodes(handle, edges, bwt, sequence) VALUES (?1, ?2, ?3, ?4)"
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
                insert.execute((forward_id, edge_bytes, bwt_bytes, encoded_sequence))?;
                inserted += 1;
        
                // Reverse orientation.
                let reverse_id = support::encode_node(node_id, Orientation::Reverse);
                let record_id = index.node_to_record(reverse_id);
                let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
                let sequence = support::reverse_complement(sequence);
                let encoded_sequence = utils::encode_sequence(&sequence);
                insert.execute((reverse_id, edge_bytes, bwt_bytes, encoded_sequence))?;
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

// FIXME document (requirements for GBWT, GAF), examples, test, query interface
// FIXME remember to test with empty paths
/// A database connection to a GAF-base database.
#[derive(Debug)]
pub struct GAFBase {
    connection: Connection,
    version: String,
    nodes: usize,
    alignments: usize,
    sequences: usize,
    prefix: String,
    // Not in use yet
    _quality_encoder: QualityEncoder,
}

/// Using the database.
impl GAFBase {
    // Key for database version.
    const KEY_VERSION: &'static str = "version";

    /// Current database version.
    pub const VERSION: &'static str = "GAF-base v0.1.0";

    // Key for node count.
    const KEY_NODES: &'static str = "nodes";

    // Key for alignment count.
    const KEY_ALIGNMENTS: &'static str = "alignments";

    // Key for sequence count.
    const KEY_SEQUENCES: &'static str = "sequences";

    // Key for the common prefix of sequence names.
    const KEY_PREFIX: &'static str = "name_prefix";

    // Key for the alphabet of base quality values.
    const KEY_QUALITY_VALUES: &'static str = "quality_values";

    // Key for the canonical Huffman code lengths of base quality values.
    const KEY_QUALITY_LENGTHS: &'static str = "quality_lengths";

    // Key for the common prefix of sample names in GBWT tags.
    // Will be stored in the database using `KEY_PREFIX`.
    const KEY_GBWT_PREFIX: &'static str = "sample_prefix";

    // Key for the alphabet of base quality values in GBWT tags.
    // Will be stored in the database using `KEY_QUALITY_VALUES`.
    const KEY_GBWT_QUALITY_VALUES: &'static str = "quality_values";

    // Key for the canonical Huffman code lengths of base quality values in GBWT tags.
    // Will be stored in the database using `KEY_QUALITY_LENGTHS`.
    const KEY_GBWT_QUALITY_LENGTHS: &'static str = "quality_lengths";

    // Parses the canonical Huffman code lengths.
    fn parse_quality_lengths(value: &str) -> Result<Vec<(usize, usize)>, String> {
        let mut result: Vec<(usize, usize)> = Vec::new();
        for field in value.split(',') {
            let code_len = field.parse::<usize>().map_err(|_|
                format!("Invalid Huffman code length: {}", field)
            )?;
            if let Some(prev) = result.last_mut() {
                if code_len == prev.0 {
                    prev.1 += 1;
                    continue;
                }
            }
            result.push((code_len, 1));
        }
        Ok(result)
    }

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
        let alignments = get_numeric_value(&mut get_tag, Self::KEY_ALIGNMENTS)?;
        let sequences = get_numeric_value(&mut get_tag, Self::KEY_SEQUENCES)?;
        let prefix = get_string_value(&mut get_tag, Self::KEY_PREFIX)?;

        let alphabet = get_string_value(&mut get_tag, Self::KEY_QUALITY_VALUES)?;
        let code_lengths = get_string_value(&mut get_tag, Self::KEY_QUALITY_LENGTHS)?;
        let dictionary = Self::parse_quality_lengths(&code_lengths)?;
        let quality_encoder = QualityEncoder::new(alphabet.as_bytes(), &dictionary).ok_or(
            String::from("Could not build a quality string encoder")
        )?;

        drop(get_tag);
        Ok(GAFBase {
            connection,
            version,
            nodes, alignments, sequences,
            prefix,
            _quality_encoder: quality_encoder,
        })
    }

    /// Returns the filename of the database.
    pub fn filename(&self) -> Result<&str, String> {
        self.connection.path().ok_or("No filename for the database".to_string())
    }

    /// Returns the size of the database file in a human-readable format.
    pub fn file_size(&self) -> Option<String> {
        let filename = self.filename();
        if filename.is_err() {
            return None;
        }
        utils::file_size(filename.unwrap())
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

    /// Returns the number of sequences in the database.
    ///
    /// Each sequence corresponds to one or more alignments.
    pub fn sequences(&self) -> usize {
        self.sequences
    }

    /// Returns the longest common prefix of sequence names.
    ///
    /// This prefix is stored only once to save space.
    pub fn prefix(&self) -> &str {
        &self.prefix
    }
}

//-----------------------------------------------------------------------------

// Statistics for the encoded alignments in the database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
 struct AlignmentStats {
    pub alignments: usize,
    pub name_bytes: usize,
    pub number_bytes: usize,
    pub quality_bytes: usize,
    pub difference_bytes: usize,
    pub pair_bytes: usize,
}

/// Creating the database.
impl GAFBase {
    /// Creates a new database from the GBWT index in file `gbwt_file` and stores the database in file `db_file`.
    ///
    /// # Arguments
    ///
    /// * `gaf_file`: GAF file storing the alignments. Can be gzip-compressed.
    /// * `gbwt_file`: GBWT file storing the target paths.
    /// * `db_file`: Output database file.
    ///
    /// # Errors
    ///
    /// Returns an error if the input files do not exist or the database already exists.
    /// Also returns an error if the GBWT index does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create_from_files<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(gaf_file: P, gbwt_file: Q, db_file: R) -> Result<(), String> {
        eprintln!("Loading GBWT index {}", gbwt_file.as_ref().display());
        let index: Arc<GBWT> = Arc::new(serialize::load_from(&gbwt_file).map_err(|x| x.to_string())?);
        Self::create(gaf_file, index, db_file)
    }

    /// Checks if the GBWT index contains sufficient metadata.
    ///
    /// Returns an error if there is no metadata.
    /// Also returns an error if any of the following is missing:
    ///
    /// * Path names.
    /// * Sample names.
    /// * Tag storing the longest common prefix of sample names.
    /// * Tag storing the alphabet of base quality values.
    pub fn check_gbwt_metadata(index: &GBWT) -> Result<(), String> {
        let metadata = index.metadata().ok_or(
            String::from("The GBWT index does not contain metadata")
        )?;
        if !metadata.has_path_names() {
            return Err("The metadata does not contain path names".to_string());
        }
        if !metadata.has_sample_names() {
            return Err("The metadata does not contain sample names".to_string());
        }

        let tags = index.tags();
        if !tags.contains_key(Self::KEY_GBWT_PREFIX) {
            return Err("The GBWT index does not contain a sample name prefix".to_string());
        }
        if !tags.contains_key(Self::KEY_GBWT_QUALITY_VALUES) {
            return Err("The GBWT index does not contain base quality values".to_string());
        }
        if !tags.contains_key(Self::KEY_GBWT_QUALITY_LENGTHS) {
            return Err("The GBWT index does not contain Huffman code lengths for base quality values".to_string());
        }

        Ok(())
    }

    /// Creates a new database in file `filename` from the given GBWT index.
    ///
    /// # Arguments
    ///
    /// * `gaf_file`: GAF file storing the alignments. Can be gzip-compressed.
    /// * `index`: GBWT index storing the target paths.
    /// * `db_file`: Output database file.
    ///
    /// # Errors
    ///
    /// Returns an error if the GAF file does not exist or if the database already exists.
    /// Also returns an error if the GBWT index does not contain sufficient metadata.
    /// Passes through any database errors.
    pub fn create<P: AsRef<Path>, Q: AsRef<Path>>(gaf_file: P, index: Arc<GBWT>, db_file: Q) -> Result<(), String> {
        eprintln!("Creating database {}", db_file.as_ref().display());
        if utils::file_exists(&db_file) {
            return Err(format!("Database {} already exists", db_file.as_ref().display()));
        }
        Self::check_gbwt_metadata(&index)?;

        let mut connection = Connection::open(&db_file).map_err(|x| x.to_string())?;
        let nodes = Self::insert_nodes(&index, &mut connection).map_err(|x| x.to_string())?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));
        let (prefix, quality_encoder) = Self::insert_tags(&index, nodes, &mut connection)?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));

        // `insert_alignments` consumes the connection, as it is moved to another thread.
        let gaf_file = utils::open_file(gaf_file)?;
        Self::insert_alignments(gaf_file, index, connection, prefix.len(), quality_encoder)?;
        eprintln!("Database size: {}", utils::file_size(&db_file).unwrap_or(String::from("unknown")));

        Ok(())
    }

    // Returns (prefix, quality encoder).
    fn insert_tags(index: &GBWT, nodes: usize, connection: &mut Connection) -> Result<(String, QualityEncoder), String> {
        eprintln!("Inserting header and tags");

        // Create the tags table.
        connection.execute(
            "CREATE TABLE Tags (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            ) STRICT",
            (),
        ).map_err(|x| x.to_string())?;

        let tags = index.tags();
        let prefix = tags.get(Self::KEY_GBWT_PREFIX).unwrap();
        let quality_values = tags.get(Self::KEY_GBWT_QUALITY_VALUES).unwrap();
        let quality_lengths = tags.get(Self::KEY_GBWT_QUALITY_LENGTHS).unwrap();
        let dictionary = Self::parse_quality_lengths(quality_lengths)?;
        let quality_encoder = QualityEncoder::new(quality_values.as_bytes(), &dictionary).ok_or(
            String::from("Could not build a quality string encoder")
        )?;

        // Insert header and selected tags.
        let mut inserted = 0;
        let transaction = connection.transaction().map_err(|x| x.to_string())?;
        {
            let mut insert = transaction.prepare(
                "INSERT INTO Tags(key, value) VALUES (?1, ?2)"
            ).map_err(|x| x.to_string())?;

            let metadata = index.metadata().unwrap();
            insert.execute((Self::KEY_VERSION, Self::VERSION)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_NODES, nodes)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_ALIGNMENTS, metadata.paths())).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_SEQUENCES, metadata.samples())).map_err(|x| x.to_string())?;
            inserted += 4;

            insert.execute((Self::KEY_PREFIX, prefix)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_QUALITY_VALUES, quality_values)).map_err(|x| x.to_string())?;
            insert.execute((Self::KEY_QUALITY_LENGTHS, quality_lengths)).map_err(|x| x.to_string())?;
            inserted += 3;
        }
        transaction.commit().map_err(|x| x.to_string())?;

        eprintln!("Inserted {} key-value pairs", inserted);


        Ok((prefix.clone(), quality_encoder))
    }

    // Returns the number of nodes in the graph.
    fn insert_nodes(index: &GBWT, connection: &mut Connection) -> rusqlite::Result<usize> {
        eprintln!("Inserting nodes");

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
                if record_id == gbwt::ENDMARKER {
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

    fn insert_alignments(
        gaf_file: Box<dyn BufRead>, index: Arc<GBWT>, connection: Connection,
        lcp: usize, quality_encoder: QualityEncoder
    ) -> Result<(), String> {
        eprintln!("Building structures for mapping alignments to GBWT paths");
        let metadata = index.metadata().unwrap();
        let paths_by_sample = Alignment::paths_by_sample(metadata)?;
        let mut used_paths = RawVector::with_len(metadata.paths(), false);

        eprintln!("Inserting alignments");
        let mut gaf_file = gaf_file;
        let mut connection = connection;

        // TODO: optional tags
        // The primary key is the 0-based line number in the GAF file.
        // `start_node` is a foreign key to `Nodes`, but we do not enforce that for performance reasons.
        connection.execute(
            "CREATE TABLE Alignments (
                handle INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                start_node INTEGER NOT NULL,
                numbers BLOB NOT NULL,
                quality BLOB,
                difference BLOB,
                pair BLOB
            ) STRICT",
            (),
        ).map_err(|x| x.to_string())?;

        // The main thread parses the GAF file and sends the alignments to another thread that sets the relative information.
        // That thread sends the alignments to a third thread that inserts them into the database.
        // If something fails within a thread, it sends an error message to the next thread and stops.
        // If a thread receives an error message, it passes it through and stops.
        let (to_relative, from_parser) = mpsc::sync_channel(4);
        let (to_insert, from_relative) = mpsc::sync_channel(4);
        let (to_report, from_insert) = mpsc::sync_channel(1);

        // Relative information thread.
        let gbwt_index = index.clone();
        let relative_thread = thread::spawn(move || {
            let mut line_num: usize = 1;
            loop {
                // This can only fail if the sender is disconnected.
                let block: Result<Vec<Alignment>, String> = from_parser.recv().unwrap();
                match block {
                    Ok(mut block) => {
                        if block.is_empty() {
                            // An empty block indicates that we are done.
                            let _ = to_insert.send(Ok(block));
                            return;
                        }
                        for aln in block.iter_mut() {
                            let result = aln.set_relative_information(
                                &gbwt_index, lcp, &paths_by_sample, Some(&mut used_paths)
                            );
                            if let Err(message) = result {
                                let _ = to_insert.send(Err(format!("Failed to match alignment {} with a GBWT path: {}", line_num, message)));
                                return;
                            };
                            line_num += 1;
                        }
                        let _ = to_insert.send(Ok(block));
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
                    Alignments(handle, name, start_node, numbers, quality, difference, pair)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)"
            ).map_err(|x| x.to_string());
            if let Err(message) = insert {
                let _ = to_report.send(Err(message));
                return;
            }
            let mut insert = insert.unwrap();

            let mut handle: usize = 0;
            let mut statistics = AlignmentStats {
                alignments: 0,
                name_bytes: 0,
                number_bytes: 0,
                quality_bytes: 0,
                difference_bytes: 0,
                pair_bytes: 0,
            };
            loop {
                // This can only fail if the sender is disconnected.
                let block: Result<Vec<Alignment>, String> = from_relative.recv().unwrap();
                match block {
                    Ok(block) => {
                        if block.is_empty() {
                            // An empty block indicates that we are done.
                            break;
                        }
                        for aln in block.iter() {
                            let name = &aln.name[lcp..];
                            let start_node = if let TargetPath::StartPosition(start) = aln.path {
                                start.node
                            } else { unreachable!() };
                            let numbers = aln.encode_numbers();
                            let quality = aln.encode_base_quality(&quality_encoder);
                            let difference = aln.encode_difference();
                            let pair = aln.encode_pair(lcp);
                            statistics.alignments += 1;
                            statistics.name_bytes += name.len();
                            statistics.number_bytes += numbers.len();
                            if let Some(quality) = &quality {
                                statistics.quality_bytes += quality.len();
                            }
                            if let Some(difference) = &difference {
                                statistics.difference_bytes += difference.len();
                            }
                            if let Some(pair) = &pair {
                                statistics.pair_bytes += pair.len();
                            }
                            let result = insert.execute((
                                handle, name, start_node,
                                numbers, quality, difference, pair
                            )).map_err(|x| x.to_string());
                            if let Err(message) = result {
                                let _ = to_report.send(Err(message));
                                return;
                            }
                            handle += 1;
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

        // Main thread that parses the alignments.
        const BLOCK_SIZE: usize = 1000;
        let mut line_num: usize = 1;
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
                    let _ = to_relative.send(Err(message));
                    failed = true;
                    break;
                },
            };
            let aln = Alignment::from_gaf(&buf).map_err(
                |x| format!("Failed to parse the alignment on line {}: {}", line_num, x)
            );
            match aln {
                Ok(aln) => {
                    block.push(aln);
                    if block.len() >= BLOCK_SIZE {
                        let _ = to_relative.send(Ok(block));
                        block = Vec::new();
                    }
                },
                Err(message) => {
                    let _ = to_relative.send(Err(message));
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
                let _ = to_relative.send(Ok(block));
            }
            let _ = to_relative.send(Ok(Vec::new()));
        }
        let _ = relative_thread.join();
        let _ = insert_thread.join();
        let statistics = from_insert.recv().unwrap()?;

        eprintln!("Inserted information on {} alignments", statistics.alignments);
        if statistics.alignments != metadata.paths() {
            eprintln!("Warning: Expected {} alignments", metadata.paths());
        }
        eprintln!(
            "Field sizes: name {}, numbers {}, quality {}, difference {}, pair {}",
            utils::human_readable_size(statistics.name_bytes),
            utils::human_readable_size(statistics.number_bytes),
            utils::human_readable_size(statistics.quality_bytes),
            utils::human_readable_size(statistics.difference_bytes),
            utils::human_readable_size(statistics.pair_bytes)
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GBZRecord {
    handle: usize,
    edges: Vec<Pos>,
    bwt: Vec<u8>,
    sequence: Vec<u8>,
}

impl GBZRecord {
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

        Some(GBZRecord { handle, edges, bwt, sequence })
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
    pub unsafe fn from_raw_parts(handle: usize, edges: Vec<Pos>, bwt: Vec<u8>, sequence: Vec<u8>) -> Self {
        GBZRecord { handle, edges, bwt, sequence }
    }

    /// Returns a GBWT record based on this record.
    ///
    /// The lifetime of the returned record is tied to this record.
    ///
    /// # Panics
    ///
    /// Will panic if the record would be empty.
    /// This should never happen with a valid database.
    pub fn to_gbwt_record(&self) -> Record {
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
        self.edges.iter().filter(|x| x.node != gbwt::ENDMARKER).map(|x| x.node)
    }

    /// Returns an iterator over the outgoing edges from the record.
    ///
    /// Each edge is a pair consisting of a destination node handle and a offset in the corresponding GBWT record.
    /// The edges are sorted by destination node.
    /// This iterator does not list the possible edge to [`gbwt::ENDMARKER`], as it only exists for technical purposes.
    #[inline]
    pub fn edges(&self) -> impl Iterator<Item = Pos> + '_ {
        self.edges.iter().filter(|x| x.node != gbwt::ENDMARKER).copied()
    }

    /// Returns the slice of edges stored in the record.
    ///
    /// Each edge is a pair consisting of a destination node handle and a offset in the corresponding GBWT record.
    /// The edges are sorted by destination node.
    /// This slice does not list the possible edge to [`gbwt::ENDMARKER`], as it only exists for technical purposes.
    pub(crate) fn edges_slice(&self) -> &[Pos] {
        if !self.edges.is_empty() && self.edges[0].node == gbwt::ENDMARKER {
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
/// use gbwt::{Orientation, FullPathName};
/// use gbwt::support;
/// use simple_sds::serialize;
/// use std::fs;
///
/// // Create the database.
/// let gbz_file = support::get_test_data("example.gbz");
/// let db_file = serialize::temp_file_name("graph-interface");
/// let result = GBZBase::create_from_file(&gbz_file, &db_file);
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

        let get_record = database.connection.prepare(
            "SELECT edges, bwt, sequence FROM Nodes WHERE handle = ?1"
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
            get_tag,
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
                Ok(GBZRecord { handle, edges, bwt, sequence })
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

// Helper functions for using the databases.

/// Returns the starting position of the given path in the given orientation.
///
/// Returns `(gbwt::ENDMARKER, 0)` if the path is empty or does not exist.
pub fn path_start(index: &GBWT, path_id: usize, orientation: Orientation) -> Pos {
    index.start(support::encode_path(path_id, orientation)).unwrap_or(Pos::new(gbwt::ENDMARKER, 0))
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

