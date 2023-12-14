//! # GBZ-base: an immutable pangenome graph stored in a SQLite database
//!
//! This is a prototype for storing a GBZ graph in a SQLite database.
//! It is intended for interactive applications that need immediate access to the graph.
//! In such applications, the overhead from loading the GBZ graph into memory can be significant (e.g. 20 seconds for a human graph).
//! As long as the application needs only a fraction of the entire graph (e.g. 1 Mbp context in a human graph), using the database is faster than loading the graph.
//! This assumes that the database is stored on a local SSD.
//!
//! The prototype builds on the [`gbwt`] crate.
//! Once the implementation has stabilized, it will become an optional feature in the [`gbwt`] crate.
//!
//! ### Basic concepts
//!
//! Nodes are accessed by handles, which are [`GBWT`] node identifiers.
//! A handle encodes both the identifier of the node in the underlying graph and its orientation.
//! Each node record corresponds to a row in table `Nodes`, with the handle as its primary key.
//!
//! Paths are accessed by handles, which are path identifiers in the original graph.
//! Each path record corresponds to a row in table `Paths`, with the handle as its primary key.
//! The record contains information for both orientations of the path.
//!
//! Paths can be indexed for random access, which can be useful for e.g. finding a graph region by its reference coordinates.
//! Indexing is based on storing the sequence offset and the GBWT position at the start of a node once every ~1000 bp.
//! The indexed positions are stored in table `ReferenceIndex`.
//! By default, only generic paths (sample name `_gbwt_ref`) and reference paths (sample name listed in GBWT tag `reference_samples`) are indexed.
//! The database can become excessively large if all paths are indexed.

use std::collections::HashSet;
use std::path::Path;

use rusqlite::{Connection, OpenFlags, OptionalExtension, Row, Statement};

use gbwt::{GBWT, GBZ, Orientation, Pos, REF_SAMPLE, REFERENCE_SAMPLES_KEY};
use gbwt::bwt::{BWT, Record};
use gbwt::support;

use simple_sds::serialize;

pub mod gfa;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// A database connection.
///
/// This structure stores a database connection and some header information.
/// In multi-threaded applications, each thread should have its own connection.
/// Graph operations are supported through the [`GraphInterface`] structure.
///
/// # Examples
///
/// ```
/// use gbz_base::GBZBase;
/// use gbwt::support;
/// use simple_sds::serialize;
/// use std::fs;
///
/// // Create the database.
/// let gbz_file = support::get_test_data("example.gbz");
/// let db_file = serialize::temp_file_name("gbz-base");
/// assert!(!GBZBase::exists(&db_file));
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
    // Index positions at the start of a node on reference paths approximately
    // every this many base pairs.
    const INDEX_INTERVAL: usize = 1000;

    // Key for database version.
    const KEY_VERSION: &'static str = "version";

    /// Current database version.
    pub const VERSION: &'static str = "0.2.0";

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
        let version = Self::get_string_value(&mut get_tag, Self::KEY_VERSION)?;
        if version != Self::VERSION {
            return Err(format!("Unsupported database version: {} (expected {})", version, Self::VERSION));
        }
        let nodes = Self::get_numeric_value(&mut get_tag, Self::KEY_NODES)?;
        let samples = Self::get_numeric_value(&mut get_tag, Self::KEY_SAMPLES)?;
        let haplotypes = Self::get_numeric_value(&mut get_tag, Self::KEY_HAPLOTYPES)?;
        let contigs = Self::get_numeric_value(&mut get_tag, Self::KEY_CONTIGS)?;
        let paths = Self::get_numeric_value(&mut get_tag, Self::KEY_PATHS)?;
        drop(get_tag);

        Ok(GBZBase {
            connection,
            version,
            nodes, samples, haplotypes, contigs, paths,
        })
    }

    /// Returns `true` if the database `filename` exists.
    pub fn exists<P: AsRef<Path>>(filename: P) -> bool {
        let flags = OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX;
        let connection = Connection::open_with_flags(filename, flags);
        connection.is_ok()
    }

    /// Returns the filename of the database.
    pub fn filename(&self) -> Result<&str, String> {
        self.connection.path().ok_or("No filename for the database".to_string())
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

    fn get_numeric_value(statement: &mut Statement, key: &str) -> Result<usize, String> {
        let value = Self::get_string_value(statement, key)?;
        value.parse::<usize>().map_err(|x| x.to_string())
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
        if !graph.has_metadata() {
            return Err("The graph does not contain metadata".to_string());
        }
        let metadata = graph.metadata().unwrap();

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
        if Self::exists(&filename) {
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
            )",
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
            )",
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
                let encoded_sequence = encode_sequence(sequence);
                insert.execute((forward_id, edge_bytes, bwt_bytes, encoded_sequence))?;
                inserted += 1;
        
                // Reverse orientation.
                let reverse_id = support::encode_node(node_id, Orientation::Reverse);
                let record_id = index.node_to_record(reverse_id);
                let (edge_bytes, bwt_bytes) = bwt.compressed_record(record_id).unwrap();
                let sequence = support::reverse_complement(sequence);
                let encoded_sequence = encode_sequence(&sequence);
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
            )",
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
            for (handle, name) in metadata.path_iter().enumerate() {
                let fw_start = index.start(support::encode_node(handle, Orientation::Forward)).unwrap();
                let rev_start = index.start(support::encode_node(handle, Orientation::Reverse)).unwrap();
                insert.execute((
                    handle,
                    fw_start.node, fw_start.offset,
                    rev_start.node, rev_start.offset,
                    metadata.sample_name(name.sample()),
                    metadata.contig_name(name.contig()),
                    name.phase(), name.fragment()
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
            )",
            (),
        )?;

        // Determine reference samples.
        eprint!("Reference samples:");
        let metadata = graph.metadata().unwrap();
        let mut ref_samples: HashSet<usize> = HashSet::new();
        if let Some(id) = metadata.sample_id(REF_SAMPLE) {
            ref_samples.insert(id);
            eprint!(" {}", REF_SAMPLE);
        }
        let index: &GBWT = graph.as_ref();
        let tags = index.tags();
        if let Some(samples) = tags.get(REFERENCE_SAMPLES_KEY) {
            for sample_name in samples.split(' ') {
                if let Some(id) = metadata.sample_id(sample_name) {
                    ref_samples.insert(id);
                    eprint!(" {}", sample_name);
                }
            }
        }
        eprintln!();

        // Determine and index reference paths.
        let mut indexed_paths = 0;
        let mut indexed_positions = 0;
        let transaction = connection.transaction()?;
        {
            let mut set_as_indexed = transaction.prepare(
                "UPDATE Paths SET is_indexed = TRUE WHERE handle = ?1"
            )?;
            let mut insert_position = transaction.prepare(
                "INSERT INTO ReferenceIndex(path_handle, path_offset, node_handle, node_offset)
                VALUES (?1, ?2, ?3, ?4)"
            )?;
            for (path_handle, name) in metadata.path_iter().enumerate() {
                if ref_samples.contains(&name.sample()) {
                    set_as_indexed.execute((path_handle,))?;
                    let mut path_offset = 0;
                    let mut next = 0;
                    let sequence_id = support::encode_path(path_handle, Orientation::Forward);
                    let mut pos = index.start(sequence_id);
                    while let Some(p) = pos {
                        if path_offset >= next {
                            insert_position.execute((path_handle, path_offset, p.node, p.offset))?;
                            indexed_positions += 1;
                            next = path_offset + Self::INDEX_INTERVAL;
                        }
                        path_offset += graph.sequence_len(support::node_id(p.node)).unwrap();
                        pos = index.forward(p);
                    }
                    indexed_paths += 1;
                }
            }
        }
        transaction.commit()?;

        eprintln!("Indexed {} positions on {} reference paths", indexed_positions, indexed_paths);
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
            panic!();
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

    /// Sample name.
    pub sample: String,

    /// Contig name.
    pub contig: String,

    /// Haplotype number.
    pub haplotype: usize,

    /// Fragment number or starting offset of the fragment.
    pub fragment: usize,

    /// Has this path been indexed for random access with [`GraphInterface::indexed_position`]?
    pub is_indexed: bool,
}

impl GBZPath {
    /// Returns a string representation of the name of the path.
    ///
    /// The format is `sample#haplotype#contig@fragment`.
    /// If the path is based on a GFA W-line, `fragment` is the starting offset of this path fragment.
    pub fn name(&self) -> String {
        format!("{}#{}#{}@{}", self.sample, self.haplotype, self.contig, self.fragment)
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
/// use gbwt::Orientation;
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
/// let path = interface.find_path("_gbwt_ref", "B", 0, 0).unwrap().unwrap();
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
            "SELECT * FROM Paths WHERE sample = ?1 AND contig = ?2 AND haplotype = ?3 AND fragment = ?4"
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
                let sequence: Vec<u8> = decode_sequence(&encoded_sequence);
                Ok(GBZRecord { handle, edges, bwt, sequence })
            }
        ).optional().map_err(|x| x.to_string())
    }

    fn row_to_gbz_path(row: &Row) -> rusqlite::Result<GBZPath> {
        let handle = row.get(0)?;
        let fw_start = Pos::new(row.get(1)?, row.get(2)?);
        let rev_start = Pos::new(row.get(3)?, row.get(4)?);
        let sample = row.get(5)?;
        let contig = row.get(6)?;
        let haplotype = row.get(7)?;
        let fragment = row.get(8)?;
        let is_indexed = row.get(9)?;
        Ok(GBZPath {
            handle,
            fw_start, rev_start,
            sample, contig, haplotype, fragment,
            is_indexed,
        })
    }

    /// Returns the path with the given handle, or [`None`] if the path does not exist.
    pub fn get_path(&mut self, handle: usize) -> Result<Option<GBZPath>, String> {
        self.get_path.query_row((handle,), Self::row_to_gbz_path).optional().map_err(|x| x.to_string())
    }

    /// Returns the path with the given metadata, or [`None`] if the path does not exist.
    ///
    /// # Arguments
    ///
    /// * `sample`: Sample name.
    /// * `contig`: Contig name.
    /// * `haplotype`: Haplotype number.
    /// * `fragment`: Fragment number or starting offset of the fragment.
    pub fn find_path(&mut self, sample: &str, contig: &str, haplotype: usize, fragment: usize) -> Result<Option<GBZPath>, String> {
        self.find_path.query_row(
            (sample, contig, haplotype, fragment),
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

const DECODE: [u8; 6] = [0, b'A', b'C', b'G', b'T', b'N'];

fn decode_sequence(encoding: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(3 * encoding.len() - 1);

    for byte in encoding {
        let mut value = *byte as usize;
        for _ in 0..3 {
            let decoded = DECODE[value % DECODE.len()];
            if decoded == 0 {
                return result;
            }
            value /= DECODE.len();
            result.push(decoded);
        }
    }

    result
}

const fn generate_encoding() -> [u8; 256] {
    let mut result = [5; 256];
    result[b'a' as usize] = 1; result[b'A' as usize] = 1;
    result[b'c' as usize] = 2; result[b'C' as usize] = 2;
    result[b'g' as usize] = 3; result[b'G' as usize] = 3;
    result[b't' as usize] = 4; result[b'T' as usize] = 4;
    result
}

const ENCODE: [u8; 256] = generate_encoding();

fn encode_sequence(sequence: &[u8]) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(sequence.len() / 3 + 1);

    let mut offset = 0;
    while offset + 3 <= sequence.len() {
        let byte = ENCODE[sequence[offset] as usize] +
            6 * ENCODE[sequence[offset + 1] as usize] +
            36 * ENCODE[sequence[offset + 2] as usize];
        result.push(byte);
        offset += 3;
    }
    let byte = match sequence.len() - offset {
        0 => 0,
        1 => ENCODE[sequence[offset] as usize],
        _ => ENCODE[sequence[offset] as usize] + 6 * ENCODE[sequence[offset + 1] as usize],
    };
    result.push(byte);

    result
}

//-----------------------------------------------------------------------------
