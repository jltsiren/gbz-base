//! GAF file sorting using multiway external memory merge sort.
//!
//! This module provides functionality to sort GAF (Graph Alignment Format) files
//! using an external memory merge sort algorithm, similar to GNU sort.

// FIXME: multithreading, proper tests
// FIXME: polish

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use gbz::{support, Orientation};

use crate::formats;
use crate::utils;

//-----------------------------------------------------------------------------

/// Types of keys that can be derived from GAF records.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KeyType {
    /// Sort by (minimum handle, maximum handle) in the path.
    /// Handles encode both node ID and orientation.
    NodeInterval,
    /// Sort by hash of the value for random shuffling.
    Hash,
}

/// A GAF record with a sorting key.
#[derive(Clone, Debug)]
struct GAFRecord {
    /// Integer key for sorting.
    key: u64,
    /// Original GAF line.
    value: Vec<u8>,
}

impl GAFRecord {
    /// Missing key value. Records without a key are sorted to the end.
    const MISSING_KEY: u64 = u64::MAX;

    /// 0-based field number for the path in a GAF line.
    const PATH_FIELD: usize = 5;

    /// Creates a new record from a line and sets the key.
    fn new(value: Vec<u8>, key_type: KeyType) -> Self {
        let mut record = Self {
            key: Self::MISSING_KEY,
            value,
        };
        record.set_key(key_type);
        record
    }

    /// Sets the key based on the key type.
    fn set_key(&mut self, key_type: KeyType) {
        self.key = match key_type {
            KeyType::NodeInterval => self.extract_node_interval_key(),
            KeyType::Hash => self.extract_hash_key(),
        };
    }

    /// Extracts (min_handle, max_handle) from the path field.
    fn extract_node_interval_key(&self) -> u64 {
        let path = match self.get_field(Self::PATH_FIELD) {
            Some(p) => p,
            None => return Self::MISSING_KEY,
        };

        let mut min_handle: u32 = u32::MAX;
        let mut max_handle: u32 = 0;

        // Parse the path field: e.g., ">1>2>3" or "<5>6<7"
        let mut i = 0;
        while i < path.len() {
            // Parse orientation character
            let orientation = if path[i] == b'>' {
                i += 1;
                Orientation::Forward
            } else if path[i] == b'<' {
                i += 1;
                Orientation::Reverse
            } else {
                return Self::MISSING_KEY;
            };

            // Parse node id
            let start = i;
            while i < path.len() && path[i].is_ascii_digit() {
                i += 1;
            }

            if start < i {
                if let Ok(id_str) = std::str::from_utf8(&path[start..i]) {
                    if let Ok(id) = id_str.parse::<usize>() {
                        let handle = support::encode_node(id, orientation) as u32;
                        min_handle = min_handle.min(handle);
                        max_handle = max_handle.max(handle);
                    } else {
                        return Self::MISSING_KEY;
                    }
                } else {
                    return Self::MISSING_KEY;
                }
            }
        }

        if min_handle == u32::MAX {
            Self::MISSING_KEY
        } else {
            ((min_handle as u64) << 32) | (max_handle as u64)
        }
    }

    /// Extracts a hash of the value for random shuffling.
    fn extract_hash_key(&self) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut hasher = DefaultHasher::new();
        self.value.hash(&mut hasher);
        hasher.finish()
    }

    /// Gets the nth field (0-indexed) from the GAF line.
    fn get_field(&self, field_index: usize) -> Option<&[u8]> {
        self.value.split(|&b| b == b'\t').nth(field_index)
    }

    /// Serializes the record to a writer.
    fn serialize<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.key.to_le_bytes())?;
        let len = self.value.len() as u64;
        writer.write_all(&len.to_le_bytes())?;
        writer.write_all(&self.value)?;
        Ok(())
    }

    /// Writes a GAF line to the output writer.
    fn write_gaf_line<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.value)?; // The line already includes the newline character.
        Ok(())
    }

    /// Deserializes a record from a reader.
    fn deserialize<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut key_bytes = [0u8; 8];
        reader.read_exact(&mut key_bytes)?;
        let key = u64::from_le_bytes(key_bytes);

        let mut len_bytes = [0u8; 8];
        reader.read_exact(&mut len_bytes)?;
        let len = u64::from_le_bytes(len_bytes) as usize;

        let mut value = vec![0u8; len];
        reader.read_exact(&mut value)?;

        Ok(Self { key, value })
    }

    /// Flips the key to reverse the order (for priority queue).
    fn flip_key(&mut self) {
        self.key = u64::MAX - self.key;
    }
}

impl PartialEq for GAFRecord {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl Eq for GAFRecord {}

impl PartialOrd for GAFRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GAFRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

//-----------------------------------------------------------------------------

/// A temporary file for storing sorted GAF records.
struct TempFile {
    path: Option<PathBuf>,
    records: usize,
}

impl TempFile {
    /// Creates a new temporary file.
    fn create() -> io::Result<Self> {
        let path = std::env::temp_dir().join(format!("gaf-sort-{}.zst", uuid()));
        Ok(Self { path: Some(path), records: 0 })
    }

    /// Opens the file for writing.
    fn writer(&self) -> io::Result<BufWriter<zstd::Encoder<'static, File>>> {
        let file = File::create(self.path.as_ref().unwrap())?;
        let encoder = zstd::Encoder::new(file, 3)?; // Compression level 3
        Ok(BufWriter::new(encoder))
    }

    /// Opens the file for reading.
    fn reader(&self) -> io::Result<zstd::Decoder<'static, BufReader<std::fs::File>>> {
        let file = std::fs::File::open(self.path.as_ref().unwrap())?;
        let decoder = zstd::Decoder::new(file)?;
        Ok(decoder)
    }

    /// Returns a new object for the same file and clears this object.
    fn take(&mut self) -> Self {
        let path = self.path.take();
        let records = self.records;
        self.records = 0;
        Self { path, records }
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        if let Some(path) = &self.path {
            let _ = fs::remove_file(path);
        }
    }
}

/// Generates a simple UUID-like string for temporary files.
fn uuid() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    format!("{:x}", now)
}

//-----------------------------------------------------------------------------

/// Parameters for GAF sorting.
#[derive(Clone, Debug)]
pub struct SortParameters {
    /// Key type used for sorting.
    pub key_type: KeyType,
    /// Number of records per file in the initial sort.
    pub records_per_file: usize,
    /// Number of files to merge at once.
    pub files_per_merge: usize,
    /// Buffer size for reading and writing records.
    pub buffer_size: usize,
    /// Use stable sorting.
    pub stable: bool,
    /// Print progress information to stderr.
    pub progress: bool,
}

impl SortParameters {
    /// Default for `records_per_file`.
    pub const DEFAULT_RECORDS_PER_FILE: usize = 1_000_000;
    /// Default for `files_per_merge`.
    pub const DEFAULT_FILES_PER_MERGE: usize = 32;
    /// Default for `buffer_size`.
    pub const DEFAULT_BUFFER_SIZE: usize = 1000;
}

impl Default for SortParameters {
    fn default() -> Self {
        Self {
            key_type: KeyType::NodeInterval,
            records_per_file: Self::DEFAULT_RECORDS_PER_FILE,
            files_per_merge: Self::DEFAULT_FILES_PER_MERGE,
            buffer_size: Self::DEFAULT_BUFFER_SIZE,
            stable: false,
            progress: false,
        }
    }
}

//-----------------------------------------------------------------------------

/// Sorts a GAF file using multiway external memory merge sort.
///
/// The function reads the input GAF file, sorts it in batches, and writes
/// the sorted output. Temporary files are created and cleaned up automatically.
///
/// # Arguments
///
/// * `input_file`: Path to the input GAF file (possibly gzip-compressed). Use "-" for stdin.
/// * `output_file`: Path to the output GAF file. Use "-" for stdout.
/// * `params`: Sorting parameters.
///
/// # Errors
///
/// Returns an error if file I/O fails or if the GAF records are malformed.
///
/// # Examples
///
/// ```no_run
/// use gbz_base::gaf_sort::{sort_gaf, SortParameters};
///
/// let params = SortParameters::default();
/// sort_gaf("input.gaf.gz", "output.gaf", &params).unwrap();
/// ```
pub fn sort_gaf<P: AsRef<Path>, Q: AsRef<Path>>(
    input_file: P,
    output_file: Q,
    params: &SortParameters,
) -> Result<(), String> {
    let start_time = std::time::Instant::now();
    let mut total_records = 0;

    // Open input file and read header lines
    let mut reader = utils::open_file(input_file.as_ref())?;
    let header_lines = formats::read_gaf_header_lines(&mut reader)
        .map_err(|e| format!("Failed to read GAF header: {}", e))?;

    // Check if we can sort directly to output (single batch)
    let mut peek_lines = Vec::new();
    let mut line = Vec::new();
    for _ in 0..params.records_per_file {
        line.clear();
        match reader.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {
                if !line.is_empty() && line[0] != b'@' {
                    peek_lines.push(line.clone());
                }
            }
            Err(e) => return Err(format!("Failed to read input: {}", e)),
        }
    }

    // Check if there's more data
    line.clear();
    let has_more = match reader.read_until(b'\n', &mut line) {
        Ok(0) => false,
        Ok(_) => !line.is_empty(),
        Err(e) => return Err(format!("Failed to read input: {}", e)),
    };

    if !has_more {
        // Single batch - sort directly to output
        if params.progress {
            eprintln!("Sorting directly to the final output");
        }
        total_records = peek_lines.len();
        sort_to_output(peek_lines, &header_lines, output_file.as_ref(), params)?;
        if params.progress {
            let elapsed = start_time.elapsed().as_secs_f64();
            eprintln!("Sorted {} records in {:.2} seconds", total_records, elapsed);
        }
        return Ok(());
    }

    // Multi-batch sort - need temporary files
    let mut temp_files = Vec::new();

    // Process first batch
    total_records += peek_lines.len();
    let temp = sort_to_temp(&peek_lines, params)?;
    temp_files.push(temp);

    // Add the peeked line
    peek_lines.clear();
    peek_lines.push(line.clone());

    // Read and sort remaining batches
    loop {
        line.clear();
        for _ in 1..params.records_per_file {
            line.clear();
            match reader.read_until(b'\n', &mut line) {
                Ok(0) => break,
                Ok(_) => {
                    if !line.is_empty() && line[0] != b'@' {
                        peek_lines.push(line.clone());
                    }
                }
                Err(e) => return Err(format!("Failed to read input: {}", e)),
            }
        }

        if peek_lines.is_empty() {
            break;
        }

        total_records += peek_lines.len();
        let temp = sort_to_temp(&peek_lines, params)?;
        temp_files.push(temp);
        peek_lines.clear();
    }

    if params.progress {
        eprintln!(
            "Initial sort finished with {} records in {} files",
            total_records,
            temp_files.len()
        );
    }

    // Merge rounds - continue until we have files_per_merge or fewer files
    let mut round = 0;
    while temp_files.len() > params.files_per_merge {
        if params.progress {
            eprintln!("Round {}: merging {} files", round, temp_files.len());
        }

        let mut next_files = Vec::new();
        let mut i = 0;
        while i < temp_files.len() {
            let end = (i + params.files_per_merge).min(temp_files.len());
            if end - i == 1 {
                next_files.push(temp_files[i].take());
            } else {
                let merged = merge_files(&temp_files[i..end], params)?;
                next_files.push(merged);
            }

            i = end;
        }

        temp_files = next_files;
        if params.progress {
            eprintln!("Round {} finished with {} files", round, temp_files.len());
        }
        round += 1;
    }

    // Final merge to output
    if params.progress {
        eprintln!("Starting the final merge");
    }

    if !temp_files.is_empty() {
        merge_to_output(&temp_files, &header_lines, output_file.as_ref(), params)?;
    } else {
        // No files - create empty output
        create_empty_output(&header_lines, output_file.as_ref())?;
    }

    if params.progress {
        let elapsed = start_time.elapsed().as_secs_f64();
        eprintln!("Sorted {} records in {:.2} seconds", total_records, elapsed);
    }

    Ok(())
}

//-----------------------------------------------------------------------------

/// Creates a writer for the output file, or stdout if the path is "-".
fn create_output_writer(output_file: &Path) -> Result<Box<dyn Write>, String> {
    if output_file == Path::new("-") {
        Ok(Box::new(BufWriter::new(io::stdout())))
    } else {
        let file = File::create(output_file)
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Sorts lines directly to output file (for single batch).
fn sort_to_output(
    lines: Vec<Vec<u8>>,
    header_lines: &[String],
    output_file: &Path,
    params: &SortParameters,
) -> Result<(), String> {
    let mut records: Vec<GAFRecord> = lines
        .into_iter()
        .map(|line| GAFRecord::new(line, params.key_type))
        .collect();

    if params.stable {
        records.sort();
    } else {
        records.sort_unstable();
    }

    let mut writer = create_output_writer(output_file)?;

    // Write header
    for line in header_lines {
        writeln!(writer, "{}", line)
            .map_err(|e| format!("Failed to write header: {}", e))?;
    }

    // Write sorted records
    for record in records {
        record.write_gaf_line(&mut writer)
            .map_err(|e| format!("Failed to write output: {}", e))?;
    }

    writer.flush()
        .map_err(|e| format!("Failed to flush output: {}", e))?;

    Ok(())
}

/// Sorts lines to a temporary file.
fn sort_to_temp(lines: &[Vec<u8>], params: &SortParameters) -> Result<TempFile, String> {
    let mut records: Vec<GAFRecord> = lines
        .iter()
        .map(|line| GAFRecord::new(line.clone(), params.key_type))
        .collect();

    if params.stable {
        records.sort();
    } else {
        records.sort_unstable();
    }

    let mut temp = TempFile::create()
        .map_err(|e| format!("Failed to create temporary file: {}", e))?;
    temp.records = records.len();

    let mut writer = temp.writer()
        .map_err(|e| format!("Failed to open temporary file for writing: {}", e))?;

    for record in records {
        record.serialize(&mut writer)
            .map_err(|e| format!("Failed to write to temporary file: {}", e))?;
    }

    writer.into_inner()
        .map_err(|e| format!("Failed to finish compression: {}", e))?
        .finish()
        .map_err(|e| format!("Failed to finish compression: {}", e))?;

    Ok(temp)
}

/// Merges multiple temporary files into a new temporary file.
fn merge_files(inputs: &[TempFile], params: &SortParameters) -> Result<TempFile, String> {
    let mut output = TempFile::create()
        .map_err(|e| format!("Failed to create temporary file: {}", e))?;

    // Open all input files
    let mut readers: Vec<_> = inputs
        .iter()
        .map(|temp| temp.reader())
        .collect::<Result<_, _>>()
        .map_err(|e| format!("Failed to open temporary file for reading: {}", e))?;

    let mut buffers: Vec<VecDeque<GAFRecord>> = vec![VecDeque::new(); readers.len()];
    let mut remaining: Vec<usize> = inputs.iter().map(|t| t.records).collect();

    // Helper to read a buffer
    let read_buffer = |reader_idx: usize,
                       readers: &mut Vec<_>,
                       buffers: &mut Vec<VecDeque<GAFRecord>>,
                       remaining: &mut Vec<usize>|
     -> Result<(), String> {
        let count = remaining[reader_idx].min(params.buffer_size);
        if count > 0 {
            buffers[reader_idx].clear();
            for _ in 0..count {
                let mut record = GAFRecord::deserialize(&mut readers[reader_idx])
                    .map_err(|e| format!("Failed to read from temporary file: {}", e))?;
                record.flip_key(); // Flip for priority queue
                buffers[reader_idx].push_back(record);
            }
            remaining[reader_idx] -= count;
        }
        Ok(())
    };

    // Read initial buffers
    for i in 0..readers.len() {
        read_buffer(i, &mut readers, &mut buffers, &mut remaining)?;
    }

    // Priority queue for merge (max heap, but we flipped keys)
    let mut heap = BinaryHeap::new();
    for (i, buffer) in buffers.iter_mut().enumerate() {
        if let Some(record) = buffer.pop_front() {
            heap.push((record, i));
        }
    }

    // Open output
    let mut writer = output.writer()
        .map_err(|e| format!("Failed to open output temporary file for writing: {}", e))?;
    let mut out_buffer = Vec::new();

    // Merge loop
    while let Some((mut record, source)) = heap.pop() {
        record.flip_key(); // Restore original key
        out_buffer.push(record);

        // Write buffer if full
        if out_buffer.len() >= params.buffer_size {
            for rec in out_buffer.drain(..) {
                rec.serialize(&mut writer)
                    .map_err(|e| format!("Failed to write to output: {}", e))?;
                output.records += 1;
            }
        }

        // Refill source buffer if empty
        if buffers[source].is_empty() && remaining[source] > 0 {
            read_buffer(source, &mut readers, &mut buffers, &mut remaining)?;
        }

        // Add next record from source
        if let Some(next) = buffers[source].pop_front() {
            heap.push((next, source));
        }
    }

    // Write remaining output buffer
    for rec in out_buffer {
        rec.serialize(&mut writer)
            .map_err(|e| format!("Failed to write to output: {}", e))?;
        output.records += 1;
    }

    writer.into_inner()
        .map_err(|e| format!("Failed to finish compression: {}", e))?
        .finish()
        .map_err(|e| format!("Failed to finish compression: {}", e))?;

    Ok(output)
}

/// Merges temporary files to the final output file.
fn merge_to_output(
    inputs: &[TempFile],
    header_lines: &[String],
    output_file: &Path,
    params: &SortParameters,
) -> Result<(), String> {
    let mut writer = create_output_writer(output_file)?;

    // Write header
    for line in header_lines {
        writeln!(writer, "{}", line)
            .map_err(|e| format!("Failed to write header: {}", e))?;
    }

    // Open all input files
    let mut readers: Vec<_> = inputs
        .iter()
        .map(|temp| temp.reader())
        .collect::<Result<_, _>>()
        .map_err(|e| format!("Failed to open temporary file for reading: {}", e))?;

    let mut buffers: Vec<VecDeque<GAFRecord>> = vec![VecDeque::new(); readers.len()];
    let mut remaining: Vec<usize> = inputs.iter().map(|t| t.records).collect();

    // Helper to read a buffer
    let read_buffer = |reader_idx: usize,
                       readers: &mut Vec<_>,
                       buffers: &mut Vec<VecDeque<GAFRecord>>,
                       remaining: &mut Vec<usize>|
     -> Result<(), String> {
        let count = remaining[reader_idx].min(params.buffer_size);
        if count > 0 {
            buffers[reader_idx].clear();
            for _ in 0..count {
                let mut record = GAFRecord::deserialize(&mut readers[reader_idx])
                    .map_err(|e| format!("Failed to read from temporary file: {}", e))?;
                record.flip_key();
                buffers[reader_idx].push_back(record);
            }
            remaining[reader_idx] -= count;
        }
        Ok(())
    };

    // Read initial buffers
    for i in 0..readers.len() {
        read_buffer(i, &mut readers, &mut buffers, &mut remaining)?;
    }

    // Priority queue for merge
    let mut heap = BinaryHeap::new();
    for (i, buffer) in buffers.iter_mut().enumerate() {
        if let Some(record) = buffer.pop_front() {
            heap.push((record, i));
        }
    }

    // Merge loop
    let mut out_buffer = Vec::new();
    while let Some((mut record, source)) = heap.pop() {
        record.flip_key();
        out_buffer.push(record);

        // Write buffer if full
        if out_buffer.len() >= params.buffer_size {
            for rec in out_buffer.drain(..) {
                rec.write_gaf_line(&mut writer)
                    .map_err(|e| format!("Failed to write to output: {}", e))?;
            }
        }

        // Refill source buffer if empty
        if buffers[source].is_empty() && remaining[source] > 0 {
            read_buffer(source, &mut readers, &mut buffers, &mut remaining)?;
        }

        // Add next record from source
        if let Some(next) = buffers[source].pop_front() {
            heap.push((next, source));
        }
    }

    // Write remaining output buffer
    for rec in out_buffer {
        rec.write_gaf_line(&mut writer)
            .map_err(|e| format!("Failed to write to output: {}", e))?;
    }

    writer.flush()
        .map_err(|e| format!("Failed to flush output: {}", e))?;

    Ok(())
}

/// Creates an empty output file with just headers.
fn create_empty_output(header_lines: &[String], output_file: &Path) -> Result<(), String> {
    let mut writer = create_output_writer(output_file)?;

    for line in header_lines {
        writeln!(writer, "{}", line)
            .map_err(|e| format!("Failed to write header: {}", e))?;
    }

    writer.flush()
        .map_err(|e| format!("Failed to flush output: {}", e))?;

    Ok(())
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests;
