//! Utility functions and structures.

use std::collections::HashMap;
use std::fs::File;
use std::ops::{Range, RangeInclusive};
use std::path::{Path, PathBuf};
use std::io::{self, BufRead, BufReader};

use flate2::read::MultiGzDecoder;

use gbz::{GBWT, Pos, ENDMARKER};
use pggname::GraphName;
use simple_sds::binaries;

//-----------------------------------------------------------------------------

/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

//-----------------------------------------------------------------------------

/// Returns a human-readable string representation of a size in bytes.
///
/// Reports the size using three decimal places.
pub fn human_readable_size(size: usize) -> String {
    let (size, unit) = binaries::human_readable_size(size);
    format!("{:.3} {}", size, unit)
}

/// Returns a human-readable string representation of the file size in bytes.
///
/// Reports the size using three decimal places.
/// Returns [`None`] if the file does not exist or the size cannot be determined.
/// See [`simple_sds::binaries::file_size`] and [`simple_sds::binaries::human_readable_size`] for further information.
pub fn file_size<P: AsRef<Path>>(filename: P) -> Option<String> {
    let (size, unit) = binaries::file_size(filename)?;
    Some(format!("{:.3} {}", size, unit))
}

/// Prints the peak resident set size to stderr if it can be determined.
///
/// Reports the size using three decimal places.
pub fn report_peak_memory_usage() {
    let peak_memory = binaries::peak_memory_usage();
    if peak_memory.is_err() {
        return;
    }
    let (size, unit) = binaries::human_readable_size(peak_memory.unwrap());
    eprintln!("Peak memory usage: {:.3} {}", size, unit);
}

/// Returns `true` if the reader appears to be gzip-compressed.
///
/// # Errors
///
/// Passes through all I/O errors from the reader.
pub fn is_gzipped<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    let buffer = reader.fill_buf()?;
    let result = buffer.len() >= 2 && buffer[0..2] == [0x1F, 0x8B];
    Ok(result)
}

/// Returns a buffered reader for the file, which may be gzip-compressed.
///
/// Use `-` as the file name to read from standard input.
///
/// # Errors
///
/// Passes through any I/O errors from trying to open and read the file.
pub fn open_file<P: AsRef<Path>>(filename: P) -> Result<Box<dyn BufRead>, String> {
    let mut inner = if filename.as_ref() == Path::new("-") {
        Box::new(BufReader::new(io::stdin())) as Box<dyn BufRead>
    } else {
        let file = File::open(&filename).map_err(|x| format!("Failed to open file {}: {}", filename.as_ref().display(), x))?;
        Box::new(BufReader::new(file)) as Box<dyn BufRead>
    };
    if is_gzipped(&mut inner).map_err(|x| format!("Failed to read file {}: {}", filename.as_ref().display(), x))? {
        let gz_inner = MultiGzDecoder::new(inner);
        Ok(Box::new(BufReader::new(gz_inner)))
    } else {
        Ok(inner)
    }
}

//-----------------------------------------------------------------------------

// Working with `Vec<u8>` buffers.

/// Appends an unsigned integer a string represented as `Vec<u8>`.
pub fn append_usize(buffer: &mut Vec<u8>, value: usize) {
    buffer.extend_from_slice(value.to_string().as_bytes());
}

/// Appends a signed integer a string represented as `Vec<u8>`.
pub fn append_isize(buffer: &mut Vec<u8>, value: isize) {
    buffer.extend_from_slice(value.to_string().as_bytes());
}

//-----------------------------------------------------------------------------

// Sequence encoding and decoding.

// TODO: Precompute the decoding table for a byte.
const DECODE: [u8; 6] = [0, b'A', b'C', b'G', b'T', b'N'];

/// Decodes a single base encoded with [`encode_base`].
///
/// # Panics
///
/// Panics if `encoded > 5`.
#[inline]
pub fn decode_base(encoded: usize) -> u8 {
    DECODE[encoded]
}

/// Decodes a sequence encoded with [`encode_sequence`].
pub fn decode_sequence(encoded: &[u8]) -> Vec<u8> {
    let capacity = if encoded.is_empty() { 0 } else { 3 * encoded.len() };
    let mut result = Vec::with_capacity(capacity);

    for byte in encoded {
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

/// Encodes a single base.
///
/// Use [`decode_base`] to decode.
#[inline]
pub fn encode_base(base: u8) -> usize {
    ENCODE[base as usize] as usize
}

/// Encodes a DNA sequence into a byte array, storing three bases in a byte.
///
/// Values outside `acgtACGT` are encoded as `N`.
/// The last encoded symbol may be a special 0 character in order to preserve the length.
/// This sentinel is not used when the length is a multiple of 3.
/// Use [`decode_sequence`] to decode the sequence.
pub fn encode_sequence(sequence: &[u8]) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(encoded_length(sequence.len()));

    let mut offset = 0;
    while offset + 3 <= sequence.len() {
        let byte = ENCODE[sequence[offset] as usize] +
            6 * ENCODE[sequence[offset + 1] as usize] +
            36 * ENCODE[sequence[offset + 2] as usize];
        result.push(byte);
        offset += 3;
    }
    if sequence.len() - offset == 1 {
        let byte = ENCODE[sequence[offset] as usize];
        result.push(byte);
    } else if sequence.len() - offset == 2 {
        let byte = ENCODE[sequence[offset] as usize] + 6 * ENCODE[sequence[offset + 1] as usize];
        result.push(byte);
    }

    result
}

/// Returns the length of the encoding for a sequence of the given length.
pub fn encoded_length(sequence_length: usize) -> usize {
    sequence_length.div_ceil(3)
}

//-----------------------------------------------------------------------------

/// Returns an error if the given graph is not a valid reference for the given alignments.
///
/// The comparison is based on the provided [`GraphName`] objects.
/// If either graph name is missing, no error is returned.
/// Otherwise the graph name for the alignments must be a subgraph of the reference graph.
pub fn require_valid_reference(alignments: &GraphName, reference: &GraphName) -> Result<(), String> {
    if !alignments.has_name() || !reference.has_name() {
        return Ok(());
    }
    if !alignments.is_subgraph_of(reference) {
        let description = alignments.describe_relationship(reference, "alignments", "reference graph");
        return Err(format!("The graph is not a valid reference for the alignments:\n{}", description));
    }
    Ok(())
}

//-----------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct NodeIdCluster {
    // Inclusive range of node ids in the cluster.
    node_id_range: RangeInclusive<usize>,
    // Range of indices in the original node id array.
    array_range: Range<usize>,
    // Array offset after the largest gap.
    max_gap_offset: Option<usize>,
}

impl NodeIdCluster {
    // Returns a new cluster covering the given range in the node id array.
    // Assumes sorted and deduplicated node ids.
    fn new(node_ids: &[usize], array_range: Range<usize>) -> Option<Self> {
        if node_ids.is_empty() {
            return None;
        }
        if array_range.is_empty() || array_range.end > node_ids.len() {
            return None;
        }

        let first = node_ids[array_range.start];
        let last = node_ids[array_range.end - 1];
        let node_id_range = first..=last;

        let mut max_gap_length = 0;
        let mut max_gap_offset = None;
        for i in (array_range.start + 1)..array_range.end {
            let gap = node_ids[i] - node_ids[i - 1];
            if gap > max_gap_length {
                max_gap_length = gap;
                max_gap_offset = Some(i);
            }
        }

        Some(Self {
            node_id_range,
            array_range,
            max_gap_offset,
        })
    }

    fn max_gap_length(&self, node_ids: &[usize]) -> Option<usize> {
        let offset = self.max_gap_offset?;
        Some(node_ids[offset] - node_ids[offset - 1])
    }

    // Splits the cluster into two at the largest gap, if any.
    // The return values are the cluster before the gap and the cluster after the gap.
    fn split(self, node_ids: &[usize]) -> (Option<Self>, Option<Self>) {
        if self.max_gap_offset.is_none() {
            return (Some(self), None);
        }

        let offset = self.max_gap_offset.unwrap();
        let left = NodeIdCluster::new(node_ids, self.array_range.start..offset);
        let right = NodeIdCluster::new(node_ids, offset..self.array_range.end);

        (left, right)
    }
}

// TODO: If we stick to a constant threshold, we could determine the final clusters in a single pass.
/// Returns a set of closed ranges that cover all node identifiers in the given set.
///
/// Initially there is a single cluster containing all node ids.
/// Each cluster is recursively split at the longest gap between successive identifiers.
/// The recursion stops when the length of the longest gap is at most `threshold`.
/// This can be useful for partitioning a [`crate::Subgraph`] into multiple ranges before querying [`crate::GAFBase`].
///
/// # Examples
///
/// ```
/// use gbz_base::utils;
///
/// let node_ids = vec![1, 2, 4, 6, 30, 31, 35];
/// let threshold = 10;
/// let clusters = utils::cluster_node_ids(node_ids, threshold);
/// assert_eq!(clusters.len(), 2);
/// assert_eq!(clusters[0], 1..=6);
/// assert_eq!(clusters[1], 30..=35);
/// ```
pub fn cluster_node_ids(node_ids: Vec<usize>, threshold: usize) -> Vec<RangeInclusive<usize>> {
    let mut node_ids = node_ids;
    node_ids.sort_unstable();
    node_ids.dedup();

    let mut stack: Vec<NodeIdCluster> = Vec::new();
    let mut result: Vec<RangeInclusive<usize>> = Vec::new();
    let initial = NodeIdCluster::new(&node_ids, 0..node_ids.len());
    if initial.is_none() {
        return result;
    }
    stack.push(initial.unwrap());

    while let Some(curr) = stack.pop() {
        if let Some(len) = curr.max_gap_length(&node_ids) {
            if len > threshold {
                let (left, right) = curr.split(&node_ids);
                if let Some(right) = right {
                    stack.push(right);
                }
                if let Some(left) = left {
                    stack.push(left);
                }
            } else {
                result.push(curr.node_id_range);
            }
        } else {
            result.push(curr.node_id_range);
        }
    }

    result
}

//-----------------------------------------------------------------------------

/// A structure that determines the GBWT starting positions for paths in a graph.
///
/// The starting positions are iterated in order.
/// If a GBWT index is provided, the starting positions are determined from the paths in the index.
/// Otherwise the starting positions for a unidirectional index are computed on the fly.
pub enum PathStartSource<'a> {
    /// A GBWT index of the paths and the next sequence id that has not been iterated.
    Index(&'a GBWT, usize),
    /// A map storing the number of paths starting from each node so far.
    Map(HashMap<usize, usize>),
}

impl<'a> PathStartSource<'a> {
    /// Returns a new source that computes the starting positions on the fly.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the starting position of the next path.
    ///
    /// If the source is a GBWT index, the provided node identifier is ignored.
    /// Returns [`None`] if the path is empty or all paths in the index have been iterated.
    ///
    /// If the source is a map, returns the starting position that would be assigned to the next path starting from the given node.
    /// Returns [`None`] if the path is empty (if the node identifier is [`ENDMARKER`]).
    pub fn next(&mut self, node_id: usize) -> Option<Pos> {
        match self {
            Self::Index(index, seq_id) => {
                let result = index.start(*seq_id);
                *seq_id += if index.is_bidirectional() { 2 } else { 1 };
                result
            }
            Self::Map(map) => {
                if node_id == ENDMARKER {
                    return None;
                }
                let count = map.entry(node_id).or_insert(0);
                let result = Pos::new(node_id, *count);
                *count += 1;
                Some(result)
            }
        }
    }
}

impl<'a> Default for PathStartSource<'a> {
    fn default() -> Self {
        Self::Map(HashMap::new())
    }
}

impl<'a> From<&'a GBWT> for PathStartSource<'a> {
    fn from(index: &'a GBWT) -> Self {
        Self::Index(index, 0)
    }
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use gbz::support::{self, Orientation};
    use simple_sds::serialize;

    #[test]
    fn sequence_encoding() {
        let full_sequence = b"GATTACACACCAGATNNNNNACATTGAACCTTACACAGTCTGAC";
        for i in 0..full_sequence.len() {
            let sequence = &full_sequence[0..i];
            let encoded = encode_sequence(sequence);
            let decoded = decode_sequence(&encoded);
            assert_eq!(decoded, sequence, "Wrong sequence encoding for length {}", i);
        }
    }

    fn test_cluster(node_ids: Vec<usize>, expected: Vec<RangeInclusive<usize>>, gap_threshold: usize, test_case: &str) {
        let clusters = cluster_node_ids(node_ids, gap_threshold);
        assert_eq!(clusters.len(), expected.len(), "Wrong number of clusters for {}", test_case);
        for (i, cluster) in clusters.iter().enumerate() {
            assert_eq!(cluster, &expected[i], "Wrong cluster {} for {}", i, test_case);
        }
    }

    #[test]
    fn cluster_node_ids_test() {
        let node_ids = Vec::new();
        let expected = Vec::new();
        test_cluster(node_ids, expected, 10, "empty");

        let node_ids = vec![5];
        let expected = vec![5..=5];
        test_cluster(node_ids, expected, 10, "single node");

        let node_ids = vec![5, 6, 7, 8, 9];
        let expected = vec![5..=9];
        test_cluster(node_ids, expected, 10, "continuous nodes");

        let node_ids = vec![6, 9, 7, 5, 8];
        let expected = vec![5..=9];
        test_cluster(node_ids, expected, 10, "continuous nodes unsorted");

        let node_ids = vec![5, 7, 9];
        let expected = vec![5..=9];
        test_cluster(node_ids, expected, 10, "equal gaps");

        let node_ids = vec![5, 6, 7, 20, 21, 22];
        let expected = vec![5..=7, 20..=22];
        test_cluster(node_ids, expected, 10, "two clusters");

        let node_ids = vec![1, 50, 52, 53, 63, 64, 200];
        let expected = vec![1..=1, 50..=64, 200..=200];
        test_cluster(node_ids, expected, 10, "one cluster and outliers");

        let node_ids = vec![1, 50, 52, 53, 73, 74, 200];
        let expected = vec![1..=1, 50..=53, 73..=74, 200..=200];
        test_cluster(node_ids, expected, 10, "two clusters and outliers");
    }

    #[test]
    fn path_start_source() {
        let gbwt_files = vec![
            get_test_data("micb-kir3dl1_HG003.gbwt"),
            get_test_data("bidirectional.gbwt"),
            get_test_data("empty.gbwt"),
            support::get_test_data("example.gbwt"),
            support::get_test_data("translation.gbwt"),
            support::get_test_data("with-empty.gbwt"),
        ];

        for gbwt_file in gbwt_files.iter() {
            let index = serialize::load_from(gbwt_file);
            assert!(index.is_ok(), "Failed to load GBWT index from {}: {}", gbwt_file.display(), index.unwrap_err());
            let index: GBWT = index.unwrap();

            let paths = if index.is_bidirectional() {
                index.sequences() / 2
            } else {
                index.sequences()
            };
            let mut index_source = PathStartSource::from(&index);
            let mut map_source = PathStartSource::new();
            for path_id in 0..paths {
                let seq_id = if index.is_bidirectional() { support::encode_path(path_id, Orientation::Forward) } else { path_id };
                let truth = index.start(seq_id);
                let node_id = truth.map(|pos| pos.node).unwrap_or(ENDMARKER);
                let index_pos = index_source.next(0);
                assert_eq!(index_pos, truth, "Wrong path start from index for path {} in {}", path_id, gbwt_file.display());
                if !index.is_bidirectional() {
                    let map_pos = map_source.next(node_id);
                    assert_eq!(map_pos, truth, "Wrong path start from map for path {} in {}", path_id, gbwt_file.display());
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
