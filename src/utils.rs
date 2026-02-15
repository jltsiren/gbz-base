//! Utility functions and structures.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::ops::{Range, RangeInclusive};
use std::path::{Path, PathBuf};
use std::io::{self, BufRead, BufReader, Read, Error, ErrorKind};

use flate2::read::MultiGzDecoder;

use gbz::{support, Orientation};
use pggname::GraphName;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Vector, Access};
use simple_sds::serialize::Serialize;

//-----------------------------------------------------------------------------

/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

//-----------------------------------------------------------------------------

// Utilities for working with files.

const SIZE_UNITS: [(f64, &str); 6] = [
    (1.0, "B"),
    (1024.0, "KiB"),
    (1024.0 * 1024.0, "MiB"),
    (1024.0 * 1024.0 * 1024.0, "GiB"),
    (1024.0 * 1024.0 * 1024.0 * 1024.0, "TiB"),
    (1024.0 * 1024.0 * 1024.0 * 1024.0 * 1024.0, "PiB"),
];

/// Returns a human-readable representation of the given number of bytes.
pub fn human_readable_size(bytes: usize) -> String {
    let mut unit = 0;
    let value = bytes as f64;
    while unit + 1 < SIZE_UNITS.len() && value >= SIZE_UNITS[unit + 1].0 {
        unit += 1;
    }
    format!("{:.3} {}", value / SIZE_UNITS[unit].0, SIZE_UNITS[unit].1)
}

/// Returns a human-readable size of the file.
pub fn file_size<P: AsRef<Path>>(filename: P) -> Option<String> {
    let metadata = fs::metadata(filename).map_err(|x| x.to_string());
    if metadata.is_err() {
        return None;
    }
    Some(human_readable_size(metadata.unwrap().len() as usize))
}

/// Returns `true` if the file exists.
pub fn file_exists<P: AsRef<Path>>(filename: P) -> bool {
    fs::metadata(filename).is_ok()
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

// TODO: Find chains from a graph:
// 1. Find all weakly connected components.
// 2. For each component (in parallel?):
//    a. Build a tree of biconnected components.
//    b. Find the longest path in the tree, upweighting the nodes on reference paths.
//    c. Use the longest path as a chain.
// TODO: Move to gbwt-rs?
/// A set of top-level chains represented as links between boundary nodes.
///
/// Top-level chains provide a linear high-level structure for each weakly connected component in the graph.
/// A chain is a sequence of nodes and snarls.
/// Boundary nodes bordering the snarls form a sketch of graph topology.
/// Given a pair of boundary nodes, the graph region between them is either a unary path or a snarl.
/// In both cases, no path can leave the region without visiting one of the boundary nodes.
///
/// This representation is based on storing links between successive boundary nodes.
/// Each link is stored twice, once in each orientation.
///
/// # Examples
///
/// ```
/// use gbz_base::Chains;
/// use gbz_base::utils;
/// use gbz::support::{self, Orientation};
///
/// let filename = utils::get_test_data("micb-kir3dl1.chains");
/// let chains = Chains::load_from(&filename);
/// assert!(chains.is_ok());
/// let chains = chains.unwrap();
///
/// assert_eq!(chains.len(), 2);
/// assert_eq!(chains.links(), 925);
/// let handle = support::encode_node(44, Orientation::Forward);
/// assert!(chains.has_handle(handle));
/// let next = support::encode_node(47, Orientation::Forward);
/// assert_eq!(chains.next(handle), Some(next));
/// ```
pub struct Chains {
    chains: usize,
    next: BTreeMap<usize, usize>,
}

impl Chains {
    // Reads the serialized chains representation.
    fn read_data<R: Read>(reader: &mut R) -> io::Result<Vec<IntVector>> {
        let chains = usize::load(reader)?;
        let mut data: Vec<IntVector> = Vec::with_capacity(chains);
        for _ in 0..chains {
            let vec = IntVector::load(reader)?;
            data.push(vec);
        }
        Ok(data)
    }

    // Converts the chains to a bidirectional link map.
    fn link_map(data: Vec<IntVector>) -> io::Result<BTreeMap<usize, usize>> {
        let mut next = BTreeMap::new();
        for chain in data {
            for i in 1..chain.len() {
                let from = chain.get(i - 1) as usize;
                if next.contains_key(&from) {
                    let msg = format!("Duplicate link from {}", from);
                    return Err(Error::new(ErrorKind::InvalidData, msg));
                }
                let to = chain.get(i) as usize;
                next.insert(from, to);

                let rev_from = support::flip_node(to);
                if next.contains_key(&rev_from) {
                    let msg = format!("Duplicate link from {}", rev_from);
                    return Err(Error::new(ErrorKind::InvalidData, msg));
                }
                let rev_to = support::flip_node(from);
                next.insert(rev_from, rev_to);
            }
        }
        Ok(next)
    }

    /// Creates an empty set of chains.
    pub fn new() -> Self {
        Self {
            chains: 0,
            next: BTreeMap::new(),
        }
    }

    /// Reads the chains from a reader in binary format.
    ///
    /// # Errors
    ///
    /// Passes through all deserialization errors.
    /// Returns an error if a handle occurs in multiple chains.
    pub fn deserialize<R: Read>(reader: &mut R) -> io::Result<Self> {
        let data = Self::read_data(reader)?;
        let chains = data.len();
        let next = Self::link_map(data)?;
        Ok(Self { chains, next })
    }

    /// Reads the chains from a binary file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or [`Self::deserialize`] fails.
    pub fn load_from(filename: &Path) -> Result<Self, String> {
        let mut file = File::open(filename).map_err(|x|
            format!("Failed to open chains file {}: {}", filename.display(), x)
        )?;
        Self::deserialize(&mut file).map_err(|x|
            format!("Failed to read chains from file {}: {}", filename.display(), x)
        )
    }

    /// Returns the number of chains.
    pub fn len(&self) -> usize {
        self.chains
    }

    /// Returns `true` if there are no chains.
    pub fn is_empty(&self) -> bool {
        self.chains == 0
    }

    /// Returns the total number of links in the chains.
    pub fn links(&self) -> usize {
        self.next.len() / 2
    }

    /// Returns the successor for the given handle in the chains, or [`None`] if there is no successor.
    pub fn next(&self, handle: usize) -> Option<usize> {
        self.next.get(&handle).copied()
    }

    /// Returns `true` if the given node is a boundary node in one of the chains.
    pub fn has_node(&self, node_id: usize) -> bool {
        let fw_handle = support::encode_node(node_id, Orientation::Forward);
        let rev_handle = support::encode_node(node_id, Orientation::Reverse);
        self.next.contains_key(&fw_handle) || self.next.contains_key(&rev_handle)
    }

    /// Returns `true` if the given handle refers to a boundary node.
    pub fn has_handle(&self, handle: usize) -> bool {
        let rev_handle = support::flip_node(handle);
        self.next.contains_key(&handle) || self.next.contains_key(&rev_handle)
    }

    /// Returns an iterator over the links, ordered by source handle.
    ///
    /// Filter using [`support::encoded_edge_is_canonical`] to visit each link in a single orientation.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.next.iter().map(|(k, v)| (*k, *v))
    }
}

impl Default for Chains {
    fn default() -> Self {
        Self::new()
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::HashSet;

    use gbz::GBZ;
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

    fn load_chains(filename: &PathBuf) -> Vec<IntVector> {
        let file = File::open(&filename);
        assert!(file.is_ok(), "Failed to open chains file {}", filename.display());
        let mut file = file.unwrap();
        let data = Chains::read_data(&mut file);
        assert!(data.is_ok(), "Failed to read chains from {}", filename.display());
        data.unwrap()
    }

    fn check_chains(chains: &Chains, data: &[IntVector]) {
        assert_eq!(chains.len(), data.len(), "Wrong number of chains");
        let mut expected_links = 0;
        for chain in data.iter() {
            if chain.len() > 1 {
                expected_links += chain.len() - 1;
            }
        }
        assert_eq!(chains.links(), expected_links, "Wrong number of links");

        // Handles and nodes should be present.
        for chain in data.iter() {
            for handle in chain.iter().map(|x| x as usize) {
                assert!(chains.has_handle(handle), "Missing handle {}", handle);
                let rev_handle = support::flip_node(handle);
                assert!(chains.has_handle(rev_handle), "Missing reverse handle {}", rev_handle);
                let (node_id, _) = support::decode_node(handle);
                assert!(chains.has_node(node_id), "Missing node {}", node_id);
            }
        }

        // Missing handles and nodes should not be present.
        let max_handle = data.iter().flat_map(|x| x.iter().map(|y| y as usize)).max().unwrap();
        assert!(!chains.has_handle(max_handle + 2), "Unexpected handle {}", max_handle + 1);
        let missing_node = support::decode_node(max_handle).0 + 1;
        assert!(!chains.has_node(missing_node), "Unexpected node {}", missing_node);
    }

    fn check_region(graph: &GBZ, chains: &Chains, from: usize, to: usize) {
        // Active handles. We proceed to their successors but not predecessors.
        let mut active = vec![from, support::flip_node(to)];
        // Visited node identifiers.
        let mut visited = HashSet::new();
        visited.insert(support::decode_node(from).0);
        visited.insert(support::decode_node(to).0);

        while !active.is_empty() {
            let curr = active.pop().unwrap();
            let (node_id, orientation) = support::decode_node(curr);
            for (next_id, next_o) in graph.successors(node_id, orientation).unwrap() {
                if visited.contains(&next_id) {
                    continue;
                }
                let fw_handle = support::encode_node(next_id, next_o);
                let rev_handle = support::flip_node(fw_handle);
                assert!(!chains.has_node(next_id), "Reached boundary node {} ({}, {}) from region {}..{}", fw_handle, next_id, next_o, from, to);
                active.push(fw_handle); active.push(rev_handle);
                visited.insert(next_id);
            }
        }
    }

    #[test]
    fn chains_empty() {
        let chains = Chains::new();
        assert_eq!(chains.len(), 0, "Expected empty chains");
        assert_eq!(chains.links(), 0, "Expected no links");
    }

    #[test]
    fn chains_nonempty() {
        let chains_file = get_test_data("micb-kir3dl1.chains");
        let data = load_chains(&chains_file);
        let chains = Chains::load_from(&chains_file);
        if let Err(msg) = chains {
            panic!("Failed to read chains from {}: {}", chains_file.display(), msg);
        }
        let chains = chains.unwrap();
        check_chains(&chains, &data);

        let graph_file = get_test_data("micb-kir3dl1.gbz");
        let graph: GBZ = serialize::load_from(&graph_file).unwrap();
        for (from, to) in chains.iter() {
            check_region(&graph, &chains, from, to);
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
}

//-----------------------------------------------------------------------------
