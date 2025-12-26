//! And index for random access to GBZ paths by sequence offsets.
//!
//! This index extends the functionality of a [`GBZ`] graph to match [`crate::GraphInterface`].

use gbwt::{GBZ, Pos, FullPathName};

use simple_sds::sparse_vector::{SparseVector, SparseBuilder};
use simple_sds::ops::PredSucc;

use std::collections::HashMap;

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

// TODO: Should this be in the `gbwt` crate?
/// An index for random access to reference and generic paths in a GBZ graph.
///
/// Indexed paths are identified by their offsets in the index.
/// The offsets range from 0 to `path_count() - 1`.
///
/// The combination of [`GBZ`] and [`PathIndex`] is functionally similar to [`crate::GraphInterface`] but tens of times faster.
/// An in-memory graph is better for batch operations, where the loading time is a negligible fraction of the total time.
/// The database is better for interactive applications, where the user works with relatively small subgraphs.
/// For a human graph, the database should be faster than the in-memory graph for subgraphs of up to 1 Mbp.
///
/// # Examples
///
/// ```
/// use gbz_base::PathIndex;
/// use gbwt::{GBZ, FullPathName, Pos, REF_SAMPLE};
/// use gbwt::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("example.gbz");
/// let graph: GBZ = serialize::load_from(&filename).unwrap();
///
/// // Create a path index with 3 bp intervals.
/// let path_index = PathIndex::new(&graph, 3, false);
/// assert!(path_index.is_ok());
/// let path_index = path_index.unwrap();
///
/// // We have two components with one generic path in each.
/// assert_eq!(path_index.path_count(), 2);
/// assert_eq!(path_index.path_length(0), Some(5));
/// assert_eq!(path_index.path_length(1), Some(4));
///
/// // Consider the generic path in component A.
/// let path_name = FullPathName::generic("A");
/// let index_offset = path_index.find_path(&graph, &path_name);
/// assert_eq!(index_offset, Some(0));
/// let index_offset = index_offset.unwrap();
///
/// // There should be two indexed positions for the path.
/// let first_sample = path_index.indexed_position(index_offset, 2);
/// assert_eq!(first_sample, Some((0, Pos::new(22, 0))));
/// let second_sample = path_index.indexed_position(index_offset, 5);
/// assert_eq!(second_sample, Some((3, Pos::new(30, 0))));
/// let next_sample = path_index.indexed_position(index_offset, 100);
/// assert_eq!(next_sample, second_sample);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathIndex {
    // Maps path identifiers to index offsets.
    path_to_offset: HashMap<usize, usize>,

    // Maps index offsets to path identifiers.
    offset_to_path: Vec<usize>,

    // Sequence lengths for each path in bp.
    path_lengths: Vec<usize>,

    // Indexed sequence positions for each path.
    sequence_positions: Vec<SparseVector>,

    // GBWT positions corresponding to the indexed sequence positions.
    gbwt_positions: Vec<Vec<Pos>>,
}

//-----------------------------------------------------------------------------

impl PathIndex {
    /// Creates a new path index for the given GBZ graph.
    ///
    /// The index is built for all reference and generic paths.
    ///
    /// # Arguments
    ///
    /// * `graph`: A GBZ graph.
    /// * `interval`: Approximate distance between indexed positions (in bp).
    /// * `verbose`: Print progress information to stderr.
    pub fn new(graph: &GBZ, interval: usize, verbose: bool) -> Result<Self, String> {
        if verbose {
            eprintln!("Building path index");
        }

        // NOTE: For fragmented haplotypes, the reference positions are relative
        // to the start of the fragment.
        let reference_paths = graph.reference_positions(interval, verbose);
        if reference_paths.is_empty() {
            return Err(String::from("No reference paths to index"));
        }

        let mut path_to_offset: HashMap<usize, usize> = HashMap::with_capacity(reference_paths.len());
        let mut offset_to_path: Vec<usize> = vec![0; reference_paths.len()];
        let mut path_lengths: Vec<usize> = Vec::with_capacity(reference_paths.len());
        let mut sequence_positions = Vec::with_capacity(reference_paths.len());
        let mut gbwt_positions = Vec::with_capacity(reference_paths.len());
        for (offset, ref_path) in reference_paths.iter().enumerate() {
            path_to_offset.insert(ref_path.id, offset);
            offset_to_path[offset] = ref_path.id;
            path_lengths.push(ref_path.len);
            let mut sequence = SparseBuilder::new(ref_path.len, ref_path.positions.len())?;
            let mut gbwt = Vec::with_capacity(ref_path.positions.len());
            for (sequence_pos, gbwt_pos) in ref_path.positions.iter() {
                sequence.set(*sequence_pos);
                gbwt.push(*gbwt_pos);
            }
            sequence_positions.push(SparseVector::try_from(sequence)?);
            gbwt_positions.push(gbwt);
        }

        Ok(PathIndex { path_to_offset, offset_to_path, path_lengths, sequence_positions, gbwt_positions })
    }

    /// Returns the number of indexed paths.
    #[inline]
    pub fn path_count(&self) -> usize {
        self.sequence_positions.len()
    }

    /// Returns the index offset for the indexed path with the given identifier.
    ///
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    #[inline]
    pub fn path_to_offset(&self, path_id: usize) -> Option<usize> {
        self.path_to_offset.get(&path_id).cloned()
    }

    /// Returns the path identifier for the indexed path with the given index offset.
    ///
    /// Returns [`None`] if the path does not exist.
    #[inline]
    pub fn offset_to_path(&self, index_offset: usize) -> Option<usize> {
        self.offset_to_path.get(index_offset).cloned()
    }

    /// Returns the index offset for the path with the given metadata.
    ///
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    pub fn find_path(&self, graph: &GBZ, path_name: &FullPathName) -> Option<usize> {
        let metadata = graph.metadata()?;
        let path_id = metadata.find_path(path_name)?;
        self.path_to_offset(path_id)
    }

    /// Returns the length of the indexed path with the given index offset.
    ///
    /// Returns [`None`] if the path does not exist.
    #[inline]
    pub fn path_length(&self, index_offset: usize) -> Option<usize> {
        self.path_lengths.get(index_offset).cloned()
    }

    /// Returns the last indexed position at or before `offset` on the path with name `path_name`.
    ///
    /// The return value consists of a sequence offset and a GBWT position.
    /// Returns [`None`] if the path does not exist or it has not been indexed.
    /// This is similar to [`crate::GraphInterface::find_path`] followed by [`crate::GraphInterface::indexed_position`].
    ///
    /// # Arguments
    ///
    /// * `index_offset`: Offset of the path in the index.
    /// * `query_offset`: Sequence position in the path (in bp).
    pub fn indexed_position(&self, index_offset: usize, query_offset: usize) -> Option<(usize, Pos)> {
        let mut iter = self.sequence_positions[index_offset].predecessor(query_offset);
        if let Some((sample_offset, sequence_offset)) = iter.next() {
            let gbwt_pos = self.gbwt_positions[index_offset][sample_offset];
            Some((sequence_offset, gbwt_pos))
        } else {
            None
        }
    }
}

//-----------------------------------------------------------------------------
