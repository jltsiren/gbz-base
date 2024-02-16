use super::*;

use simple_sds::serialize;

//-----------------------------------------------------------------------------

fn gbz_and_path_index(filename: &'static str, interval: usize) -> (GBZ, PathIndex) {
    let gbz_file = support::get_test_data(filename);
    let graph: GBZ = serialize::load_from(gbz_file).unwrap();
    let path_index = PathIndex::new(&graph, interval, false);
    if let Err(err) = path_index {
        panic!("Failed to create path index with interval {}: {}", interval, err);
    }
    (graph, path_index.unwrap())
}

//-----------------------------------------------------------------------------

#[test]
fn path_index() {
    for interval in 0..10 {
        let (graph, path_index) = gbz_and_path_index("example.gbz", interval);
        let metadata = graph.metadata().unwrap();

        // Use the reference positions as the ground truth.
        let reference_paths = graph.reference_positions(interval, false);
        assert_eq!(path_index.path_count(), reference_paths.len(), "Path count mismatch for interval {}", interval);
        for (index_offset, path) in reference_paths.iter().enumerate() {
            assert_eq!(path_index.path_to_offset(path.id), Some(index_offset), "Wrong index offset for path {} with interval {}", path.id, interval);
            assert_eq!(path_index.offset_to_path(index_offset), Some(path.id), "Wrong path for index offset {} with interval {}", index_offset, interval);
            let path_name = FullPathName::from_metadata(&metadata, path.id).unwrap();
            assert_eq!(path_index.find_path(&graph, &path_name), Some(index_offset), "Path not found for name {} with interval {}", path_name, interval);
            assert_eq!(path_index.path_length(index_offset), Some(path.len), "Wrong path length for index offset {} with interval {}", index_offset, interval);

            let mut pos_offset = 0;
            for query_offset in 0..=path.len {
                while pos_offset + 1 < path.positions.len() && path.positions[pos_offset + 1].0 <= query_offset {
                    pos_offset += 1;
                }
                let truth = Some(path.positions[pos_offset]);
                assert_eq!(path_index.indexed_position(index_offset, query_offset), truth, "Wrong indexed position for query offset {} with interval {}", query_offset, interval);
            }
        }

        // Now try things that should not exist.
        assert_eq!(path_index.path_to_offset(graph.paths()), None, "Found an index offset for a nonexistent path with interval {}", interval);
        assert_eq!(path_index.offset_to_path(reference_paths.len()), None, "Found a path for a nonexistent index offset with interval {}", interval);
        let path_name = FullPathName::generic("nonexistent");
        assert_eq!(path_index.find_path(&graph, &path_name), None, "Found a nonexistent path by name with interval {}", interval);
        assert_eq!(path_index.path_length(reference_paths.len()), None, "Found a length for a nonexistent index offset with interval {}", interval);
    }
}

//-----------------------------------------------------------------------------

// FIXME: Subgraph from GBZ

//-----------------------------------------------------------------------------

// FIXME: Subgraph from DB

//-----------------------------------------------------------------------------

// FIXME: GFA / JSON output from Subgraph

//-----------------------------------------------------------------------------
