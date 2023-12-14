use crate::{GBZRecord, GBZPath, GraphInterface};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};

use gbwt::{Orientation, Pos, ENDMARKER};

use gbwt::support;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GraphPosition {
    pub node: usize,
    pub orientation: Orientation,
    pub offset: usize,
}

//-----------------------------------------------------------------------------

// Returns the graph position and the GBWT position for the given offset.
pub fn query_position(graph: &mut GraphInterface, path: &GBZPath, query_offset: usize) -> Result<(GraphPosition, Pos), String> {
    let result = graph.indexed_position(path.handle, query_offset)?;
    let (mut path_offset, mut pos) = result.ok_or(format!("Path {} is not indexed", path.name()))?;

    let mut graph_pos: Option<GraphPosition> = None;
    let mut gbwt_pos: Option<Pos> = None;
    loop {
        let handle = pos.node;
        let record = graph.get_record(handle)?;
        let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
        if path_offset + record.sequence_len() > query_offset {
            graph_pos = Some(GraphPosition {
                node: support::node_id(handle),
                orientation: support::node_orientation(handle),
                offset: query_offset - path_offset,
            });
            gbwt_pos = Some(pos);
            break;
        }
        path_offset += record.sequence_len();
        let gbwt_record = record.to_gbwt_record();
        let next = gbwt_record.lf(pos.offset);
        if next.is_none() {
            break;
        }
        pos = next.unwrap();
    }

    let graph_pos = graph_pos.ok_or(
        format!("Path {} does not contain offset {}", path.name(), query_offset)
    )?;
    let gbwt_pos = gbwt_pos.unwrap();
    Ok((graph_pos, gbwt_pos))
}

//-----------------------------------------------------------------------------

pub fn distance_to_end(record: &GBZRecord, orientation: Orientation, offset: usize) -> usize {
    if orientation == support::node_orientation(record.handle()) {
        record.sequence_len() - offset
    } else {
        offset + 1
    }
}

pub fn extract_context(
    graph: &mut GraphInterface,
    from: GraphPosition,
    context: usize
) -> Result<BTreeMap<usize, GBZRecord>, String> {
    // Start graph traversal from the initial node.
    let mut active: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::new(); // (distance, node id)
    active.push(Reverse((0, from.node)));

    // Traverse in both directions.
    let mut selected: BTreeMap<usize, GBZRecord> = BTreeMap::new();
    while !active.is_empty() {
        let (distance, node_id) = active.pop().unwrap().0;
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let handle = support::encode_node(node_id, orientation);
            if selected.contains_key(&handle) {
                continue;
            }
            let record = graph.get_record(handle)?;
            let record = record.ok_or(format!("The graph does not contain handle {}", handle))?;
            let next_distance = if node_id == from.node {
                distance_to_end(&record, from.orientation, from.offset)
            } else {
                distance + record.sequence_len()
            };
            if next_distance <= context {
                for successor in record.successors() {
                    if !selected.contains_key(&successor) {
                        active.push(Reverse((next_distance, support::node_id(successor))));
                    }
                }
            }
            selected.insert(handle, record);
        }
    }

    Ok(selected)
}

//-----------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct PathInfo {
    pub path: Vec<usize>,
    pub len: usize,
    pub weight: Option<usize>,
}

impl PathInfo {
    pub fn new(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: None }
    }

    pub fn weighted(path: Vec<usize>, len: usize) -> Self {
        PathInfo { path, len, weight: Some(1) }
    }

    fn increment_weight(&mut self) {
        if let Some(weight) = self.weight {
            self.weight = Some(weight + 1);
        }
    }
}

fn next_pos(pos: Pos, successors: &BTreeMap<usize, Vec<(Pos, bool)>>) -> Option<Pos> {
    if let Some(v) = successors.get(&pos.node) {
        let (next, _) = v[pos.offset];
        if next.node == ENDMARKER || !successors.contains_key(&next.node) {
            None
        } else {
            Some(next)
        }
    } else {
        None
    }
}

// Extract all paths in the subgraph. The second return value is
// (offset in result, offset on that path) for the handle corresponding to `ref_pos`.
pub fn extract_paths(
    subgraph: &BTreeMap<usize, GBZRecord>,
    ref_pos: Pos
) -> Result<(Vec<PathInfo>, (usize, usize)), String> {
    // Decompress the GBWT node records for the subgraph.
    let mut keys: Vec<usize> = Vec::new();
    let mut successors: BTreeMap<usize, Vec<(Pos, bool)>> = BTreeMap::new();
    for (handle, gbz_record) in subgraph.iter() {
        let gbwt_record = gbz_record.to_gbwt_record();
        let decompressed: Vec<(Pos, bool)> = gbwt_record.decompress().into_iter().map(|x| (x, false)).collect();
        keys.push(*handle);
        successors.insert(*handle, decompressed);
    }

    // Mark the positions that have predecessors in the subgraph.
    for handle in keys.iter() {
        let decompressed = successors.get(handle).unwrap().clone();
        for (pos, _) in decompressed.iter() {
            if let Some(v) = successors.get_mut(&pos.node) {
                v[pos.offset].1 = true;
            }
        }
    }

    // FIXME: Check for infinite loops.
    // Extract all paths and note if one of them passes through `ref_pos`.
    let mut result: Vec<PathInfo> = Vec::new();
    let mut ref_id_offset: Option<(usize, usize)> = None;
    for (handle, positions) in successors.iter() {
        for (offset, (_, has_predecessor)) in positions.iter().enumerate() {
            if *has_predecessor {
                continue;
            }
            let mut curr = Some(Pos::new(*handle, offset));
            let mut is_ref = false;
            let mut path: Vec<usize> = Vec::new();
            let mut len = 0;
            while let Some(pos) = curr {
                if pos == ref_pos {
                    ref_id_offset = Some((result.len(), path.len()));
                    is_ref = true;
                }
                path.push(pos.node);
                len += subgraph.get(&pos.node).unwrap().sequence_len();
                curr = next_pos(pos, &successors);
            }
            if is_ref {
                if !support::encoded_path_is_canonical(&path) {
                    eprintln!("Warning: the reference path is not in canonical orientation");
                }
                result.push(PathInfo::new(path, len));
            } else if support::encoded_path_is_canonical(&path) {
                result.push(PathInfo::new(path, len));
            }
        }
    }

    let ref_id_offset = ref_id_offset.ok_or("Could not find the reference path".to_string())?;
    Ok((result, ref_id_offset))
}

pub fn distance_to(subgraph: &BTreeMap<usize, GBZRecord>, path: &[usize], path_offset: usize, node_offset: usize) -> usize {
    let mut result = node_offset;
    for handle in path.iter().take(path_offset) {
        result += subgraph.get(handle).unwrap().sequence_len();
    }
    result
}

// Returns all distinct paths and uses the weight field for storing their counts.
// Also updates `ref_id`.
pub fn distinct_paths(
    paths: Vec<PathInfo>,
    ref_id: usize
) -> (Vec<PathInfo>, usize) {
    let ref_path = paths[ref_id].path.clone();
    let mut paths = paths;
    paths.sort_unstable();

    let mut new_paths: Vec<PathInfo> = Vec::new();
    let mut ref_id = 0;
    for info in paths.iter() {
        if new_paths.is_empty() || new_paths.last().unwrap().path != info.path {
            if info.path == ref_path {
                ref_id = new_paths.len();
            }
            new_paths.push(PathInfo::weighted(info.path.clone(), info.len));
        } else {
            new_paths.last_mut().unwrap().increment_weight();
        }
    }

    (new_paths, ref_id)
}

//-----------------------------------------------------------------------------
