use super::*;

use crate::GBZBase;
use crate::formats;

use simple_sds::serialize;

use std::fs;

//-----------------------------------------------------------------------------

// Synthetic tests for Subgraph internals.

fn subgraph_from_sequences(nodes: &[(usize, Vec<u8>)]) -> Subgraph {
    let mut records: BTreeMap<usize, GBZRecord> = BTreeMap::new();
    for (handle, sequence) in nodes.iter() {
        let record = unsafe {
            GBZRecord::from_raw_parts(*handle, Vec::new(), Vec::new(), sequence.clone())
        };
        records.insert(*handle, record);
    }
    Subgraph {
        records,
        paths: Vec::new(),
        ref_id: None,
        ref_path: None,
        ref_interval: None,
    }
}

fn create_subgraph() -> Subgraph {
    let nodes = vec![
        (1, b"A".to_vec()),
        (2, b"B".to_vec()),
        (3, b"AA".to_vec()),
        (4, b"BB".to_vec()),
        (5, b"ABA".to_vec()),
        (6, b"BAB".to_vec()),
    ];
    subgraph_from_sequences(&nodes)
}

fn check_edits(subgraph: &Subgraph, path: &[usize], ref_path: &[usize], truth: &[(EditOperation, usize)], name: &str) {
    let mut edits = Vec::new();
    subgraph.align(path, ref_path, &mut edits);
    assert_eq!(edits.len(), truth.len(), "Wrong number of edits for {}", name);
    for (i, (edit, truth_edit)) in edits.iter().zip(truth.iter()).enumerate() {
        assert_eq!(edit, truth_edit, "Wrong edit {} for {}", i, name);
    }
}

#[test]
fn align_special_cases() {
    let subgraph = create_subgraph();
    let empty = Vec::new();
    let non_empty = vec![5, 6];

    // (empty, empty)
    {
        let truth = Vec::new();
        check_edits(&subgraph, &empty, &empty, &truth, "empty paths");
    }

    // (empty, non-empty)
    {
        let truth = vec![(EditOperation::Deletion, 6)];
        check_edits(&subgraph, &empty, &non_empty, &truth, "empty vs. non-empty paths");
    }

    // (non-empty, empty)
    {
        let truth = vec![(EditOperation::Insertion, 6)];
        check_edits(&subgraph, &non_empty, &empty, &truth, "non-empty vs. empty paths");
    }

    // (non-empty, non-empty)
    {
        let truth = vec![(EditOperation::Match, 6)];
        check_edits(&subgraph, &non_empty, &non_empty, &truth, "identical paths");
    }

    // FIXME identical bases, different paths

    // FIXME prefix + suffix length exceeds the length of the shorter path
}

#[test]
fn align_no_prefix_no_suffix() {
    let subgraph = create_subgraph();

    // Insertion and deletion are in `align_special_cases`.

    // Mismatch + insertion.
    {
        let path = vec![1, 2, 5]; // ABABA
        let ref_path = vec![4]; // BB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Insertion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + insertion");
    }

    // Mismatch + deletion.
    {
        let path = vec![3]; // AA
        let ref_path = vec![2, 1, 6]; // BABAB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Deletion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + deletion");
    }

    // Insertion + deletion.
    {
        let path = vec![1, 2, 5, 3]; // ABABAAA
        let ref_path = vec![4, 2, 1, 6]; // BBBABAB
        let truth = vec![
            (EditOperation::Insertion, 7),
            (EditOperation::Deletion, 7),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion + deletion");
    }
}

#[test]
fn align_no_prefix_with_suffix() {
    let subgraph = create_subgraph();

    // Insertion.
    {
        let path = vec![1, 2, 5]; // ABABA
        let ref_path = vec![2, 1]; // BA
        let truth = vec![
            (EditOperation::Insertion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion");
    }

    // Deletion.
    {
        let path = vec![1, 2]; // AB
        let ref_path = vec![2, 1, 6]; // BABAB
        let truth = vec![
            (EditOperation::Deletion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "deletion");
    }

    // Mismatch + insertion.
    {
        let path = vec![5, 2, 5]; // ABABABA
        let ref_path = vec![4, 2, 1]; // BBBA
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Insertion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + insertion");
    }

    // Mismatch + deletion.
    {
        let path = vec![3, 1, 2]; // AAAB
        let ref_path = vec![6, 1, 6]; // BABABAB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Deletion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + deletion");
    }

    // Insertion + deletion.
    {
        let path = vec![3, 2, 5, 5]; // AABABAABA
        let ref_path = vec![4, 2, 1, 6, 1]; // BBBABABA
        let truth = vec![
            (EditOperation::Insertion, 6),
            (EditOperation::Deletion, 5),
            (EditOperation::Match, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion + deletion");
    }
}

#[test]
fn align_with_prefix_no_suffix() {
    let subgraph = create_subgraph();

    // Insertion.
    {
        let path = vec![5, 2, 1]; // ABABA
        let ref_path = vec![1, 2]; // AB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Insertion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion");
    }

    // Deletion.
    {
        let path = vec![2, 1]; // BA
        let ref_path = vec![6, 1, 2]; // BABAB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Deletion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "deletion");
    }

    // Mismatch + insertion.
    {
        let path = vec![5, 2, 5]; // ABABABA
        let ref_path = vec![1, 2, 4]; // ABBB
        let truth = vec![
            (EditOperation::Match, 4),
            (EditOperation::Insertion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + insertion");
    }

    // Mismatch + deletion.
    {
        let path = vec![2, 1, 3]; // BAAA
        let ref_path = vec![6, 1, 6]; // BABABAB
        let truth = vec![
            (EditOperation::Match, 4),
            (EditOperation::Deletion, 3),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + deletion");
    }

    // Insertion + deletion.
    {
        let path = vec![3, 2, 5, 5]; // AABABAABA
        let ref_path = vec![1, 1, 4, 2, 1, 6, 2]; // AABBBABABB
        let truth = vec![
            (EditOperation::Match, 3),
            (EditOperation::Insertion, 6),
            (EditOperation::Deletion, 7),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion + deletion");
    }
}

// FIXME prefix and suffix

//-----------------------------------------------------------------------------

// TODO: We should also have a graph with reference paths for testing.
// TODO: We should also have a graph with fragmented reference paths for testing.
// TODO: We should also have a graph with longer nodes for testing.

fn gbz_and_path_index(filename: &'static str, interval: usize) -> (GBZ, PathIndex) {
    let gbz_file = support::get_test_data(filename);
    let graph: GBZ = serialize::load_from(gbz_file).unwrap();
    let path_index = PathIndex::new(&graph, interval, false);
    if let Err(err) = path_index {
        panic!("Failed to create path index with interval {}: {}", interval, err);
    }
    (graph, path_index.unwrap())
}

fn queries_and_counts() -> (Vec<SubgraphQuery>, Vec<(usize, usize)>) {
    let path_a = FullPathName::generic("A");
    let path_b = FullPathName::generic("B");
    let queries = vec![
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::All),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::node(14, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::ReferenceOnly),
        SubgraphQuery::path_offset(&path_b, 2, 1, HaplotypeOutput::Distinct),
    ];
    let counts = vec![
        (5, 3),
        (5, 2),
        (5, 2),
        (5, 1),
        (4, 2),
    ];
    (queries, counts)
}

fn queries_and_gfas(cigar: bool) -> (Vec<SubgraphQuery>, Vec<Vec<String>>){
    let path_a = FullPathName::generic("A");
    let path_b = FullPathName::generic("B");
    let queries = vec![
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::All),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::node(14, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::ReferenceOnly),
        SubgraphQuery::path_offset(&path_b, 2, 1, HaplotypeOutput::Distinct),
    ];
    let gfas = vec![
        vec![
            String::from("H\tVN:Z:1.1\tRS:Z:_gbwt_ref"),
            String::from("S\t12\tA"),
            String::from("S\t13\tT"),
            String::from("S\t14\tT"),
            String::from("S\t15\tA"),
            String::from("S\t16\tC"),
            String::from("L\t12\t+\t14\t+\t0M"),
            String::from("L\t13\t+\t14\t+\t0M"),
            String::from("L\t14\t+\t15\t+\t0M"),
            String::from("L\t14\t+\t16\t+\t0M"),
            String::from("W\t_gbwt_ref\t0\tA\t1\t4\t>12>14>15"),
            format!("W\tunknown\t1\tA\t0\t3\t>12>14>15{}", if cigar { "\tCG:Z:3M" } else { "" }),
            format!("W\tunknown\t2\tA\t0\t3\t>13>14>16{}", if cigar { "\tCG:Z:3M" } else { "" }),
        ],
        vec![
            String::from("H\tVN:Z:1.1\tRS:Z:_gbwt_ref"),
            String::from("S\t12\tA"),
            String::from("S\t13\tT"),
            String::from("S\t14\tT"),
            String::from("S\t15\tA"),
            String::from("S\t16\tC"),
            String::from("L\t12\t+\t14\t+\t0M"),
            String::from("L\t13\t+\t14\t+\t0M"),
            String::from("L\t14\t+\t15\t+\t0M"),
            String::from("L\t14\t+\t16\t+\t0M"),
            String::from("W\t_gbwt_ref\t0\tA\t1\t4\t>12>14>15\tWT:i:2"),
            format!("W\tunknown\t1\tA\t0\t3\t>13>14>16\tWT:i:1{}", if cigar { "\tCG:Z:3M" } else { "" }),
        ],
        vec![
            String::from("H\tVN:Z:1.1"),
            String::from("S\t12\tA"),
            String::from("S\t13\tT"),
            String::from("S\t14\tT"),
            String::from("S\t15\tA"),
            String::from("S\t16\tC"),
            String::from("L\t12\t+\t14\t+\t0M"),
            String::from("L\t13\t+\t14\t+\t0M"),
            String::from("L\t14\t+\t15\t+\t0M"),
            String::from("L\t14\t+\t16\t+\t0M"),
            String::from("W\tunknown\t1\tunknown\t0\t3\t>12>14>15\tWT:i:2"),
            String::from("W\tunknown\t2\tunknown\t0\t3\t>13>14>16\tWT:i:1"),
        ],
        vec![
            String::from("H\tVN:Z:1.1\tRS:Z:_gbwt_ref"),
            String::from("S\t12\tA"),
            String::from("S\t13\tT"),
            String::from("S\t14\tT"),
            String::from("S\t15\tA"),
            String::from("S\t16\tC"),
            String::from("L\t12\t+\t14\t+\t0M"),
            String::from("L\t13\t+\t14\t+\t0M"),
            String::from("L\t14\t+\t15\t+\t0M"),
            String::from("L\t14\t+\t16\t+\t0M"),
            String::from("W\t_gbwt_ref\t0\tA\t1\t4\t>12>14>15"),
        ],
        vec![
            String::from("H\tVN:Z:1.1\tRS:Z:_gbwt_ref"),
            String::from("S\t22\tA"),
            String::from("S\t23\tT"),
            String::from("S\t24\tT"),
            String::from("S\t25\tA"),
            String::from("L\t22\t+\t24\t+\t0M"),
            String::from("L\t23\t+\t24\t-\t0M"),
            String::from("L\t24\t+\t25\t+\t0M"),
            String::from("W\t_gbwt_ref\t0\tB\t1\t4\t>22>24>25\tWT:i:2"),
            format!("W\tunknown\t1\tB\t0\t3\t>22>24<23\tWT:i:1{}", if cigar { "\tCG:Z:3M" } else { "" }),
        ]
    ];
    (queries, gfas)
}

// Fewer queries here, because JSON literals are so inconvenient to generate.
fn queries_and_jsons(cigar: bool) -> (Vec<SubgraphQuery>, Vec<String>){
    let path_a = FullPathName::generic("A");
    //let path_b = FullPathName::generic("B");
    let queries = vec![
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::All),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::node(14, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::ReferenceOnly),
    ];

    let mut nodes: Vec<JSONValue> = Vec::new();
    for (id, sequence) in [("12", "A"), ("13", "T"), ("14", "T"), ("15", "A"), ("16", "C")] {
        nodes.push(JSONValue::Object(vec![
            (String::from("id"), JSONValue::String(String::from(id))),
            (String::from("sequence"), JSONValue::String(String::from(sequence))),
        ]));
    }

    let mut edges: Vec<JSONValue> = Vec::new();
    for (from, from_is_reverse, to, to_is_reverse) in [
        ("12", false, "14", false),
        ("13", false, "14", false),
        ("14", false, "15", false),
        ("14", false, "16", false),
    ] {
        edges.push(JSONValue::Object(vec![
            (String::from("from"), JSONValue::String(String::from(from))),
            (String::from("from_is_reverse"), JSONValue::Boolean(from_is_reverse)),
            (String::from("to"), JSONValue::String(String::from(to))),
            (String::from("to_is_reverse"), JSONValue::Boolean(to_is_reverse)),
        ]));
    }

    let mut jsons: Vec<String> = Vec::new();
    let all: Vec<(Vec<usize>, WalkMetadata, Option<usize>, Option<String>)> = vec![
        (vec![24, 28, 30], WalkMetadata::path_interval(&path_a, 1..4.clone()), None, None),
        (vec![24, 28, 30], WalkMetadata::anonymous(1, "A", 3), None, Some(String::from("3M"))),
        (vec![26, 28, 32], WalkMetadata::anonymous(2, "A", 3), None, Some(String::from("3M"))),
    ];
    let distinct_path: Vec<(Vec<usize>, WalkMetadata, Option<usize>, Option<String>)> = vec![
        (vec![24, 28, 30], WalkMetadata::path_interval(&path_a, 1..4.clone()), Some(2), None),
        (vec![26, 28, 32], WalkMetadata::anonymous(1, "A", 3), Some(1), Some(String::from("3M"))),
    ];
    let distinct_node: Vec<(Vec<usize>, WalkMetadata, Option<usize>, Option<String>)> = vec![
        (vec![24, 28, 30], WalkMetadata::anonymous(1, "unknown", 3), Some(2), None),
        (vec![26, 28, 32], WalkMetadata::anonymous(2, "unknown", 3), Some(1), None),
    ];
    let ref_only: Vec<(Vec<usize>, WalkMetadata, Option<usize>, Option<String>)> = vec![
        (vec![24, 28, 30], WalkMetadata::path_interval(&path_a, 1..4.clone()), None, None),
    ];
    for paths in [all, distinct_path, distinct_node, ref_only] {
        let mut json_paths: Vec<JSONValue> = Vec::new();
        for (path, mut metadata, weight, cigar_string) in paths {
            metadata.add_weight(weight);
            if cigar {
                metadata.add_cigar(cigar_string);
            }
            let json_path = formats::json_path(&path, &metadata);
            json_paths.push(json_path);
        }
        jsons.push(JSONValue::Object(
            vec![
                (String::from("nodes"), JSONValue::Array(nodes.clone())),
                (String::from("edges"), JSONValue::Array(edges.clone())),
                (String::from("paths"), JSONValue::Array(json_paths)),
            ]
        ).to_string());
    }

    (queries, jsons)
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

#[test]
fn subgraph_from_gbz() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    let (queries, nodes_and_paths) = queries_and_counts();
    for (query, (node_count, path_count)) in queries.iter().zip(nodes_and_paths.iter()) {
        let subgraph = Subgraph::from_gbz(&graph, Some(&path_index), query);
        if let Err(err) = subgraph {
            panic!("Failed to create subgraph for query {}: {}", query, err);
        }
        let subgraph = subgraph.unwrap();
        assert_eq!(subgraph.node_count(), *node_count, "Wrong node count for query {}", query);
        assert_eq!(subgraph.path_count(), *path_count, "Wrong path count for query {}", query);
    }
}

#[test]
fn subgraph_from_db() {
    let gbz_file = support::get_test_data("example.gbz");
    let db_file = serialize::temp_file_name("subgraph-from-db");
    let result = GBZBase::create_from_file(&gbz_file, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    let mut database = GBZBase::open(&db_file).unwrap();
    let mut graph = GraphInterface::new(&mut database).unwrap();

    let (queries, nodes_and_paths) = queries_and_counts();
    for (query, (node_count, path_count)) in queries.iter().zip(nodes_and_paths.iter()) {
        let subgraph = Subgraph::from_db(&mut graph, query);
        if let Err(err) = subgraph {
            panic!("Failed to create subgraph for query {}: {}", query, err);
        }
        let subgraph = subgraph.unwrap();
        assert_eq!(subgraph.node_count(), *node_count, "Wrong node count for query {}", query);
        assert_eq!(subgraph.path_count(), *path_count, "Wrong path count for query {}", query);
    }

    drop(graph);
    drop(database);
    fs::remove_file(&db_file).unwrap();
}

//-----------------------------------------------------------------------------

#[test]
fn gfa_output() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    for cigar in [false, true] {
        let (queries, gfas) = queries_and_gfas(cigar);
        for (query, truth) in queries.iter().zip(gfas.iter()) {
            let subgraph = Subgraph::from_gbz(&graph, Some(&path_index), query).unwrap();
            let mut output = Vec::new();
            let result = subgraph.write_gfa(&mut output, cigar);
            assert!(result.is_ok(), "Failed to write GFA for query {}: {}", query, result.unwrap_err());
            let gfa = String::from_utf8(output).unwrap();
            let lines: Vec<&str> = gfa.lines().collect();
            assert_eq!(lines.len(), truth.len(), "Wrong number of lines in GFA output for query {}", query);
            for (line_num, (line, truth_line)) in lines.iter().zip(truth.iter()).enumerate() {
                assert_eq!(line, truth_line, "Wrong line {} in GFA output for query {}", line_num + 1, query);
            }
        }
    }
}

#[test]
fn json_output() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    for cigar in [false, true] {
        let (queries, jsons) = queries_and_jsons(cigar);
        for (query, truth) in queries.iter().zip(jsons.iter()) {
            let subgraph = Subgraph::from_gbz(&graph, Some(&path_index), query).unwrap();
            let mut output = Vec::new();
            let result = subgraph.write_json(&mut output, cigar);
            assert!(result.is_ok(), "Failed to write JSON for query {}: {}", query, result.unwrap_err());
            let json = String::from_utf8(output).unwrap();
            assert_eq!(json, *truth, "Wrong JSON output for query {}", query);
        }
    }
}

//-----------------------------------------------------------------------------
