use super::*;

use crate::GBZBase;
use crate::formats;

use simple_sds::serialize;

use rand::Rng;

use std::fs;
use std::vec;

//-----------------------------------------------------------------------------

// Alignment.

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

    // Identical bases, different paths.
    {
        let path = vec![1, 5, 3]; // AABAAA
        let ref_path = vec![3, 2, 1, 3]; // AABAAA
        let truth = vec![
            (EditOperation::Match, 6),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "identical bases, different paths");
    }

    // Prefix + suffix length exceeds the length of the shorter path.
    {
        let short = vec![5, 5]; // ABAABA
        let long = vec![1, 2, 3, 5, 3, 2, 1]; // ABAAABAAABA
        let short_path = vec![
            (EditOperation::Match, 4),
            (EditOperation::Deletion, 5),
            (EditOperation::Match, 2),
        ];
        let short_ref = vec![
            (EditOperation::Match, 4),
            (EditOperation::Insertion, 5),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &short, &long, &short_path, "Prefix + suffix length exceeds path length");
        check_edits(&subgraph, &long, &short, &short_ref, "Prefix + suffix length exceeds reference length");
    }
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

#[test]
fn align_with_prefix_with_suffix() {
    let subgraph = create_subgraph();

    // Insertion.
    {
        let path = vec![5, 2, 5]; // ABABABA
        let ref_path = vec![1, 4, 1]; // ABBA
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Insertion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion");
    }

    // Deletion.
    {
        let path = vec![2, 3, 2]; // BAAB
        let ref_path = vec![6, 1, 6]; // BABABAB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Deletion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "deletion");
    }

    // Mismatch + insertion.
    {
        let path = vec![1, 2, 4, 6, 2, 1]; // ABBBBABBA
        let ref_path = vec![5, 5]; // ABAABA
        let truth = vec![
            (EditOperation::Match, 4),
            (EditOperation::Insertion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + insertion");
    }

    // Mismatch + deletion.
    {
        let path = vec![2, 3, 3, 2]; // BAAAAB
        let ref_path = vec![2, 5, 3, 6]; // BABAAABAB
        let truth = vec![
            (EditOperation::Match, 4),
            (EditOperation::Deletion, 3),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "mismatch + deletion");
    }

    // Insertion + deletion.
    {
        let path = vec![1, 5, 4, 3, 4]; // AABABBAABB
        let ref_path = vec![3, 1, 5, 1, 4, 2]; // AAAABAABBB
        let truth = vec![
            (EditOperation::Match, 2),
            (EditOperation::Insertion, 6),
            (EditOperation::Deletion, 6),
            (EditOperation::Match, 2),
        ];
        check_edits(&subgraph, &path, &ref_path, &truth, "insertion + deletion");
    }
}

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

fn check_subgraph(graph: &GBZ, subgraph: &Subgraph, true_nodes: &[usize], path_count: usize, test_case: &str) {
    // Counts.
    assert_eq!(subgraph.nodes(), true_nodes.len(), "Wrong number of nodes for {}", test_case);
    assert_eq!(subgraph.paths(), path_count, "Wrong number of paths for {}", test_case);

    // Minimum and maximum node ids, assuming that `true_nodes` is sorted.
    if !true_nodes.is_empty() {
        assert_eq!(subgraph.min_node(), true_nodes.first().copied(), "Wrong minimum node id for {}", test_case);
        assert_eq!(subgraph.max_node(), true_nodes.last().copied(), "Wrong maximum node id for {}", test_case);
    }

    // Node ids and sequences.
    assert!(subgraph.node_iter().eq(true_nodes.iter().copied()), "Wrong node ids for {}", test_case);
    for node_id in graph.node_iter() {
        if true_nodes.contains(&node_id) {
            assert!(subgraph.has_node(node_id), "Subgraph {} does not contain node {}", test_case, node_id);
            let forward = graph.sequence(node_id).unwrap();
            let reverse = support::reverse_complement(forward);
            assert_eq!(subgraph.sequence(node_id), Some(forward), "Subgraph {} has wrong sequence for node {}", test_case, node_id);
            assert_eq!(subgraph.oriented_sequence(node_id, Orientation::Forward), Some(forward), "Subgraph {} has wrong forward sequence for node {}", test_case, node_id);
            assert_eq!(subgraph.oriented_sequence(node_id, Orientation::Reverse), Some(reverse.as_slice()), "Subgraph {} has wrong reverse sequence for node {}", test_case, node_id);
            assert_eq!(subgraph.sequence_len(node_id), Some(forward.len()), "Subgraph {} has wrong sequence length for node {}", test_case, node_id);
        } else {
            assert!(!subgraph.has_node(node_id), "Subgraph {} contains node {}", test_case, node_id);
            assert!(subgraph.sequence(node_id).is_none(), "Subgraph {} contains sequence for node {}", test_case, node_id);
            assert!(subgraph.oriented_sequence(node_id, Orientation::Forward).is_none(), "Subgraph {} contains forward sequence for node {}", test_case, node_id);
            assert!(subgraph.oriented_sequence(node_id, Orientation::Reverse).is_none(), "Subgraph {} contains reverse sequence for node {}", test_case, node_id);
            assert!(subgraph.sequence_len(node_id).is_none(), "Subgraph {} contains sequence length for node {}", test_case, node_id);
        }
    }

    // Edges.
    for node_id in graph.node_iter() {
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            if subgraph.has_node(node_id) {
                let graph_succ = graph.successors(node_id, orientation).unwrap();
                let subgraph_succ = subgraph.successors(node_id, orientation);
                assert!(subgraph_succ.is_some(), "Subgraph {} does not contain successors for node {} ({})", test_case, node_id, orientation);
                let subgraph_succ = subgraph_succ.unwrap();
                assert!(graph_succ.filter(|(id, _)| subgraph.has_node(*id)).eq(subgraph_succ), "Subgraph {} has wrong successors for node {} ({})", test_case, node_id, orientation);

                let graph_pred = graph.predecessors(node_id, orientation).unwrap();
                let subgraph_pred = subgraph.predecessors(node_id, orientation);
                assert!(subgraph_pred.is_some(), "Subgraph {} does not contain predecessors for node {} ({})", test_case, node_id, orientation);
                let subgraph_pred = subgraph_pred.unwrap();
                assert!(graph_pred.filter(|(id, _)| subgraph.has_node(*id)).eq(subgraph_pred), "Subgraph {} has wrong predecessors for node {} ({})", test_case, node_id, orientation);
            } else {
                assert!(subgraph.successors(node_id, orientation).is_none(), "Subgraph {} contains successors for node {} ({})", test_case, node_id, orientation);
                assert!(subgraph.predecessors(node_id, orientation).is_none(), "Subgraph {} contains predecessors for node {} ({})", test_case, node_id, orientation);
            }
        }
    }
}

//-----------------------------------------------------------------------------

// TODO: Add a multi-node query.
// For each query, returns (true nodes, path count).
fn queries_and_truth() -> (Vec<SubgraphQuery>, Vec<(Vec<usize>, usize)>) {
    let path_a = FullPathName::generic("A");
    let path_b = FullPathName::generic("B");
    let queries = vec![
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::All),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::nodes([14], 1, HaplotypeOutput::Distinct),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::ReferenceOnly),
        SubgraphQuery::path_offset(&path_b, 2, 1, HaplotypeOutput::Distinct),
    ];
    let truth = vec![
        (vec![12, 13, 14, 15, 16], 3),
        (vec![12, 13, 14, 15, 16], 2),
        (vec![12, 13, 14, 15, 16], 2),
        (vec![12, 13, 14, 15, 16], 1),
        (vec![22, 23, 24, 25], 2),
    ];
    (queries, truth)
}

fn queries_and_gfas(cigar: bool) -> (Vec<SubgraphQuery>, Vec<Vec<String>>){
    let path_a = FullPathName::generic("A");
    let path_b = FullPathName::generic("B");
    let queries = vec![
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::All),
        SubgraphQuery::path_offset(&path_a, 2, 1, HaplotypeOutput::Distinct),
        SubgraphQuery::nodes([14], 1, HaplotypeOutput::Distinct),
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
        SubgraphQuery::nodes([14], 1, HaplotypeOutput::Distinct),
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

// Construction and operations.

#[test]
fn random_nodes() {
    let gbz_file = support::get_test_data("example.gbz");
    let graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let min_node = support::node_id(graph.min_node());
    let max_node = support::node_id(graph.max_node());

    let mut selected: BTreeSet<usize> = BTreeSet::new();
    let mut rng = rand::thread_rng();
    let mut subgraph = Subgraph::new();
    for _ in 0..100 {
        let node_id = rng.gen_range(min_node..=max_node);
        if !graph.has_node(node_id) {
            continue;
        }
        if subgraph.has_node(node_id) {
            selected.remove(&node_id);
            subgraph.remove_node(node_id);
        } else {
            selected.insert(node_id);
            let result = subgraph.add_node(node_id, &mut |handle| {
                GBZRecord::from_gbz(&graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
            });
            if let Err(err) = result {
                panic!("Failed to add node {}: {}", node_id, err);
            }
        }
    }

    // We have a subgraph with known nodes and no paths.
    let true_nodes: Vec<usize> = selected.iter().copied().collect();
    check_subgraph(&graph, &subgraph, &true_nodes, 0, "(random nodes)");
}

#[test]
fn subgraph_from_gbz() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    let (queries, truth) = queries_and_truth();
    for (query, (true_nodes, path_count)) in queries.iter().zip(truth.iter()) {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_gbz(&graph, Some(&path_index), query);
        if let Err(err) = result {
            panic!("Failed to create subgraph for query {}: {}", query, err);
        }
        check_subgraph(&graph, &subgraph, &true_nodes, *path_count, &query.to_string());
    }
}

#[test]
fn subgraph_from_db() {
    let gbz_file = support::get_test_data("example.gbz");
    let gbz_graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let db_file = serialize::temp_file_name("subgraph-from-db");
    let result = GBZBase::create_from_file(&gbz_file, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    let mut database = GBZBase::open(&db_file).unwrap();
    let mut graph = GraphInterface::new(&mut database).unwrap();

    let (queries, truth) = queries_and_truth();
    for (query, (true_nodes, path_count)) in queries.iter().zip(truth.iter()) {
        let mut subgraph = Subgraph::new();
        let result = subgraph.from_db(&mut graph, query);
        if let Err(err) = result {
            panic!("Failed to create subgraph for query {}: {}", query, err);
        }
        check_subgraph(&gbz_graph, &subgraph, &true_nodes, *path_count, &query.to_string());
    }

    drop(graph);
    drop(database);
    fs::remove_file(&db_file).unwrap();
}

#[test]
fn manual_gbz_queries() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    let (queries, truth) = queries_and_truth();
    for (query, (true_nodes, path_count)) in queries.iter().zip(truth.iter()) {
        let mut subgraph = Subgraph::new();
        let mut reference_path = None;
        match query.query_type() {
            QueryType::PathOffset(query_pos) => {
                let result = subgraph.path_pos_from_gbz(&graph, &path_index, query_pos);
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
                reference_path = Some(result.unwrap());
                let graph_pos = reference_path.as_ref().unwrap().0.graph_pos();
                let result = subgraph.around_position(graph_pos, query.context, &mut |handle| {
                    GBZRecord::from_gbz(&graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
            QueryType::PathInterval((query_pos, len)) => {
                let result = subgraph.path_pos_from_gbz(&graph, &path_index, query_pos);
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
                reference_path = Some(result.unwrap());
                let start_pos = reference_path.as_ref().unwrap().0;
                let result = subgraph.around_interval(start_pos, *len, query.context, &mut |handle| {
                    GBZRecord::from_gbz(&graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
            QueryType::Nodes(nodes) => {
                let result = subgraph.around_nodes(nodes, query.context, &mut |handle| {
                    GBZRecord::from_gbz(&graph, handle).ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
        }

        // We do not have paths yet.
        assert_eq!(subgraph.paths(), 0, "Subgraph {} has paths", query);
        let result = subgraph.extract_paths(reference_path, query.output());
        if let Err(err) = result {
            panic!("Path extraction for query {} failed: {}", query, err);
        }
        check_subgraph(&graph, &subgraph, &true_nodes, *path_count, &query.to_string());
    }
}

#[test]
fn manual_db_queries() {
    let gbz_file = support::get_test_data("example.gbz");
    let gbz_graph: GBZ = serialize::load_from(&gbz_file).unwrap();
    let db_file = serialize::temp_file_name("subgraph-from-db");
    let result = GBZBase::create_from_file(&gbz_file, &db_file);
    assert!(result.is_ok(), "Failed to create database: {}", result.unwrap_err());
    let mut database = GBZBase::open(&db_file).unwrap();
    let mut graph = GraphInterface::new(&mut database).unwrap();

    let (queries, truth) = queries_and_truth();
    for (query, (true_nodes, path_count)) in queries.iter().zip(truth.iter()) {
        let mut subgraph = Subgraph::new();
        let mut reference_path = None;
        match query.query_type() {
            QueryType::PathOffset(query_pos) => {
                let result = subgraph.path_pos_from_db(&mut graph, query_pos);
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
                reference_path = Some(result.unwrap());
                let graph_pos = reference_path.as_ref().unwrap().0.graph_pos();
                let result = subgraph.around_position(graph_pos, query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
            QueryType::PathInterval((query_pos, len)) => {
                let result = subgraph.path_pos_from_db(&mut graph, query_pos);
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
                reference_path = Some(result.unwrap());
                let start_pos = reference_path.as_ref().unwrap().0;
                let result = subgraph.around_interval(start_pos, *len, query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
            QueryType::Nodes(nodes) => {
                let result = subgraph.around_nodes(nodes, query.context, &mut |handle| {
                    let record = graph.get_record(handle)?;
                    record.ok_or(format!("The graph does not contain handle {}", handle))
                });
                if let Err(err) = result {
                    panic!("Query {} failed: {}", query, err);
                }
            }
        }

        // We do not have paths yet.
        assert_eq!(subgraph.paths(), 0, "Subgraph {} has paths", query);
        let result = subgraph.extract_paths(reference_path, query.output());
        if let Err(err) = result {
            panic!("Path extraction for query {} failed: {}", query, err);
        }
        check_subgraph(&gbz_graph, &subgraph, &true_nodes, *path_count, &query.to_string());
    }

    drop(graph);
    drop(database);
    fs::remove_file(&db_file).unwrap();
}

//-----------------------------------------------------------------------------

// GFA/JSON output.

#[test]
fn gfa_output() {
    let (graph, path_index) = gbz_and_path_index("example.gbz", GBZBase::INDEX_INTERVAL);
    for cigar in [false, true] {
        let (queries, gfas) = queries_and_gfas(cigar);
        for (query, truth) in queries.iter().zip(gfas.iter()) {
            let mut subgraph = Subgraph::new();
            let _ = subgraph.from_gbz(&graph, Some(&path_index), query);
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
            let mut subgraph = Subgraph::new();
            let _ = subgraph.from_gbz(&graph, Some(&path_index), query);
            let mut output = Vec::new();
            let result = subgraph.write_json(&mut output, cigar);
            assert!(result.is_ok(), "Failed to write JSON for query {}: {}", query, result.unwrap_err());
            let json = String::from_utf8(output).unwrap();
            assert_eq!(json, *truth, "Wrong JSON output for query {}", query);
        }
    }
}

//-----------------------------------------------------------------------------
