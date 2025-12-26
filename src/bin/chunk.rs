// FIXME: plan
// * Generate a mapping from node ids to chunk ids
//   * Do boundary nodes go to multiple chunks?
// * Process weakly connected components separately (in parallel?)
//   * For each chunk in the component, create a subgraph and write it as GFA to a file
//   * For each path in the component, trace it and write the corresponding fragments to the chunk GFA files
// * Do a single pass over the GAF file and partition the alignments between chunks
//   * Trace each alignment and clip it into pieces corresponding to the chunks
//   * Write the clipped pieces to the appropriate chunk GAF files

fn main() -> Result<(), String> {
    Ok(())
}