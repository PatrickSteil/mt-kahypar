// mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

namespace mt_kahypar {
namespace filtering {

struct TinyCutParams {
  NodeWeight U;
  NodeWeight tau = 5;
};

// Result of a bounded-traversal connected-components computation that
// excludes a small edge set `excluded_edges` (a 2-edge-cut equivalence
// class of size m, m >= 2). Every genuine, non-bridge-containing class this
// module produces splits the graph into EXACTLY m components when all its
// edges are removed (proof sketch: viewed through the lens of just this
// class's edges, contracting everything else, the class's edges must form
// a simple cycle rather than a simple path -- a path's two end edges would
// individually be bridges, which the class-finding step already excludes;
// a cycle of m edges over m "virtual" super-nodes yields exactly m pieces
// when fully cut). This lets us find all but the single LARGEST resulting
// component explicitly, each touching only ~2x its own size of work (the
// paper's "two-at-a-time" bound, design spec section 4.4 step 8), and infer
// the largest one's weight by subtraction from the graph's total weight
// instead of an O(n) scan -- never touching the majority of the graph.
struct BoundedComponentsResult {
  // One entry per explicitly-found (i.e. not the inferred-by-elimination
  // largest) component.
  std::vector<std::vector<NodeID>> small_component_members;
  std::vector<NodeWeight> small_component_weight;
  // total_graph_weight - sum(small_component_weight). The single remaining
  // component's weight, known without ever visiting it.
  NodeWeight leftover_weight;
  // Only populated (via a fallback O(n) scan) in the rare case that even
  // the leftover component's weight is <= U and its exact membership is
  // therefore needed; empty otherwise, since it is normally far too heavy
  // to be a contraction candidate and materializing it would defeat the
  // whole point of this function.
  std::vector<NodeID> leftover_members;
};

// `excluded_edges` is the class's edge set; `seeds` must contain both
// endpoints of every edge in `excluded_edges` (duplicates are fine and are
// deduplicated internally). `total_graph_weight` is the sum of all of
// `graph`'s vertex weights, computed once by the caller across all classes,
// not per call.
BoundedComponentsResult compute_bounded_components_excluding(
    const FilterGraph& graph, const std::vector<EdgeID>& excluded_edges,
    const std::vector<NodeID>& seeds, NodeWeight total_graph_weight, NodeWeight U);

// Part 1, pass 1 (design spec section 4.2, corrected per the "Correction
// from the spec" note in this plan's Global Constraints): builds the bridge
// tree (nodes = 2-edge-connected blocks, edges = bridges), roots it at its
// heaviest block, and contracts any subtree with total weight <= U top-down.
// A contracted subtree of weight <= tau is additionally folded into its
// parent block if the merged weight is still <= U.
ContractionResult contract_component_tree(const FilterGraph& graph, const TinyCutParams& params);

// Part 1, pass 2 (design spec section 4.3): contracts each maximal chain of
// degree-2 vertices into a single vertex, provided the chain's total weight
// is <= U (whole-chain-or-nothing; a chain that's too big is left as-is). A
// connected component made entirely of degree-2 vertices (a pure cycle, with
// no anchor of different degree) is treated as one chain.
ContractionResult contract_degree2_chains(const FilterGraph& graph, NodeWeight U);

// Part 1, pass 3 (design spec section 4.4): finds 2-edge-cut equivalence
// classes and, for each class S, contracts every connected component of
// (V, E \ S) whose total weight is <= U. Classes are processed in parallel
// (see the implementation comment), and each class's components are found
// via the bounded "two-at-a-time" traversal (compute_bounded_components_
// excluding above) rather than a full O(n) graph scan -- real road-network
// graphs can produce tens of thousands of classes (Task 19's DIMACS
// validation found ~38K on a 264K-vertex instance), and a full scan per
// class, even parallelized across classes, still does O(classes * n) total
// work and does not finish in reasonable time.
ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U);

// Runs all three tiny-cut passes in sequence (design spec section 4),
// composing their node mappings so the result maps the ORIGINAL input
// graph's vertices directly to the final, smallest output graph's vertices.
// When `verbose` is set, prints each pass's name, timing, and resulting
// vertex count to stderr as it runs.
ContractionResult run_tiny_cut_detection(const FilterGraph& graph, const TinyCutParams& params,
                                          bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
