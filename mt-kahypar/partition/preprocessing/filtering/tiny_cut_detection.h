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
// (V, E \ S) whose total weight is <= U.
ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U);

}  // namespace filtering
}  // namespace mt_kahypar
