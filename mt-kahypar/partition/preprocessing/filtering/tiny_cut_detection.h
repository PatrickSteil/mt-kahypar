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

}  // namespace filtering
}  // namespace mt_kahypar
