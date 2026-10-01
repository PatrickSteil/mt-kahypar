#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

struct EdgeSignatures {
  std::vector<unsigned __int128> signature;
  std::vector<char> is_tree_edge;
};

// Pritchard & Thurimella: every non-tree edge gets a random 128-bit label, a
// tree edge the XOR of the labels of all non-tree edges covering it.
EdgeSignatures compute_edge_signatures(const FilterGraph& graph, const SpanningForest& forest);

// Bridges are the tree edges with signature 0.
std::vector<char> compute_bridges(const FilterGraph& graph, const EdgeSignatures& sigs);

// 2-edge-cut equivalence classes: non-bridge edges with equal signatures.
// Any two edges of a class form a 2-edge cut.
std::vector<std::vector<EdgeID>> find_two_edge_cut_classes(const EdgeSignatures& sigs);

}  // namespace filtering
}  // namespace mt_kahypar
