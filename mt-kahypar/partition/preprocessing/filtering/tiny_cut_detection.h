#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

namespace mt_kahypar {
namespace filtering {

// PUNCH filtering, part 1: contract everything separated by tiny cuts (1-cuts,
// degree-2 chains, 2-cuts) whose small side weighs at most U.

struct TinyCutParams {
  NodeWeight U;
  NodeWeight tau = 5;
};

// Removing all m edges of a 2-edge-cut class splits the graph into exactly m
// components. All but the largest are found by a bounded traversal; the
// weight of the largest one is inferred by subtraction.
struct BoundedComponentsResult {
  std::vector<std::vector<NodeID>> small_component_members;
  std::vector<NodeWeight> small_component_weight;
  NodeWeight leftover_weight;
  // Only filled if leftover_weight <= U.
  std::vector<NodeID> leftover_members;
};

// `seeds` must contain both endpoints of every excluded edge.
BoundedComponentsResult compute_bounded_components_excluding(
    const FilterGraph& graph, const std::vector<EdgeID>& excluded_edges,
    const std::vector<NodeID>& seeds, NodeWeight total_graph_weight, NodeWeight U);

// Bridge tree rooted at its heaviest block: contracts subtrees of weight <= U.
// Subtrees of weight <= tau are merged into their parent if it stays <= U.
ContractionResult contract_component_tree(const FilterGraph& graph, const TinyCutParams& params);

// Contracts maximal chains of degree-2 vertices of total weight <= U.
ContractionResult contract_degree2_chains(const FilterGraph& graph, NodeWeight U);

// For each 2-edge-cut class S, contracts the components of (V, E \ S) of
// weight <= U.
ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U);

// Runs all three passes; the mapping refers to the input graph.
ContractionResult run_tiny_cut_detection(const FilterGraph& graph, const TinyCutParams& params,
                                          bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
