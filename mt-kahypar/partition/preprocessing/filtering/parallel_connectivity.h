// mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct SpanningForest {
  std::vector<NodeID> parent;       // BFS-tree parent; roots are their own parent
  std::vector<EdgeID> parent_edge;  // canonical edge id to parent; kInvalidEdge for roots
  std::vector<NodeID> bfs_order;    // visitation order, non-decreasing depth within each root
  std::vector<NodeID> roots;
};

// Computes connected components of `graph` in parallel via TBB parallel_for
// unioning every adjacency pair through an atomic union-find. Component ids
// are union-find representatives (small vertex ids), not a dense range. See
// design spec section 4.2 for why this simple approach (rather than
// Shiloach-Vishkin pointer jumping) is adequate here.
std::vector<NodeID> parallel_connected_components(const FilterGraph& graph);

// Same as above, but treats every edge e with excluded_edge[e] == true as
// absent. excluded_edge must have size graph.numEdges().
std::vector<NodeID> parallel_connected_components_excluding(
    const FilterGraph& graph, const std::vector<char>& excluded_edge);

// Builds a BFS spanning forest given a graph and its component labels (as
// produced by either function above). If `roots` is empty, one root per
// distinct component value is chosen automatically (its lowest-numbered
// vertex); otherwise `roots` must contain exactly one vertex per distinct
// component value. BFS runs sequentially per component, but components are
// processed in parallel via TBB (safe: components touch disjoint index
// ranges of the output arrays).
SpanningForest build_spanning_forest(const FilterGraph& graph,
                                      const std::vector<NodeID>& component,
                                      const std::vector<NodeID>& roots);

}  // namespace filtering
}  // namespace mt_kahypar
