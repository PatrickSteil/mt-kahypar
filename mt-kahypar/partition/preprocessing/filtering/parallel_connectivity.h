#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct SpanningForest {
  std::vector<NodeID> parent;       // roots are their own parent
  std::vector<EdgeID> parent_edge;  // kInvalidEdge for roots
  std::vector<NodeID> bfs_order;
  std::vector<NodeID> roots;
};

// Component ids are union-find representatives, not a dense range.
std::vector<NodeID> parallel_connected_components(const FilterGraph& graph);

std::vector<NodeID> parallel_connected_components_excluding(
    const FilterGraph& graph, const std::vector<char>& excluded_edge);

// BFS forest, one tree per component (processed in parallel). `component`
// must come from parallel_connected_components() on the same graph. If
// `roots` is empty, the smallest vertex of each component is used.
SpanningForest build_spanning_forest(const FilterGraph& graph,
                                      const std::vector<NodeID>& component,
                                      const std::vector<NodeID>& roots);

}  // namespace filtering
}  // namespace mt_kahypar
