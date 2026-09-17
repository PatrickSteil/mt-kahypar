#pragma once

#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "mt-kahypar/datastructures/array.h"

namespace mt_kahypar {
namespace filtering {

using NodeID = uint32_t;
using EdgeID = uint32_t;
using NodeWeight = uint64_t;
using EdgeWeight = int64_t;

constexpr NodeID kInvalidNode = std::numeric_limits<NodeID>::max();
constexpr EdgeID kInvalidEdge = std::numeric_limits<EdgeID>::max();

// Plain CSR representation of an undirected weighted graph, independent of
// Mt-KaHyPar's hypergraph/partitioning types (see design spec section 3).
struct FilterGraph {
  ds::Array<EdgeID> node_begin;       // size numNodes() + 1
  ds::Array<NodeID> adj;              // size 2 * numEdges(), neighbor per half-edge
  ds::Array<EdgeID> adj_edge;         // size 2 * numEdges(), half-edge -> canonical edge id
  ds::Array<NodeWeight> node_weight;  // size numNodes()
  ds::Array<EdgeWeight> edge_weight;  // size numEdges()

  size_t numNodes() const { return node_weight.size(); }
  size_t numEdges() const { return edge_weight.size(); }
  size_t degree(NodeID v) const { return node_begin[v + 1] - node_begin[v]; }
};

// One undirected edge, used while building a FilterGraph from an edge list.
struct EdgeListEntry {
  NodeID u;
  NodeID v;
  EdgeWeight weight;
};

// Builds a FilterGraph in CSR form. node_weights.size() determines the number
// of nodes; every edge must reference node ids in [0, node_weights.size()).
// Callers must not pass self-loops (u == v) or duplicate edges between the
// same pair; graph_contraction.cpp is responsible for merging/dropping those
// when a contraction would create them.
FilterGraph build_csr_from_edge_list(const std::vector<EdgeListEntry>& edges,
                                      const std::vector<NodeWeight>& node_weights);

// Recovers the (u, v) endpoint pair for every canonical edge id, by scanning
// the adjacency structure once (O(numNodes() + numEdges())).
std::vector<std::pair<NodeID, NodeID>> compute_edge_endpoints(const FilterGraph& graph);

}  // namespace filtering
}  // namespace mt_kahypar
