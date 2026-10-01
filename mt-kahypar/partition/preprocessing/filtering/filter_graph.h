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

// Undirected weighted graph in CSR format.
struct FilterGraph {
  ds::Array<EdgeID> node_begin;
  ds::Array<NodeID> adj;
  ds::Array<EdgeID> adj_edge;  // half-edge -> edge id
  ds::Array<NodeWeight> node_weight;
  ds::Array<EdgeWeight> edge_weight;

  size_t numNodes() const { return node_weight.size(); }
  size_t numEdges() const { return edge_weight.size(); }
  size_t degree(NodeID v) const { return node_begin[v + 1] - node_begin[v]; }
};

struct EdgeListEntry {
  NodeID u;
  NodeID v;
  EdgeWeight weight;
};

// Edges must not contain self-loops or parallel edges.
FilterGraph build_csr_from_edge_list(
    const std::vector<EdgeListEntry>& edges,
    const std::vector<NodeWeight>& node_weights);

std::vector<std::pair<NodeID, NodeID>> compute_edge_endpoints(
    const FilterGraph& graph);

}  // namespace filtering
}  // namespace mt_kahypar
