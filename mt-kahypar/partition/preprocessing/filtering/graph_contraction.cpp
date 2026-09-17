#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

#include <map>

namespace mt_kahypar {
namespace filtering {

ContractionResult contract_graph(const FilterGraph& graph, const AtomicUnionFind& uf) {
  const size_t n = graph.numNodes();

  std::vector<NodeID> mapping(n, kInvalidNode);
  NodeID next_id = 0;
  for (size_t v = 0; v < n; ++v) {
    const uint32_t rep = uf.find(static_cast<uint32_t>(v));
    if (mapping[rep] == kInvalidNode) mapping[rep] = next_id++;
  }
  for (size_t v = 0; v < n; ++v) {
    mapping[v] = mapping[uf.find(static_cast<uint32_t>(v))];
  }

  std::vector<NodeWeight> node_weights(next_id, 0);
  for (size_t v = 0; v < n; ++v) node_weights[mapping[v]] += graph.node_weight[v];

  // Merge parallel edges (sum weights), drop self-loops created by contraction.
  // std::map keeps this simple/correct; can become a parallel hash
  // aggregation later if profiling shows it matters (tiny-cut detection is a
  // cheap stage of the pipeline overall -- see design spec section 4).
  std::map<std::pair<NodeID, NodeID>, EdgeWeight> merged_edges;
  for (NodeID v = 0; v < static_cast<NodeID>(n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;  // visit every undirected edge exactly once
      const NodeID mu = mapping[v];
      const NodeID mv = mapping[other];
      if (mu == mv) continue;  // self-loop created by contraction
      const std::pair<NodeID, NodeID> key = mu < mv ? std::make_pair(mu, mv) : std::make_pair(mv, mu);
      merged_edges[key] += graph.edge_weight[graph.adj_edge[pos]];
    }
  }

  std::vector<EdgeListEntry> edges;
  edges.reserve(merged_edges.size());
  for (auto& [endpoints, weight] : merged_edges) {
    edges.push_back(EdgeListEntry{endpoints.first, endpoints.second, weight});
  }

  ContractionResult result;
  result.graph = build_csr_from_edge_list(edges, node_weights);
  result.mapping = std::move(mapping);
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
