#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

FilterGraph build_csr_from_edge_list(const std::vector<EdgeListEntry>& edges,
                                      const std::vector<NodeWeight>& node_weights) {
  const size_t n = node_weights.size();
  const size_t m = edges.size();

  FilterGraph graph;
  graph.node_begin.resize(n + 1, 0);
  graph.adj.resize(2 * m, kInvalidNode);
  graph.adj_edge.resize(2 * m, kInvalidEdge);
  graph.node_weight.resize(n, 0);
  graph.edge_weight.resize(m, 0);

  for (size_t v = 0; v < n; ++v) graph.node_weight[v] = node_weights[v];
  for (size_t e = 0; e < m; ++e) graph.edge_weight[e] = edges[e].weight;

  std::vector<size_t> degree(n, 0);
  for (const EdgeListEntry& e : edges) {
    ++degree[e.u];
    ++degree[e.v];
  }

  EdgeID offset = 0;
  for (size_t v = 0; v < n; ++v) {
    graph.node_begin[v] = offset;
    offset += static_cast<EdgeID>(degree[v]);
  }
  graph.node_begin[n] = offset;

  std::vector<EdgeID> cursor(n);
  for (size_t v = 0; v < n; ++v) cursor[v] = graph.node_begin[v];

  for (EdgeID e = 0; e < static_cast<EdgeID>(m); ++e) {
    const NodeID u = edges[e].u;
    const NodeID v = edges[e].v;
    graph.adj[cursor[u]] = v;
    graph.adj_edge[cursor[u]] = e;
    ++cursor[u];
    graph.adj[cursor[v]] = u;
    graph.adj_edge[cursor[v]] = e;
    ++cursor[v];
  }

  return graph;
}

std::vector<std::pair<NodeID, NodeID>> compute_edge_endpoints(const FilterGraph& graph) {
  std::vector<std::pair<NodeID, NodeID>> endpoints(
      graph.numEdges(), {kInvalidNode, kInvalidNode});
  for (size_t v = 0; v < graph.numNodes(); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const EdgeID e = graph.adj_edge[pos];
      if (endpoints[e].first == kInvalidNode) {
        endpoints[e].first = static_cast<NodeID>(v);
      } else {
        endpoints[e].second = static_cast<NodeID>(v);
      }
    }
  }
  return endpoints;
}

}  // namespace filtering
}  // namespace mt_kahypar
