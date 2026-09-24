#include "mt-kahypar/partition/preprocessing/filtering/metis_io.h"

#include <tbb/parallel_sort.h>

#include <fstream>
#include <stdexcept>

#include "mt-kahypar/io/hypergraph_io.h"

namespace mt_kahypar {
namespace filtering {

FilterGraph read_metis_graph(const std::string& path) {
  HyperedgeID num_edges = 0;
  HypernodeID num_nodes = 0;
  io::EdgeVector edges;
  vec<HyperedgeWeight> edge_weights;
  vec<HypernodeWeight> node_weights;
  io::readGraphFile(path, num_edges, num_nodes, edges, edge_weights, node_weights);

  // readGraphFile yields every undirected edge once, as (u, v) with u < v
  std::vector<EdgeListEntry> edge_list(edges.size());
  tbb::parallel_for(size_t(0), edges.size(), [&](const size_t e) {
    const EdgeWeight weight = edge_weights.empty() ? 1 : static_cast<EdgeWeight>(edge_weights[e]);
    edge_list[e] = EdgeListEntry{static_cast<NodeID>(edges[e].first),
                                 static_cast<NodeID>(edges[e].second), weight};
  });

  tbb::parallel_sort(edge_list.begin(), edge_list.end(), [](const EdgeListEntry& a, const EdgeListEntry& b) {
    return a.u < b.u || (a.u == b.u && a.v < b.v);
  });
  size_t num_unique = 0;
  for (size_t e = 0; e < edge_list.size(); ++e) {
    if (num_unique > 0 && edge_list[num_unique - 1].u == edge_list[e].u &&
        edge_list[num_unique - 1].v == edge_list[e].v) {
      edge_list[num_unique - 1].weight += edge_list[e].weight;
    } else {
      edge_list[num_unique++] = edge_list[e];
    }
  }
  edge_list.resize(num_unique);

  std::vector<NodeWeight> weights(num_nodes, 1);
  if (!node_weights.empty()) {
    for (size_t v = 0; v < num_nodes; ++v) weights[v] = static_cast<NodeWeight>(node_weights[v]);
  }
  return build_csr_from_edge_list(edge_list, weights);
}

void write_metis_graph(const FilterGraph& graph, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not write METIS graph file: " + path);
  out << graph.numNodes() << " " << graph.numEdges() << " 11\n";
  for (NodeID v = 0; v < static_cast<NodeID>(graph.numNodes()); ++v) {
    out << graph.node_weight[v];
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      out << " " << (graph.adj[pos] + 1) << " " << graph.edge_weight[graph.adj_edge[pos]];
    }
    out << "\n";
  }
}

}  // namespace filtering
}  // namespace mt_kahypar
