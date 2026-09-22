#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"

#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

FilteringResult assemble_fragments(const FilterGraph& graph,
                                    const std::vector<char>& keep,
                                    const std::vector<NodeID>& part1_mapping) {
  const size_t graph_n = graph.numNodes();
  const size_t orig_n = part1_mapping.size();

  AtomicUnionFind uf(graph_n);
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;
      if (!keep[graph.adj_edge[pos]]) uf.unite(v, other);
    }
  }

  std::vector<NodeID> fragment_of_graph_vertex(graph_n, kInvalidNode);
  NodeID next_fragment = 0;
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    const uint32_t rep = uf.find(v);
    if (fragment_of_graph_vertex[rep] == kInvalidNode) fragment_of_graph_vertex[rep] = next_fragment++;
    fragment_of_graph_vertex[v] = fragment_of_graph_vertex[rep];
  }

  FilteringResult result;
  result.fragment_id.resize(orig_n);
  for (size_t orig = 0; orig < orig_n; ++orig) {
    result.fragment_id[orig] = fragment_of_graph_vertex[part1_mapping[orig]];
  }

  // Fragment sizes are the sum of Part-1-output vertex weights per fragment,
  // which already equals the sum of ORIGINAL vertex weights, since Part 1's
  // contraction sums weights transitively (see graph_contraction.cpp) --
  // original weights aren't needed directly here.
  result.fragment_size.assign(next_fragment, 0);
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    result.fragment_size[fragment_of_graph_vertex[v]] += graph.node_weight[v];
  }

  std::vector<NodeID> representative(graph_n, kInvalidNode);
  for (size_t orig = 0; orig < orig_n; ++orig) {
    NodeID& rep = representative[part1_mapping[orig]];
    if (rep == kInvalidNode) rep = static_cast<NodeID>(orig);
  }

  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;
      if (keep[graph.adj_edge[pos]]) {
        result.kept_edges.emplace_back(representative[v], representative[other]);
      }
    }
  }

  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
