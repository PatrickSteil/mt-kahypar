// mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.cpp
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

#include <deque>
#include <unordered_map>

#include <tbb/parallel_for.h>

#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

namespace {
std::vector<NodeID> components_from_union_find(const FilterGraph& graph, AtomicUnionFind& uf) {
  const size_t n = graph.numNodes();
  std::vector<NodeID> component(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    component[v] = uf.find(static_cast<uint32_t>(v));
  });
  return component;
}
}  // namespace

std::vector<NodeID> parallel_connected_components(const FilterGraph& graph) {
  const size_t n = graph.numNodes();
  AtomicUnionFind uf(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      uf.unite(static_cast<uint32_t>(v), static_cast<uint32_t>(graph.adj[pos]));
    }
  });
  return components_from_union_find(graph, uf);
}

std::vector<NodeID> parallel_connected_components_excluding(
    const FilterGraph& graph, const std::vector<char>& excluded_edge) {
  const size_t n = graph.numNodes();
  AtomicUnionFind uf(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      if (excluded_edge[graph.adj_edge[pos]]) continue;
      uf.unite(static_cast<uint32_t>(v), static_cast<uint32_t>(graph.adj[pos]));
    }
  });
  return components_from_union_find(graph, uf);
}

SpanningForest build_spanning_forest(const FilterGraph& graph,
                                      const std::vector<NodeID>& component,
                                      const std::vector<NodeID>& roots_in) {
  const size_t n = graph.numNodes();

  std::vector<NodeID> roots = roots_in;
  if (roots.empty()) {
    std::unordered_map<NodeID, bool> seen;
    for (size_t v = 0; v < n; ++v) {
      if (seen.emplace(component[v], true).second) roots.push_back(static_cast<NodeID>(v));
    }
  }

  SpanningForest forest;
  forest.parent.assign(n, kInvalidNode);
  forest.parent_edge.assign(n, kInvalidEdge);
  forest.roots = roots;
  std::vector<std::vector<NodeID>> order_per_root(roots.size());

  // Each root's BFS only ever touches vertices in its own component, so
  // concurrent writes to forest.parent/forest.parent_edge from different
  // roots never touch the same index -- no data race despite plain vectors.
  tbb::parallel_for(size_t(0), roots.size(), [&](size_t i) {
    const NodeID root = roots[i];
    std::vector<NodeID>& order = order_per_root[i];
    std::deque<NodeID> queue;
    forest.parent[root] = root;
    queue.push_back(root);
    order.push_back(root);
    while (!queue.empty()) {
      const NodeID u = queue.front();
      queue.pop_front();
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (forest.parent[v] == kInvalidNode) {
          forest.parent[v] = u;
          forest.parent_edge[v] = graph.adj_edge[pos];
          queue.push_back(v);
          order.push_back(v);
        }
      }
    }
  });

  forest.bfs_order.reserve(n);
  for (const std::vector<NodeID>& order : order_per_root) {
    forest.bfs_order.insert(forest.bfs_order.end(), order.begin(), order.end());
  }
  return forest;
}

}  // namespace filtering
}  // namespace mt_kahypar
