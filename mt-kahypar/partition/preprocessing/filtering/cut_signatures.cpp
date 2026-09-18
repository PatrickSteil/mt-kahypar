#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"

#include <map>
#include <random>
#include <set>
#include <unordered_map>

#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

namespace {
unsigned __int128 random_128() {
  static thread_local std::mt19937_64 rng(std::random_device{}());
  const uint64_t hi = rng();
  const uint64_t lo = rng();
  return (static_cast<unsigned __int128>(hi) << 64) | lo;
}
}  // namespace

EdgeSignatures compute_edge_signatures(const FilterGraph& graph, const SpanningForest& forest) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();

  EdgeSignatures sigs;
  sigs.signature.assign(m, 0);
  sigs.is_tree_edge.assign(m, 0);

  for (size_t v = 0; v < n; ++v) {
    if (forest.parent_edge[v] != kInvalidEdge) sigs.is_tree_edge[forest.parent_edge[v]] = 1;
  }

  // diff[v] will hold the XOR of labels of non-tree edges whose fundamental
  // cycle's tree-path endpoint is exactly v (i.e. one endpoint of the
  // non-tree edge); accumulating diff bottom-up along the tree turns this
  // into, for each tree edge (parent(v), v), the XOR over all non-tree edges
  // whose path covers that edge.
  std::vector<unsigned __int128> diff(n, 0);

  for (size_t v = 0; v < n; ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const EdgeID e = graph.adj_edge[pos];
      const NodeID other = graph.adj[pos];
      if (sigs.is_tree_edge[e]) continue;
      if (other <= static_cast<NodeID>(v)) continue;  // process each non-tree edge once
      const unsigned __int128 label = random_128();
      sigs.signature[e] = label;  // non-tree edge's own signature is its label
      diff[v] ^= label;
      diff[other] ^= label;
    }
  }

  // Accumulate diff into per-tree-edge signatures via reverse BFS order:
  // BFS order is non-decreasing in depth, so processing in reverse guarantees
  // every child is folded into its parent before the parent itself is used.
  for (auto it = forest.bfs_order.rbegin(); it != forest.bfs_order.rend(); ++it) {
    const NodeID v = *it;
    const NodeID parent = forest.parent[v];
    if (parent == v) continue;  // root
    sigs.signature[forest.parent_edge[v]] = diff[v];
    diff[parent] ^= diff[v];
  }

  return sigs;
}

std::vector<char> compute_bridges(const FilterGraph& graph, const EdgeSignatures& sigs) {
  std::vector<char> is_bridge(graph.numEdges(), 0);
  for (EdgeID e = 0; e < static_cast<EdgeID>(graph.numEdges()); ++e) {
    is_bridge[e] = sigs.is_tree_edge[e] && sigs.signature[e] == 0;
  }
  return is_bridge;
}

namespace {
struct SignatureHash {
  size_t operator()(unsigned __int128 x) const {
    const uint64_t hi = static_cast<uint64_t>(x >> 64);
    const uint64_t lo = static_cast<uint64_t>(x);
    return std::hash<uint64_t>()(hi) ^ (std::hash<uint64_t>()(lo) * 0x9e3779b97f4a7c15ULL);
  }
};

// A class is verified if:
// 1. Removing all its edges actually disconnects at least one pair of incident
//    vertices from each other.
// 2. The class does NOT form a simple cycle within a larger component (such a
//    cycle is not a genuine 2-cut, but a degenerate case of being "trapped" in
//    a local cycle).
bool verify_class_disconnects(const FilterGraph& graph, const std::vector<EdgeID>& cls) {
  std::vector<char> excluded(graph.numEdges(), 0);
  for (EdgeID e : cls) excluded[e] = 1;
  std::vector<NodeID> component = parallel_connected_components_excluding(graph, excluded);
  auto endpoints = compute_edge_endpoints(graph);
  const NodeID u = endpoints[cls[0]].first;
  const NodeID v = endpoints[cls[0]].second;

  // Check basic disconnection
  if (component[u] == component[v]) return false;

  // Check that the class doesn't form a simple cycle (every vertex in the class
  // has degree exactly 2 within the class). A simple cycle that is a proper
  // subset of a larger component is degenerate and should be rejected.
  std::set<NodeID> class_vertices;
  std::map<NodeID, int> degree_in_class;
  for (EdgeID e : cls) {
    auto [s, t] = endpoints[e];
    class_vertices.insert(s);
    class_vertices.insert(t);
    degree_in_class[s]++;
    degree_in_class[t]++;
  }

  bool is_simple_cycle = true;
  for (NodeID v : class_vertices) {
    if (degree_in_class[v] != 2) {
      is_simple_cycle = false;
      break;
    }
  }

  // If it's a simple cycle, check if it's a proper subset of a larger component
  // by checking if there are edges from cycle vertices to non-cycle vertices.
  if (is_simple_cycle) {
    for (NodeID v : class_vertices) {
      for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
        const EdgeID e = graph.adj_edge[pos];
        const NodeID other = graph.adj[pos];
        // Skip edges that are in the class (already excluded)
        if (excluded[e]) continue;
        // If this edge connects to a vertex outside the cycle, reject the class
        if (class_vertices.find(other) == class_vertices.end()) {
          return false;  // Cycle is a proper subset of a larger component
        }
      }
    }
    // If we reach here, the cycle is not connected to anything outside, so it's
    // the entire component and is valid
  }

  return true;
}
}  // namespace

std::vector<std::vector<EdgeID>> find_two_edge_cut_classes(const FilterGraph& graph,
                                                            const EdgeSignatures& sigs) {
  std::unordered_map<unsigned __int128, std::vector<EdgeID>, SignatureHash> buckets;
  for (EdgeID e = 0; e < static_cast<EdgeID>(graph.numEdges()); ++e) {
    const bool is_bridge = sigs.is_tree_edge[e] && sigs.signature[e] == 0;
    if (is_bridge) continue;
    buckets[sigs.signature[e]].push_back(e);
  }

  std::vector<std::vector<EdgeID>> classes;
  for (auto& [signature, edges] : buckets) {
    if (edges.size() < 2) continue;
    if (verify_class_disconnects(graph, edges)) classes.push_back(std::move(edges));
  }
  return classes;
}

}  // namespace filtering
}  // namespace mt_kahypar
