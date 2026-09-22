// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                         const NaturalCutParams& params,
                                         NaturalCutScratch& scratch,
                                         std::vector<char>& covered) {
  const size_t n = graph.numNodes();
  scratch.in_tree.assign(n, 0);
  scratch.in_core.assign(n, 0);
  scratch.tree_order.clear();
  scratch.bfs_queue.clear();

  const auto target_tree_size = static_cast<NodeWeight>(params.alpha * params.U);
  const auto target_core_size = static_cast<NodeWeight>(params.alpha * params.U / params.f);

  NodeWeight tree_size = 0;
  size_t head = 0;
  scratch.bfs_queue.push_back(seed);
  scratch.in_tree[seed] = 1;
  scratch.tree_order.push_back(seed);
  tree_size += graph.node_weight[seed];
  covered[seed] = 1;

  while (head < scratch.bfs_queue.size() && tree_size < target_tree_size) {
    const NodeID u = scratch.bfs_queue[head++];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      if (tree_size >= target_tree_size) break;
      const NodeID v = graph.adj[pos];
      if (scratch.in_tree[v]) continue;
      scratch.in_tree[v] = 1;
      scratch.tree_order.push_back(v);
      scratch.bfs_queue.push_back(v);
      tree_size += graph.node_weight[v];
      covered[v] = 1;
    }
  }

  // Add vertices to the core, checking the threshold AFTER each add (not
  // before). This is deliberate, not an off-by-one: checking before adding
  // means the very first vertex is only admitted if target_core_size > 0,
  // so a small U/f combination that floors target_core_size to 0 (e.g.
  // U=3, f=10 gives alpha*U/f=0.3, truncated to 0) would leave core
  // permanently empty -- the seed itself never gets in -- which then
  // starves the local flow network's source of any identity. Checking
  // after admitting guarantees core always contains at least the seed,
  // regardless of how small target_core_size rounds down to.
  NodeWeight running = 0;
  for (NodeID v : scratch.tree_order) {
    scratch.in_core[v] = 1;
    running += graph.node_weight[v];
    if (running >= target_core_size) break;
  }

  std::vector<NodeID> ring;
  std::vector<char> in_ring(n, 0);
  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (!scratch.in_tree[v] && !in_ring[v]) {
        in_ring[v] = 1;
        ring.push_back(v);
      }
    }
  }

  // local id 0 = s (core), local id 1 = t (ring), 2.. = tree-interior.
  std::vector<NodeID> local_id(n, kInvalidNode);
  uint32_t next_local_id = 2;
  for (NodeID v : scratch.tree_order) {
    local_id[v] = scratch.in_core[v] ? 0 : next_local_id++;
  }
  for (NodeID v : ring) local_id[v] = 1;

  FlowNetwork network(next_local_id);
  network.source = 0;
  network.sink = 1;
  std::vector<std::vector<EdgeID>> local_arc_to_original(next_local_id);

  auto add_local_edge = [&](NodeID lu, NodeID lv, EdgeID original_edge, EdgeWeight weight) {
    if (lu == lv) return;  // both endpoints collapsed to the same local node
    network.add_edge(lu, lv, weight);
    local_arc_to_original[lu].push_back(original_edge);
    local_arc_to_original[lv].push_back(original_edge);
  };

  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (local_id[v] == kInvalidNode) continue;
      if (scratch.in_tree[v] && v <= u) continue;  // dedup tree-tree edges
      add_local_edge(local_id[u], local_id[v], graph.adj_edge[pos], graph.edge_weight[graph.adj_edge[pos]]);
    }
  }
  for (NodeID u : ring) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (local_id[v] == kInvalidNode || local_id[v] == 1) continue;  // outside net, or ring-ring
      if (scratch.in_tree[v]) continue;  // handled from the tree side above
      add_local_edge(local_id[u], local_id[v], graph.adj_edge[pos], graph.edge_weight[graph.adj_edge[pos]]);
    }
  }

  dinic_max_flow(network);
  std::vector<char> reachable = min_cut_reachable_from_source(network);

  std::vector<EdgeID> cut_edges;
  for (uint32_t lu = 0; lu < network.numNodes(); ++lu) {
    if (!reachable[lu]) continue;
    for (size_t i = 0; i < network.adj[lu].size(); ++i) {
      if (!reachable[network.adj[lu][i].to]) cut_edges.push_back(local_arc_to_original[lu][i]);
    }
  }
  return cut_edges;
}

}  // namespace filtering
}  // namespace mt_kahypar
