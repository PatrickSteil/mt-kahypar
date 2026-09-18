// mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

#include <unordered_map>

#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

ContractionResult contract_component_tree(const FilterGraph& graph, const TinyCutParams& params) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();

  // Step 1: bridges, via the shared spanning-tree XOR-signature technique.
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  EdgeSignatures sigs = compute_edge_signatures(graph, forest);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);

  // Step 2: blocks = connected components of the graph with bridges removed.
  std::vector<NodeID> block_rep = parallel_connected_components_excluding(graph, is_bridge);
  std::unordered_map<NodeID, NodeID> dense_block_id;
  std::vector<NodeID> block_of(n);
  for (size_t v = 0; v < n; ++v) {
    auto it = dense_block_id.find(block_rep[v]);
    if (it == dense_block_id.end()) {
      const NodeID id = static_cast<NodeID>(dense_block_id.size());
      dense_block_id.emplace(block_rep[v], id);
      block_of[v] = id;
    } else {
      block_of[v] = it->second;
    }
  }
  const NodeID num_blocks = static_cast<NodeID>(dense_block_id.size());

  std::vector<NodeWeight> block_weight(num_blocks, 0);
  for (size_t v = 0; v < n; ++v) block_weight[block_of[v]] += graph.node_weight[v];

  // Step 3: the block quotient graph. Its edges are exactly the bridges
  // (each bridge connects exactly two distinct blocks); this quotient graph
  // is a forest -- one tree per connected component of the original graph --
  // matching the design spec's "component tree T".
  std::vector<std::pair<NodeID, NodeID>> edge_endpoints = compute_edge_endpoints(graph);
  std::vector<EdgeListEntry> quotient_edges;
  for (EdgeID e = 0; e < static_cast<EdgeID>(m); ++e) {
    if (!is_bridge[e]) continue;
    quotient_edges.push_back(EdgeListEntry{
        block_of[edge_endpoints[e].first], block_of[edge_endpoints[e].second], graph.edge_weight[e]});
  }
  FilterGraph quotient = build_csr_from_edge_list(quotient_edges, block_weight);

  // Step 4: root each quotient-forest component at its heaviest block
  // ("the edge-connected component with maximum size", design spec 4.2).
  std::vector<NodeID> quotient_component = parallel_connected_components(quotient);
  std::unordered_map<NodeID, NodeID> best_root_for_component;
  std::unordered_map<NodeID, NodeWeight> best_weight_for_component;
  for (NodeID b = 0; b < num_blocks; ++b) {
    const NodeID comp = quotient_component[b];
    auto it = best_weight_for_component.find(comp);
    if (it == best_weight_for_component.end() || quotient.node_weight[b] > it->second) {
      best_weight_for_component[comp] = quotient.node_weight[b];
      best_root_for_component[comp] = b;
    }
  }
  std::vector<NodeID> roots;
  roots.reserve(best_root_for_component.size());
  for (auto& [comp, root] : best_root_for_component) roots.push_back(root);
  SpanningForest quotient_forest = build_spanning_forest(quotient, quotient_component, roots);

  // Step 5: subtree weights over the quotient tree, bottom-up via reverse
  // BFS order (BFS order is non-decreasing in depth, so every child is
  // folded into its parent before the parent itself is used).
  std::vector<NodeWeight> subtree_weight(num_blocks);
  for (NodeID b = 0; b < num_blocks; ++b) subtree_weight[b] = quotient.node_weight[b];
  for (auto it = quotient_forest.bfs_order.rbegin(); it != quotient_forest.bfs_order.rend(); ++it) {
    const NodeID b = *it;
    const NodeID parent = quotient_forest.parent[b];
    if (parent != b) subtree_weight[parent] += subtree_weight[b];
  }

  // Step 6: one representative original vertex per block, computed
  // unconditionally (a tau-merge target block may never itself be "chosen"
  // -- see below -- so every block needs a representative available).
  std::vector<NodeID> block_representative_vertex(num_blocks, kInvalidNode);
  for (size_t v = 0; v < n; ++v) {
    NodeID& rep = block_representative_vertex[block_of[v]];
    if (rep == kInvalidNode) rep = static_cast<NodeID>(v);
  }

  // Step 7: top-down selection of subtrees to contract (weight <= U), plus
  // the tau-merge into the parent block, building the vertex-level
  // union-find DIRECTLY (not via an intermediate block-level union-find --
  // see the note below on why that translation is unsound).
  //
  // `is_within_chosen_subtree[b]` marks a block whose own subtree_weight
  // triggered contraction, either directly or by inheriting from an
  // already-chosen ancestor (forward BFS order visits parents before
  // children, so this propagates correctly downward). A tau-merge target
  // (the immediate parent of a chosen small subtree) is unioned into the
  // vertex group directly, WITHOUT setting is_within_chosen_subtree on the
  // parent itself: doing so would incorrectly make every OTHER child of that
  // parent auto-absorb regardless of its own size, which is not what a
  // tau-merge means (it is a narrow, single-subtree-into-its-immediate-
  // parent extension, not "the parent is now fully contracted").
  //
  // `tau_merged_extra_weight[parent]` guards against a genuine correctness
  // bug that a naive per-merge check misses: TWO SEPARATE small subtrees
  // hanging off the SAME parent can each individually satisfy
  // "subtree + parent's own weight <= U", yet their COMBINED weight once
  // both are merged into the same parent can exceed U (e.g. parent weight 1,
  // two children of weight 4 each, U = 6: 4+1 <= 6 passes twice, but
  // 4+4+1 = 9 > 6). Tracking how much extra weight has already been folded
  // into a given parent via prior tau-merges, and including it in the check,
  // prevents this -- since Part 1's own U-cap is part of what keeps Part 2's
  // final fragments within U (see design spec section 6), not merely Part
  // 2's alpha <= 1 guarantee on its own.
  std::vector<char> is_within_chosen_subtree(num_blocks, 0);
  std::vector<NodeWeight> tau_merged_extra_weight(num_blocks, 0);
  AtomicUnionFind vertex_groups(n);
  for (NodeID b : quotient_forest.bfs_order) {
    const NodeID parent = quotient_forest.parent[b];
    if (parent != b && is_within_chosen_subtree[parent]) {
      is_within_chosen_subtree[b] = 1;
      vertex_groups.unite(block_representative_vertex[b], block_representative_vertex[parent]);
      continue;
    }
    if (subtree_weight[b] <= params.U) {
      is_within_chosen_subtree[b] = 1;
      if (parent != b && subtree_weight[b] <= params.tau &&
          subtree_weight[b] + quotient.node_weight[parent] + tau_merged_extra_weight[parent] <= params.U) {
        vertex_groups.unite(block_representative_vertex[b], block_representative_vertex[parent]);
        tau_merged_extra_weight[parent] += subtree_weight[b];
      }
    }
  }

  // Step 8: collapse every vertex of a chosen block into that block's own
  // representative (blocks that were never chosen, and never tau-merge
  // targets, keep their original internal vertices and edges untouched --
  // this is intentional: a block too big to contract must retain its
  // internal structure for pass 2/3 to potentially reduce further).
  for (size_t v = 0; v < n; ++v) {
    const NodeID b = block_of[v];
    if (is_within_chosen_subtree[b]) {
      vertex_groups.unite(block_representative_vertex[b], static_cast<NodeID>(v));
    }
  }

  return contract_graph(graph, vertex_groups);
}

namespace {
std::vector<NodeID> walk_chain_from(const FilterGraph& graph, const std::vector<char>& is_deg2,
                                     std::vector<char>& visited, NodeID prev, NodeID cur) {
  std::vector<NodeID> chain;
  while (is_deg2[cur] && !visited[cur]) {
    visited[cur] = 1;
    chain.push_back(cur);
    NodeID next = kInvalidNode;
    for (EdgeID pos = graph.node_begin[cur]; pos < graph.node_begin[cur + 1]; ++pos) {
      if (graph.adj[pos] != prev) { next = graph.adj[pos]; break; }
    }
    prev = cur;
    cur = next;
  }
  return chain;
}
}  // namespace

ContractionResult contract_degree2_chains(const FilterGraph& graph, NodeWeight U) {
  const size_t n = graph.numNodes();
  std::vector<char> is_deg2(n);
  for (size_t v = 0; v < n; ++v) is_deg2[v] = (graph.degree(static_cast<NodeID>(v)) == 2);

  std::vector<char> visited(n, 0);
  AtomicUnionFind uf(n);

  auto try_contract_chain = [&](const std::vector<NodeID>& chain) {
    if (chain.empty()) return;
    NodeWeight total = 0;
    for (NodeID v : chain) total += graph.node_weight[v];
    if (total <= U) {
      for (size_t i = 1; i < chain.size(); ++i) uf.unite(chain[0], chain[i]);
    }
  };

  // Chains anchored at both ends by a vertex of degree != 2.
  for (NodeID anchor = 0; anchor < static_cast<NodeID>(n); ++anchor) {
    if (is_deg2[anchor]) continue;
    for (EdgeID pos = graph.node_begin[anchor]; pos < graph.node_begin[anchor + 1]; ++pos) {
      const NodeID cur = graph.adj[pos];
      if (!is_deg2[cur] || visited[cur]) continue;
      try_contract_chain(walk_chain_from(graph, is_deg2, visited, anchor, cur));
    }
  }

  // Pure cycles: components with no anchor at all.
  for (NodeID v = 0; v < static_cast<NodeID>(n); ++v) {
    if (!is_deg2[v] || visited[v]) continue;
    try_contract_chain(walk_chain_from(graph, is_deg2, visited, kInvalidNode, v));
  }

  return contract_graph(graph, uf);
}

ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U) {
  const size_t n = graph.numNodes();
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  EdgeSignatures sigs = compute_edge_signatures(graph, forest);
  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(graph, sigs);

  AtomicUnionFind uf(n);
  for (const std::vector<EdgeID>& cls : classes) {
    std::vector<char> excluded(graph.numEdges(), 0);
    for (EdgeID e : cls) excluded[e] = 1;
    std::vector<NodeID> comp = parallel_connected_components_excluding(graph, excluded);

    std::unordered_map<NodeID, NodeWeight> comp_weight;
    for (size_t v = 0; v < n; ++v) comp_weight[comp[v]] += graph.node_weight[v];

    for (size_t v = 0; v < n; ++v) {
      if (comp_weight[comp[v]] <= U) uf.unite(static_cast<NodeID>(v), comp[v]);
    }
  }
  return contract_graph(graph, uf);
}

}  // namespace filtering
}  // namespace mt_kahypar
