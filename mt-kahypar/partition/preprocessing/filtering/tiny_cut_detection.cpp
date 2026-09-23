// mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

#include <unordered_map>

#include <tbb/parallel_for.h>

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

BoundedComponentsResult compute_bounded_components_excluding(
    const FilterGraph& graph, const std::vector<EdgeID>& excluded_edges,
    const std::vector<NodeID>& seeds, NodeWeight total_graph_weight, NodeWeight U) {
  const size_t n = graph.numNodes();

  // excluded_edges is always small (a 2-edge-cut class -- 2 in the
  // overwhelming common case), so a linear scan per adjacency-list visit is
  // cheaper than building/allocating an O(numEdges()) bitmap per class.
  auto is_excluded = [&](EdgeID e) {
    for (EdgeID x : excluded_edges) {
      if (x == e) return true;
    }
    return false;
  };

  std::vector<NodeID> owner(n, kInvalidNode);  // which seed-group first claimed v
  std::vector<NodeID> all_visited;             // every claimed vertex, for final bucketing
  std::vector<std::vector<NodeID>> queue;      // per seed-group BFS queue
  std::vector<size_t> cursor;                  // per seed-group read cursor into `queue`
  std::vector<NodeWeight> group_weight;

  for (NodeID s : seeds) {
    if (owner[s] != kInvalidNode) continue;  // duplicate/shared seed vertex
    const NodeID gid = static_cast<NodeID>(queue.size());
    owner[s] = gid;
    all_visited.push_back(s);
    queue.push_back({s});
    cursor.push_back(0);
    group_weight.push_back(graph.node_weight[s]);
  }
  const size_t num_groups = queue.size();

  BoundedComponentsResult result;
  if (num_groups <= 1) {
    // Degenerate: every edge endpoint collapsed to a single seed. The
    // structural guarantee above says this shouldn't happen for a genuine
    // class of size >= 2, but if it ever does, conservatively report no
    // explicit components (the whole graph is the "leftover") rather than
    // risk an out-of-bounds round-robin below.
    result.leftover_weight = total_graph_weight;
    return result;
  }

  AtomicUnionFind group_uf(num_groups);
  std::vector<char> retired(num_groups, 0);
  size_t active_count = num_groups;
  size_t rr = 0;

  while (active_count > 1) {
    NodeID g = group_uf.find(static_cast<NodeID>(rr % num_groups));
    if (retired[g]) {
      ++rr;
      continue;
    }
    if (cursor[g] >= queue[g].size()) {
      retired[g] = 1;
      --active_count;
      ++rr;
      continue;
    }
    const NodeID u = queue[g][cursor[g]++];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const EdgeID e = graph.adj_edge[pos];
      if (is_excluded(e)) continue;
      const NodeID v = graph.adj[pos];
      if (owner[v] == kInvalidNode) {
        owner[v] = g;
        all_visited.push_back(v);
        queue[g].push_back(v);
        group_weight[g] += graph.node_weight[v];
      } else {
        const NodeID other = group_uf.find(owner[v]);
        if (other != g) {
          group_uf.unite(g, other);
          const NodeID surviving = group_uf.find(g);
          const NodeID folded = (surviving == g) ? other : g;
          // The folded frontier may still have unexplored work; splice its
          // remaining (unprocessed) queue tail into the survivor's queue so
          // the merged group's traversal continues correctly. Losing this
          // would silently abandon part of the merged component, leaving
          // some of its vertices unclaimed (or worse, later misclaimed by
          // a different still-active group).
          for (size_t idx = cursor[folded]; idx < queue[folded].size(); ++idx) {
            queue[surviving].push_back(queue[folded][idx]);
          }
          queue[folded].clear();
          cursor[folded] = 0;
          group_weight[surviving] += group_weight[folded];
          group_weight[folded] = 0;
          if (!retired[folded]) {
            retired[folded] = 1;
            --active_count;
          }
          // `g` may itself have been the folded side; re-point it at the
          // surviving id so any further neighbors of `u` processed in the
          // REST of this same adjacency-list scan are still attributed to
          // the group that actually still exists (queue[g]/group_weight[g]
          // were just cleared/zeroed above if g == folded).
          g = surviving;
          if (active_count <= 1) break;
        }
        // else: same group already, nothing to do.
      }
    }
    ++rr;
  }

  // Exactly one group remains un-retired: the one whose full extent is
  // inferred by elimination rather than explored (it is never correct to
  // report it as an explicit component -- we may only have SOME of its
  // members from partial exploration before it became the sole survivor,
  // never all of them).
  NodeID leftover_root = kInvalidNode;
  for (size_t i = 0; i < num_groups; ++i) {
    const NodeID gid = static_cast<NodeID>(i);
    if (!retired[gid] && group_uf.find(gid) == gid) {
      leftover_root = gid;
      break;
    }
  }

  // Bucket every claimed vertex (bounded by all_visited.size(), i.e. ~2x the
  // total size of the non-largest components -- never the full graph) by its
  // final group representative, excluding the leftover group's own bucket.
  std::unordered_map<NodeID, std::vector<NodeID>> buckets;
  for (NodeID v : all_visited) {
    const NodeID root = group_uf.find(owner[v]);
    if (root == leftover_root) continue;
    buckets[root].push_back(v);
  }

  NodeWeight explored_weight = 0;
  for (auto& [rep, members] : buckets) {
    const NodeWeight w = group_weight[rep];
    explored_weight += w;
    result.small_component_members.push_back(std::move(members));
    result.small_component_weight.push_back(w);
  }
  result.leftover_weight = total_graph_weight - explored_weight;

  if (result.leftover_weight > 0 && result.leftover_weight <= U) {
    // Rare: even the (typically huge) leftover component qualifies for
    // contraction. We never fully explored it, so materialize its exact
    // membership with one O(n) scan -- acceptable since this only triggers
    // when U is large relative to the graph, not in the hot path this
    // function exists to avoid. A vertex belongs to the leftover iff it was
    // never claimed by an EXPLICIT (non-leftover) group -- this correctly
    // includes both truly-unvisited vertices and the leftover group's own
    // partially-explored members.
    std::vector<char> is_explicit(n, 0);
    for (NodeID v : all_visited) {
      if (group_uf.find(owner[v]) != leftover_root) is_explicit[v] = 1;
    }
    for (NodeID v = 0; v < static_cast<NodeID>(n); ++v) {
      if (!is_explicit[v]) result.leftover_members.push_back(v);
    }
  }
  return result;
}

ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U) {
  const size_t n = graph.numNodes();
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  EdgeSignatures sigs = compute_edge_signatures(graph, forest);
  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(sigs);

  // Computed once, up front, rather than per class.
  std::vector<std::pair<NodeID, NodeID>> edge_endpoints = compute_edge_endpoints(graph);
  NodeWeight total_graph_weight = 0;
  for (size_t v = 0; v < n; ++v) total_graph_weight += graph.node_weight[v];

  AtomicUnionFind uf(n);
  // Each class's processing reads only shared IMMUTABLE state (`graph`,
  // `classes`, `edge_endpoints`) and writes only to `uf`, whose unite() is
  // already safe under concurrent calls from multiple threads (see
  // union_find.h) -- so classes can be processed fully in parallel. Per
  // class, compute_bounded_components_excluding bounds work to ~2x the
  // total size of the non-largest resulting components (the paper's
  // "two-at-a-time" traversal, design spec section 4.4 step 8) instead of
  // a full O(n) scan -- fixing a real-world-validation finding (Task 19,
  // re-profiled after an earlier fix that only parallelized the old O(n)-
  // per-class scan): on real DIMACS road networks, tens of thousands of
  // classes each doing a full graph scan dominated 92-98% of total
  // pipeline runtime, independent of parallelization, since the total work
  // was O(classes * n) regardless of how many threads ran it concurrently.
  tbb::parallel_for(size_t(0), classes.size(), [&](size_t i) {
    const std::vector<EdgeID>& cls = classes[i];
    std::vector<NodeID> seeds;
    seeds.reserve(cls.size() * 2);
    for (EdgeID e : cls) {
      seeds.push_back(edge_endpoints[e].first);
      seeds.push_back(edge_endpoints[e].second);
    }

    BoundedComponentsResult res =
        compute_bounded_components_excluding(graph, cls, seeds, total_graph_weight, U);

    for (size_t k = 0; k < res.small_component_members.size(); ++k) {
      if (res.small_component_weight[k] <= U) {
        const std::vector<NodeID>& members = res.small_component_members[k];
        for (size_t idx = 1; idx < members.size(); ++idx) uf.unite(members[0], members[idx]);
      }
    }
    if (!res.leftover_members.empty()) {
      for (size_t idx = 1; idx < res.leftover_members.size(); ++idx) {
        uf.unite(res.leftover_members[0], res.leftover_members[idx]);
      }
    }
  });
  return contract_graph(graph, uf);
}

ContractionResult run_tiny_cut_detection(const FilterGraph& graph, const TinyCutParams& params) {
  ContractionResult r1 = contract_component_tree(graph, params);
  ContractionResult r2 = contract_degree2_chains(r1.graph, params.U);
  ContractionResult r3 = contract_two_edge_cuts(r2.graph, params.U);

  std::vector<NodeID> composed(graph.numNodes());
  for (size_t v = 0; v < graph.numNodes(); ++v) {
    composed[v] = r3.mapping[r2.mapping[r1.mapping[v]]];
  }

  ContractionResult result;
  result.graph = std::move(r3.graph);
  result.mapping = std::move(composed);
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
