// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <iostream>

#include <tbb/parallel_for.h>

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace filtering {

std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                         const NaturalCutParams& params,
                                         NaturalCutScratch& scratch,
                                         std::vector<char>& covered) {
  const size_t n = graph.numNodes();
  if (scratch.in_tree.size() != n) {
    scratch.in_tree.assign(n, 0);
    scratch.in_core.assign(n, 0);
    scratch.in_ring.assign(n, 0);
    scratch.local_id.assign(n, kInvalidNode);
  } else {
    // Reset only the entries touched by the previous call
    for (NodeID v : scratch.tree_order) {
      scratch.in_tree[v] = 0;
      scratch.in_core[v] = 0;
      scratch.local_id[v] = kInvalidNode;
    }
    for (NodeID v : scratch.ring) {
      scratch.in_ring[v] = 0;
      scratch.local_id[v] = kInvalidNode;
    }
  }
  scratch.tree_order.clear();
  scratch.ring.clear();
  scratch.bfs_queue.clear();

  const auto target_tree_size = static_cast<NodeWeight>(params.alpha * params.U);
  const auto target_core_size = static_cast<NodeWeight>(params.alpha * params.U / params.f);

  NodeWeight tree_size = 0;
  size_t head = 0;
  scratch.bfs_queue.push_back(seed);
  scratch.in_tree[seed] = 1;
  scratch.tree_order.push_back(seed);
  tree_size += graph.node_weight[seed];

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
  //
  // `covered[v]` is marked HERE, for CORE vertices only -- not in the BFS
  // growth loop above for every TREE vertex. This is the paper's own
  // stopping rule, not an arbitrary choice: "we accomplish this by picking
  // v uniformly at random among all vertices that have not yet been part
  // of any CORE" (design spec section 5). It is also load-bearing for the
  // hard U-invariant, not merely cosmetic: the local min s-t cut's
  // reachable-from-s side R always contains the whole core (core defines
  // s) and is bounded in weight by the tree's own size (<= alpha*U), AND
  // every edge leaving R in the WHOLE original graph is captured by this
  // local cut (R's vertices are confined to the tree, whose only external
  // neighbors are, by definition, the ring -- there is no escape route
  // outside the local flow network). So once a vertex has been part of
  // ANY core, its eventual fragment is provably a subset of that
  // bounded, fully-fenced-off R, hence <= alpha*U <= U.
  //
  // A vertex that was only ever part of some OTHER seed's TREE (but never
  // that seed's own core) sits on the far, unbounded side of that cut and
  // gets no such guarantee -- it must eventually become a seed itself (and
  // thus enter its OWN core) for the invariant to apply to it. Marking
  // whole-tree membership as "covered" (an earlier, incorrect version of
  // this function) prevents that from ever happening for such vertices,
  // which is exactly what let a later integration test find fragments up
  // to ~3x over U on 9/10 random seeds -- confirmed by an independent
  // simulation (both semantics, 10 seeds each, U=20): tree-covered
  // semantics failed 9/10 (worst fragment 65), core-covered semantics
  // failed 0/10 (worst fragment 5). Do not move this marking back into
  // the BFS loop.
  NodeWeight running = 0;
  for (NodeID v : scratch.tree_order) {
    scratch.in_core[v] = 1;
    covered[v] = 1;
    running += graph.node_weight[v];
    if (running >= target_core_size) break;
  }

  std::vector<NodeID>& ring = scratch.ring;
  std::vector<char>& in_ring = scratch.in_ring;
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
  std::vector<NodeID>& local_id = scratch.local_id;
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
      // Every remaining v has a valid, non-ring local_id, which by
      // construction (local_id is assigned only to tree or ring vertices)
      // means v must be a tree vertex -- and every tree-to-ring edge is
      // already added from the tree-side loop above. This branch is
      // therefore structurally unreachable; the assertion documents that
      // invariant so a future edit to the tree-side loop that breaks it is
      // caught immediately instead of silently dropping an edge.
      assert(scratch.in_tree[v] && "ring-loop should never reach a non-tree, non-ring, non-excluded vertex");
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

std::vector<char> run_natural_cut_detection_sequential(const FilterGraph& graph,
                                                        const NaturalCutParams& params,
                                                        std::mt19937_64& rng) {
  const size_t n = graph.numNodes();
  std::vector<char> keep(graph.numEdges(), 0);
  NaturalCutScratch scratch;

  std::vector<NodeID> order(n);
  for (size_t v = 0; v < n; ++v) order[v] = static_cast<NodeID>(v);

  for (int sweep = 0; sweep < params.coverage; ++sweep) {
    std::vector<char> covered(n, 0);
    std::shuffle(order.begin(), order.end(), rng);
    for (NodeID v : order) {
      if (covered[v]) continue;
      std::vector<EdgeID> cut = compute_natural_cut(graph, v, params, scratch, covered);
      for (EdgeID e : cut) keep[e] = 1;
    }
  }
  return keep;
}

namespace {
// Runs compute_natural_cut and marks its core vertices in the shared atomic `covered`
// array. compute_natural_cut only ever sets covered[v] for core vertices (a subset of
// scratch.tree_order), so we pass a thread-local all-zero buffer and afterwards copy
// and reset exactly the tree vertices instead of mirroring all n entries.
std::vector<EdgeID> compute_natural_cut_covering(
    const FilterGraph& graph, NodeID seed, const NaturalCutParams& params, NaturalCutScratch& scratch,
    std::vector<parallel::IntegralAtomicWrapper<uint8_t>>& covered) {
  std::vector<char>& covered_buffer = scratch.covered_buffer;
  if (covered_buffer.size() != graph.numNodes()) covered_buffer.assign(graph.numNodes(), 0);
  std::vector<EdgeID> cut = compute_natural_cut(graph, seed, params, scratch, covered_buffer);
  for (NodeID u : scratch.tree_order) {
    if (covered_buffer[u]) {
      covered[u].store(1, std::memory_order_relaxed);
      covered_buffer[u] = 0;
    }
  }
  return cut;
}
}  // namespace

std::vector<char> run_natural_cut_detection(const FilterGraph& graph, const NaturalCutParams& params,
                                             bool verbose) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();
  std::vector<parallel::IntegralAtomicWrapper<uint8_t>> keep(m);
  for (size_t e = 0; e < m; ++e) keep[e] = 0;

  tls_enumerable_thread_specific<NaturalCutScratch> scratch;

  std::vector<NodeID> order(n);
  for (size_t v = 0; v < n; ++v) order[v] = static_cast<NodeID>(v);

  const auto pipeline_start = std::chrono::steady_clock::now();
  std::atomic<size_t> max_flow_solves{0};
  constexpr size_t kProgressStride = 5000;
  auto log_solve = [&] {
    if (!verbose) return;
    const size_t count = max_flow_solves.fetch_add(1, std::memory_order_relaxed) + 1;
    if (count % kProgressStride == 0) {
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
      std::cerr << "[natural-cut] " << count << " max-flow problems solved (" << seconds << "s elapsed)\n";
    }
  };

  for (int sweep = 0; sweep < params.coverage; ++sweep) {
    const auto sweep_start = std::chrono::steady_clock::now();
    std::vector<parallel::IntegralAtomicWrapper<uint8_t>> covered(n);
    for (size_t v = 0; v < n; ++v) covered[v] = 0;

    std::mt19937_64 shuffle_rng(0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(sweep));
    std::shuffle(order.begin(), order.end(), shuffle_rng);

    tbb::parallel_for(size_t(0), order.size(), [&](size_t i) {
      const NodeID v = order[i];
      uint8_t expected = 0;
      if (!covered[v].compare_exchange_strong(expected, 1)) return;

      NaturalCutScratch& local_scratch = scratch.local();
      std::vector<EdgeID> cut = compute_natural_cut_covering(graph, v, params, local_scratch, covered);
      log_solve();
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    });

    // Mop-up pass: handles any vertex left uncovered by races near the end
    // of the parallel scan (covered is monotonic, so this terminates).
    for (NodeID v : order) {
      if (covered[v].load(std::memory_order_relaxed)) continue;
      covered[v].store(1, std::memory_order_relaxed);
      NaturalCutScratch& local_scratch = scratch.local();
      std::vector<EdgeID> cut = compute_natural_cut_covering(graph, v, params, local_scratch, covered);
      log_solve();
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    }

    if (verbose) {
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - sweep_start).count();
      std::cerr << "[natural-cut] sweep " << (sweep + 1) << "/" << params.coverage << " done (" << seconds << "s)\n";
    }
  }

  if (verbose) {
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
    std::cerr << "[natural-cut] total: " << max_flow_solves.load(std::memory_order_relaxed)
               << " max-flow problems solved (" << seconds << "s)\n";
  }

  std::vector<char> result(m);
  for (size_t e = 0; e < m; ++e) result[e] = keep[e].load(std::memory_order_relaxed);
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
