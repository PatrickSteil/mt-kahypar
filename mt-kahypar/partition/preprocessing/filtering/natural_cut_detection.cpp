#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <iostream>

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace filtering {

namespace {
// mark_covered is called before the flow computation, so that concurrent
// callers do not pick core vertices as seeds.
template <typename MarkCovered>
std::vector<EdgeID> compute_natural_cut_impl(const FilterGraph& graph,
                                             NodeID seed,
                                             const NaturalCutParams& params,
                                             NaturalCutScratch& scratch,
                                             MarkCovered&& mark_covered) {
  const size_t n = graph.numNodes();
  if (scratch.in_tree.size() != n) {
    scratch.in_tree.assign(n, 0);
    scratch.in_core.assign(n, 0);
    scratch.in_ring.assign(n, 0);
    scratch.local_id.assign(n, kInvalidNode);
  } else {
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

  const auto target_tree_size =
      static_cast<NodeWeight>(params.alpha * params.U);
  const auto target_core_size =
      static_cast<NodeWeight>(params.alpha * params.U / params.f);

  NodeWeight tree_size = 0;
  size_t head = 0;
  scratch.bfs_queue.push_back(seed);
  scratch.in_tree[seed] = 1;
  scratch.tree_order.push_back(seed);
  tree_size += graph.node_weight[seed];

  while (head < scratch.bfs_queue.size() && tree_size < target_tree_size) {
    const NodeID u = scratch.bfs_queue[head++];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1];
         ++pos) {
      if (tree_size >= target_tree_size) break;
      const NodeID v = graph.adj[pos];
      if (scratch.in_tree[v]) continue;
      scratch.in_tree[v] = 1;
      scratch.tree_order.push_back(v);
      scratch.bfs_queue.push_back(v);
      tree_size += graph.node_weight[v];
    }
  }

  // The threshold is checked after adding, so the core contains at least the
  // seed. Only core vertices (not the whole tree) are marked covered: the
  // source side of the cut contains the core and is fenced off by the cut, so
  // it bounds their fragment by alpha * U. Marking tree vertices as covered
  // breaks this guarantee.
  NodeWeight running = 0;
  for (NodeID v : scratch.tree_order) {
    scratch.in_core[v] = 1;
    mark_covered(v);
    running += graph.node_weight[v];
    if (running >= target_core_size) break;
  }

  std::vector<NodeID>& ring = scratch.ring;
  std::vector<char>& in_ring = scratch.in_ring;
  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1];
         ++pos) {
      const NodeID v = graph.adj[pos];
      if (!scratch.in_tree[v] && !in_ring[v]) {
        in_ring[v] = 1;
        ring.push_back(v);
      }
    }
  }

  // Local ids: 0 = s (core), 1 = t (ring), 2.. = rest of the tree.
  std::vector<NodeID>& local_id = scratch.local_id;
  uint32_t next_local_id = 2;
  for (NodeID v : scratch.tree_order) {
    local_id[v] = scratch.in_core[v] ? 0 : next_local_id++;
  }
  for (NodeID v : ring) local_id[v] = 1;

  FlowNetwork& network = scratch.network;
  network.reset(next_local_id);
  network.source = 0;
  network.sink = 1;

  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1];
         ++pos) {
      const NodeID v = graph.adj[pos];
      assert(local_id[v] != kInvalidNode);
      if (scratch.in_tree[v] && v <= u) continue;
      if (local_id[u] == local_id[v]) continue;
      network.add_edge(local_id[u], local_id[v],
                       graph.edge_weight[graph.adj_edge[pos]]);
    }
  }

  push_relabel_max_flow(network);

  // The source side lies within the tree, so the cut edges are the edges
  // leaving it.
  std::vector<EdgeID> cut_edges;
  auto collect_cut = [&](auto&& in_source_side) {
    for (NodeID u : scratch.tree_order) {
      if (!in_source_side(local_id[u])) continue;
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1];
           ++pos) {
        if (!in_source_side(local_id[graph.adj[pos]]))
          cut_edges.push_back(graph.adj_edge[pos]);
      }
    }
  };
  if (params.cut_side != CutSide::Sink) {
    std::vector<char>& reachable = scratch.reachable;
    min_cut_reachable_from_source(network, reachable, scratch.flow_queue);
    collect_cut([&](NodeID x) { return reachable[x] != 0; });
  }
  if (params.cut_side != CutSide::Source) {
    std::vector<char>& reaches_sink = scratch.reaches_sink;
    min_cut_reaching_sink(network, reaches_sink, scratch.flow_queue);
    collect_cut([&](NodeID x) { return reaches_sink[x] == 0; });
  }
  return cut_edges;
}
}  // namespace

std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                        const NaturalCutParams& params,
                                        NaturalCutScratch& scratch,
                                        std::vector<char>& covered) {
  return compute_natural_cut_impl(graph, seed, params, scratch,
                                  [&](NodeID v) { covered[v] = 1; });
}

std::vector<char> run_natural_cut_detection_sequential(
    const FilterGraph& graph, const NaturalCutParams& params,
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
      std::vector<EdgeID> cut =
          compute_natural_cut(graph, v, params, scratch, covered);
      for (EdgeID e : cut) keep[e] = 1;
    }
  }
  return keep;
}

std::vector<char> run_natural_cut_detection(const FilterGraph& graph,
                                            const NaturalCutParams& params,
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
    const size_t count =
        max_flow_solves.fetch_add(1, std::memory_order_relaxed) + 1;
    if (count % kProgressStride == 0) {
      const double seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        pipeline_start)
              .count();
      std::cerr << "[natural-cut] " << count << " max-flow problems solved ("
                << seconds << "s elapsed)\n";
    }
  };

  for (int sweep = 0; sweep < params.coverage; ++sweep) {
    const auto sweep_start = std::chrono::steady_clock::now();
    std::vector<parallel::IntegralAtomicWrapper<uint8_t>> covered(n);
    for (size_t v = 0; v < n; ++v) covered[v] = 0;

    std::mt19937_64 shuffle_rng(0x9e3779b97f4a7c15ULL +
                                static_cast<uint64_t>(sweep));
    std::shuffle(order.begin(), order.end(), shuffle_rng);

    auto mark_covered = [&](NodeID u) {
      covered[u].store(1, std::memory_order_relaxed);
    };
    tbb::parallel_for(size_t(0), order.size(), [&](size_t i) {
      const NodeID v = order[i];
      uint8_t expected = 0;
      if (!covered[v].compare_exchange_strong(expected, 1)) return;

      std::vector<EdgeID> cut = compute_natural_cut_impl(
          graph, v, params, scratch.local(), mark_covered);
      log_solve();
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    });

    // Safety net; normally every vertex is covered at this point.
    for (NodeID v : order) {
      if (covered[v].load(std::memory_order_relaxed)) continue;
      covered[v].store(1, std::memory_order_relaxed);
      std::vector<EdgeID> cut = compute_natural_cut_impl(
          graph, v, params, scratch.local(), mark_covered);
      log_solve();
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    }

    if (verbose) {
      const double seconds = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - sweep_start)
                                 .count();
      std::cerr << "[natural-cut] sweep " << (sweep + 1) << "/"
                << params.coverage << " done (" << seconds << "s)\n";
    }
  }

  if (verbose) {
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      pipeline_start)
            .count();
    std::cerr << "[natural-cut] total: "
              << max_flow_solves.load(std::memory_order_relaxed)
              << " max-flow problems solved (" << seconds << "s)\n";
  }

  std::vector<char> result(m);
  for (size_t e = 0; e < m; ++e)
    result[e] = keep[e].load(std::memory_order_relaxed);
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
