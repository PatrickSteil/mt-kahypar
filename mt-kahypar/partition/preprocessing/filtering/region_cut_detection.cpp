#include "mt-kahypar/partition/preprocessing/filtering/region_cut_detection.h"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <numeric>

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace filtering {

namespace {
// Multi-source BFS from the centers. Unreached vertices get new centers.
void assign_voronoi_cells(const FilterGraph& graph, std::vector<NodeID>& centers,
                          std::vector<NodeID>& region_of) {
  const size_t n = graph.numNodes();
  region_of.assign(n, kInvalidNode);
  std::vector<NodeID> queue;
  queue.reserve(n);
  auto bfs = [&](size_t head) {
    while (head < queue.size()) {
      const NodeID u = queue[head++];
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (region_of[v] != kInvalidNode) continue;
        region_of[v] = region_of[u];
        queue.push_back(v);
      }
    }
  };
  for (NodeID i = 0; i < centers.size(); ++i) {
    region_of[centers[i]] = i;
    queue.push_back(centers[i]);
  }
  bfs(0);
  for (NodeID v = 0; v < n; ++v) {
    if (region_of[v] != kInvalidNode) continue;
    const size_t head = queue.size();
    region_of[v] = static_cast<NodeID>(centers.size());
    centers.push_back(v);
    queue.push_back(v);
    bfs(head);
  }
}

// Moves every center to the vertex of its region farthest from the boundary.
void move_centers_inward(const FilterGraph& graph, const std::vector<NodeID>& region_of,
                         std::vector<NodeID>& centers) {
  const size_t n = graph.numNodes();
  std::vector<int32_t> dist(n, -1);
  std::vector<NodeID> queue;
  queue.reserve(n);
  for (NodeID u = 0; u < n; ++u) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      if (region_of[graph.adj[pos]] != region_of[u]) {
        dist[u] = 0;
        queue.push_back(u);
        break;
      }
    }
  }
  for (size_t head = 0; head < queue.size(); ++head) {
    const NodeID u = queue[head];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (dist[v] != -1 || region_of[v] != region_of[u]) continue;
      dist[v] = dist[u] + 1;
      queue.push_back(v);
    }
  }
  std::vector<int32_t> best(centers.size(), -1);
  for (NodeID v = 0; v < n; ++v) {
    const NodeID r = region_of[v];
    if (dist[v] > best[r]) {
      best[r] = dist[v];
      centers[r] = v;
    }
  }
}

// Per-thread buffers; only the entries of the previous subgraph are reset.
struct PairScratch {
  std::vector<int32_t> dist;      // hops from the r1-r2 border
  std::vector<NodeID> local_id;
  std::vector<NodeID> subgraph;   // r1 u r2 in BFS order, then margin vertices
  std::vector<NodeID> layer;
  std::vector<NodeID> next_layer;
  std::vector<char> side;         // 0 = reached from r1, 1 = from r2
  FlowNetwork network;
  std::vector<char> reachable;
  std::vector<char> reaches_sink;
  std::vector<uint32_t> flow_queue;
  std::vector<NodeID> ball;
  std::vector<int32_t> ball_dist;
};

// Shrinks S = {v : local_id[v] == id} to a ball of `fraction` of its weight
// around the vertex of S farthest from S's boundary.
void shrink_to_ball(const FilterGraph& graph, std::vector<NodeID>& local_id,
                    const std::vector<NodeID>& subgraph, NodeID id, double fraction,
                    PairScratch& scratch) {
  std::vector<NodeID>& members = scratch.ball;
  std::vector<int32_t>& dist = scratch.ball_dist;
  if (dist.size() != graph.numNodes()) dist.assign(graph.numNodes(), -1);
  members.clear();
  NodeWeight total = 0;
  for (NodeID v : subgraph) {
    if (local_id[v] != id) continue;
    members.push_back(v);
    total += graph.node_weight[v];
  }
  if (members.empty()) return;
  auto in_set = [&](NodeID v) { return local_id[v] == id; };

  // Multi-source BFS from S's boundary; the last vertex reached is the center.
  std::vector<NodeID> queue;
  queue.reserve(members.size());
  for (NodeID u : members) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      if (!in_set(graph.adj[pos])) {
        dist[u] = 0;
        queue.push_back(u);
        break;
      }
    }
  }
  if (queue.empty()) queue.push_back(members.front()), dist[members.front()] = 0;
  for (size_t head = 0; head < queue.size(); ++head) {
    const NodeID u = queue[head];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (dist[v] != -1 || !in_set(v)) continue;
      dist[v] = dist[u] + 1;
      queue.push_back(v);
    }
  }
  const NodeID center = queue.back();
  for (NodeID v : members) dist[v] = -1;

  // BFS inside S from the center; keep the first `fraction` of S's weight.
  const NodeWeight target = std::max<NodeWeight>(1, static_cast<NodeWeight>(fraction * total));
  queue.clear();
  queue.push_back(center);
  dist[center] = 0;
  NodeWeight kept = 0;
  size_t head = 0;
  for (; head < queue.size() && kept < target; ++head) {
    const NodeID u = queue[head];
    kept += graph.node_weight[u];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (dist[v] != -1 || !in_set(v)) continue;
      dist[v] = 0;
      queue.push_back(v);
    }
  }
  for (NodeID v : members) dist[v] = -1;
  for (size_t i = 0; i < head; ++i) dist[queue[i]] = 0;
  for (NodeID v : members) {
    if (dist[v] == -1) local_id[v] = kInvalidNode;
    dist[v] = -1;
  }
}

// Solves region edge (a, b) with S1 in a. Returns the flow value.
template <typename Mark, typename Cover>
int64_t cut_region_pair(const FilterGraph& graph, const std::vector<NodeID>& region_of,
                        const std::vector<NodeWeight>& region_weight, NodeID a, NodeID b,
                        const std::vector<NodeID>& pair_vertices, const RegionCutParams& params,
                        PairScratch& scratch, Mark&& mark, Cover&& cover) {
  const size_t n = graph.numNodes();
  std::vector<int32_t>& dist = scratch.dist;
  std::vector<NodeID>& local_id = scratch.local_id;
  std::vector<NodeID>& subgraph = scratch.subgraph;
  if (dist.size() != n) {
    dist.assign(n, -1);
    local_id.assign(n, kInvalidNode);
  } else {
    for (NodeID v : subgraph) {
      dist[v] = -1;
      local_id[v] = kInvalidNode;
    }
  }
  subgraph.clear();

  // BFS inside a + b from the a-b border (Center: also from the boundary to
  // other regions).
  auto in_pair = [&](NodeID v) { return region_of[v] == a || region_of[v] == b; };
  const bool center = params.terminals == TerminalChoice::Center;
  for (NodeID u : pair_vertices) {
    const NodeID other = region_of[u] == a ? b : a;
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID r = region_of[graph.adj[pos]];
      if (r == other || (center && r != region_of[u])) {
        dist[u] = 0;
        subgraph.push_back(u);
        break;
      }
    }
  }
  for (size_t head = 0; head < subgraph.size(); ++head) {
    const NodeID u = subgraph[head];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (dist[v] != -1 || !in_pair(v)) continue;
      dist[v] = dist[u] + 1;
      subgraph.push_back(v);
    }
  }
  const size_t pair_size = subgraph.size();

  // Terminals: the tail of the BFS order of each region. Local ids: 0 = s, 1 = t.
  const NodeWeight target_a = static_cast<NodeWeight>(params.terminal_fraction * region_weight[a]);
  const NodeWeight target_b = static_cast<NodeWeight>(params.terminal_fraction * region_weight[b]);
  NodeWeight weight_a = 0, weight_b = 0;
  bool full_a = false, full_b = false;
  for (size_t i = pair_size; i-- > 0 && !(full_a && full_b);) {
    const NodeID v = subgraph[i];
    if (region_of[v] == a && !full_a) {
      local_id[v] = 0;
      weight_a += graph.node_weight[v];
      full_a = weight_a >= target_a;
    } else if (region_of[v] == b && !full_b) {
      local_id[v] = 1;
      weight_b += graph.node_weight[v];
      full_b = weight_b >= target_b;
    }
  }
  if (params.terminals == TerminalChoice::Ball) {
    shrink_to_ball(graph, local_id, subgraph, 0, params.ball_fraction, scratch);
    shrink_to_ball(graph, local_id, subgraph, 1, params.ball_fraction, scratch);
  }
  uint32_t next_local_id = 2;
  for (size_t i = 0; i < pair_size; ++i) {
    if (local_id[subgraph[i]] == kInvalidNode) local_id[subgraph[i]] = next_local_id++;
  }

  // Extend by `margin` layers; margin vertices inherit the side they were
  // reached from.
  std::vector<char>& side = scratch.side;
  if (side.size() != n) side.assign(n, 0);
  for (NodeID v : subgraph) side[v] = region_of[v] == b;
  // In closed mode, one more layer is added and contracted into t.
  const int layers = params.closed ? params.margin + 1 : params.margin;
  scratch.layer.assign(subgraph.begin(), subgraph.end());
  for (int depth = 0; depth < layers && !scratch.layer.empty(); ++depth) {
    const bool outermost = (params.pin_ring || params.closed) && depth + 1 == layers;
    scratch.next_layer.clear();
    for (NodeID u : scratch.layer) {
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (local_id[v] != kInvalidNode) continue;
        side[v] = side[u];
        local_id[v] = outermost ? (params.closed ? 1 : static_cast<NodeID>(side[u])) : next_local_id++;
        subgraph.push_back(v);
        scratch.next_layer.push_back(v);
      }
    }
    std::swap(scratch.layer, scratch.next_layer);
  }

  // Edges leaving the subgraph are dropped (open boundary).
  FlowNetwork& network = scratch.network;
  network.reset(next_local_id);
  network.source = 0;
  network.sink = 1;
  for (NodeID u : subgraph) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (v <= u || local_id[v] == kInvalidNode || local_id[u] == local_id[v]) continue;
      network.add_edge(local_id[u], local_id[v], graph.edge_weight[graph.adj_edge[pos]]);
    }
  }
  const int64_t flow = push_relabel_max_flow(network);

  auto collect_cut = [&](auto&& in_source_side) {
    for (NodeID u : subgraph) {
      if (!in_source_side(local_id[u])) continue;
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (local_id[v] != kInvalidNode && !in_source_side(local_id[v])) mark(graph.adj_edge[pos]);
      }
    }
  };
  if (params.cut_side != CutSide::Sink) {
    min_cut_reachable_from_source(network, scratch.reachable, scratch.flow_queue);
    collect_cut([&](NodeID x) { return scratch.reachable[x] != 0; });
  }
  if (params.cut_side != CutSide::Source) {
    min_cut_reaching_sink(network, scratch.reaches_sink, scratch.flow_queue);
    collect_cut([&](NodeID x) { return scratch.reaches_sink[x] == 0; });
  }

  if (params.closed) {
    NodeWeight enclosed = 0;
    for (NodeID v : subgraph)
      if (local_id[v] != 1) enclosed += graph.node_weight[v];
    if (enclosed <= params.max_enclosed_weight) {
      for (NodeID v : subgraph)
        if (local_id[v] == 0) cover(v);
    }
  }
  return flow;
}
}  // namespace

Regions compute_regions(const FilterGraph& graph, NodeID num_regions, int lloyd_iterations,
                        std::mt19937_64& rng, const std::vector<NodeID>* avoid, double spread_fraction) {
  const size_t n = graph.numNodes();
  Regions regions;
  if (n == 0) return regions;
  num_regions = std::max<NodeID>(1, std::min<NodeID>(num_regions, static_cast<NodeID>(n)));

  std::vector<NodeID> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);
  if (avoid && !avoid->empty()) {
    // Spread: draw centers among the spread_fraction of vertices farthest from
    // `avoid`.
    std::vector<int32_t> dist(n, -1);
    std::vector<NodeID> queue;
    queue.reserve(n);
    for (NodeID c : *avoid) {
      if (c < n && dist[c] == -1) {
        dist[c] = 0;
        queue.push_back(c);
      }
    }
    for (size_t head = 0; head < queue.size(); ++head) {
      const NodeID u = queue[head];
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (dist[v] != -1) continue;
        dist[v] = dist[u] + 1;
        queue.push_back(v);
      }
    }
    for (NodeID v = 0; v < n; ++v)
      if (dist[v] == -1) dist[v] = std::numeric_limits<int32_t>::max();
    const size_t top = std::max<size_t>(
        num_regions, std::min<size_t>(n, static_cast<size_t>(spread_fraction * n)));
    std::nth_element(order.begin(), order.begin() + (top - 1), order.end(),
                     [&](NodeID a, NodeID b) { return dist[a] > dist[b]; });
    std::shuffle(order.begin(), order.begin() + top, rng);
  }
  std::vector<NodeID> centers(order.begin(), order.begin() + num_regions);

  assign_voronoi_cells(graph, centers, regions.region_of);
  for (int it = 0; it < lloyd_iterations; ++it) {
    move_centers_inward(graph, regions.region_of, centers);
    assign_voronoi_cells(graph, centers, regions.region_of);
  }
  regions.num_regions = static_cast<NodeID>(centers.size());
  regions.centers = std::move(centers);
  return regions;
}

std::vector<std::pair<NodeID, NodeID>> compute_region_edges(const FilterGraph& graph,
                                                            const Regions& regions) {
  std::vector<std::pair<NodeID, NodeID>> edges;
  for (NodeID u = 0; u < graph.numNodes(); ++u) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID ru = regions.region_of[u];
      const NodeID rv = regions.region_of[graph.adj[pos]];
      if (ru < rv) edges.emplace_back(ru, rv);
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

PairSubproblem solve_region_pair(const FilterGraph& graph, const Regions& regions, NodeID r1, NodeID r2,
                                 const RegionCutParams& params) {
  std::vector<NodeWeight> region_weight(regions.num_regions, 0);
  std::vector<NodeID> pair_vertices;
  for (NodeID v = 0; v < graph.numNodes(); ++v) {
    region_weight[regions.region_of[v]] += graph.node_weight[v];
    if (regions.region_of[v] == r1) pair_vertices.push_back(v);
  }
  for (NodeID v = 0; v < graph.numNodes(); ++v)
    if (regions.region_of[v] == r2) pair_vertices.push_back(v);

  PairSubproblem result;
  PairScratch scratch;
  result.flow = cut_region_pair(graph, regions.region_of, region_weight, r1, r2, pair_vertices, params, scratch,
                                [&](EdgeID e) { result.cut_edges.push_back(e); }, [](NodeID) {});
  std::sort(result.cut_edges.begin(), result.cut_edges.end());
  result.cut_edges.erase(std::unique(result.cut_edges.begin(), result.cut_edges.end()), result.cut_edges.end());

  const bool sink_side_cut = params.cut_side == CutSide::Sink;
  for (NodeID v : scratch.subgraph) {
    const NodeID id = scratch.local_id[v];
    const bool in_pair = regions.region_of[v] == r1 || regions.region_of[v] == r2;
    PairRole role = in_pair ? PairRole::Free : PairRole::Margin;
    if (id == 0) role = in_pair ? PairRole::Source : PairRole::RingSource;
    if (id == 1) role = in_pair ? PairRole::Sink : PairRole::RingSink;
    result.vertices.push_back(v);
    result.role.push_back(role);
    result.source_side.push_back(sink_side_cut ? scratch.reaches_sink[id] == 0 : scratch.reachable[id] != 0);
  }
  return result;
}

std::vector<char> run_region_cut_detection(const FilterGraph& graph, const RegionCutParams& params,
                                           RegionCutStats* stats, bool verbose) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();
  const auto start = std::chrono::steady_clock::now();
  auto elapsed = [&] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  };

  std::mt19937_64 rng(params.seed);
  const bool spread = params.centers == CenterChoice::Spread;
  Regions regions = compute_regions(graph, params.num_regions, params.lloyd_iterations, rng,
                                    spread ? params.avoid_centers : nullptr, params.spread_fraction);
  std::vector<std::pair<NodeID, NodeID>> region_edges = compute_region_edges(graph, regions);
  if (verbose) {
    std::cerr << "[region-cut] " << regions.num_regions << " regions, " << region_edges.size()
              << " region edges (" << elapsed() << "s)\n";
  }

  const NodeID k = regions.num_regions;
  std::vector<size_t> region_begin(k + 1, 0);
  std::vector<NodeWeight> region_weight(k, 0);
  for (NodeID v = 0; v < n; ++v) {
    ++region_begin[regions.region_of[v] + 1];
    region_weight[regions.region_of[v]] += graph.node_weight[v];
  }
  for (NodeID r = 0; r < k; ++r) region_begin[r + 1] += region_begin[r];
  std::vector<NodeID> region_vertices(n);
  {
    std::vector<size_t> fill(region_begin.begin(), region_begin.end() - 1);
    for (NodeID v = 0; v < n; ++v) region_vertices[fill[regions.region_of[v]]++] = v;
  }

  std::vector<parallel::IntegralAtomicWrapper<uint8_t>> keep(m);
  for (size_t e = 0; e < m; ++e) keep[e] = 0;
  // Closed cuts are not symmetric: subproblem 2i encloses the first region,
  // 2i + 1 the second.
  const size_t num_subproblems = params.closed ? 2 * region_edges.size() : region_edges.size();
  std::vector<int64_t> flow_value(num_subproblems, 0);
  std::vector<parallel::IntegralAtomicWrapper<uint8_t>> covered(params.closed ? n : 0);
  for (auto& c : covered) c = 0;

  tls_enumerable_thread_specific<PairScratch> scratch;
  tls_enumerable_thread_specific<std::vector<NodeID>> pair_vertices;
  tbb::parallel_for(size_t(0), num_subproblems, [&](size_t i) {
    auto [a, b] = region_edges[params.closed ? i / 2 : i];
    if (params.closed && i % 2 == 1) std::swap(a, b);
    std::vector<NodeID>& vertices = pair_vertices.local();
    vertices.assign(region_vertices.begin() + region_begin[a],
                    region_vertices.begin() + region_begin[a + 1]);
    vertices.insert(vertices.end(), region_vertices.begin() + region_begin[b],
                    region_vertices.begin() + region_begin[b + 1]);
    flow_value[i] = cut_region_pair(graph, regions.region_of, region_weight, a, b, vertices,
                                    params, scratch.local(),
                                    [&](EdgeID e) { keep[e].store(1, std::memory_order_relaxed); },
                                    [&](NodeID v) { covered[v].store(1, std::memory_order_relaxed); });
  });
  if (verbose) std::cerr << "[region-cut] " << num_subproblems << " cuts done (" << elapsed() << "s)\n";

  std::vector<char> result(m);
  for (size_t e = 0; e < m; ++e) result[e] = keep[e].load(std::memory_order_relaxed);
  if (stats) {
    stats->regions = std::move(regions);
    stats->region_edges = std::move(region_edges);
    stats->flow_value = std::move(flow_value);
    stats->covered.assign(n, 0);
    for (size_t v = 0; v < covered.size(); ++v) stats->covered[v] = covered[v].load(std::memory_order_relaxed);
  }
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
