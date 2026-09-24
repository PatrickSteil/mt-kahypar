#include "mt-kahypar/partition/preprocessing/filtering/local_max_flow.h"

#include <algorithm>
#include <limits>

namespace mt_kahypar {
namespace filtering {

void FlowNetwork::finalize() {
  if (finalized_) return;
  first_arc_.assign(num_nodes_ + 1, 0);
  for (const Edge& e : edges_) {
    ++first_arc_[e.u + 1];
    ++first_arc_[e.v + 1];
  }
  for (uint32_t u = 0; u < num_nodes_; ++u) first_arc_[u + 1] += first_arc_[u];

  const size_t num_arcs = 2 * edges_.size();
  head_.resize(num_arcs);
  residual_.resize(num_arcs);
  reverse_.resize(num_arcs);
  // current_arc_ doubles as the fill cursor here
  current_arc_.assign(first_arc_.begin(), first_arc_.end() - 1);
  for (const Edge& e : edges_) {
    const uint32_t a = current_arc_[e.u]++;
    const uint32_t b = current_arc_[e.v]++;
    head_[a] = e.v;
    head_[b] = e.u;
    residual_[a] = e.cap;
    residual_[b] = e.cap;
    reverse_[a] = b;
    reverse_[b] = a;
  }
  finalized_ = true;
}

namespace {
// Exact distances to the sink in the residual graph (a backward BFS), so the
// blocking-flow search only ever walks towards the sink. Nodes that cannot
// reach the sink keep level -1 and are never entered. Stops once the source
// is labeled: nodes farther away cannot lie on a shortest augmenting path.
bool compute_levels_to_sink(const FlowNetwork& net, std::vector<int32_t>& level,
                            std::vector<uint32_t>& queue) {
  level.assign(net.numNodes(), -1);
  queue.clear();
  level[net.sink] = 0;
  queue.push_back(net.sink);
  for (size_t head = 0; head < queue.size(); ++head) {
    const uint32_t v = queue[head];
    if (level[v] >= level[net.source] && level[net.source] >= 0) break;
    for (uint32_t a = net.firstArc(v); a < net.firstArc(v + 1); ++a) {
      const uint32_t w = net.head(a);
      // w can push into v iff the arc w -> v (the reverse of a) has residual capacity
      if (level[w] == -1 && net.residual(net.reverse(a)) > 0) {
        level[w] = level[v] + 1;
        queue.push_back(w);
      }
    }
  }
  return level[net.source] != -1;
}
}  // namespace

int64_t dinic_max_flow(FlowNetwork& network) {
  network.finalize();
  const uint32_t s = network.source;
  const uint32_t t = network.sink;
  std::vector<int32_t>& level = network.level_;
  std::vector<uint32_t>& current_arc = network.current_arc_;
  std::vector<uint32_t>& path = network.path_;
  std::vector<int64_t>& residual = network.residual_;
  int64_t total_flow = 0;
  network.excess_.clear();  // a flow, not a preflow: no excess nodes
  if (s == t) return 0;

  while (compute_levels_to_sink(network, level, network.queue_)) {
    current_arc.assign(network.first_arc_.begin(), network.first_arc_.end() - 1);
    path.clear();
    uint32_t u = s;
    while (true) {
      if (u == t) {
        int64_t bottleneck = std::numeric_limits<int64_t>::max();
        for (uint32_t a : path) bottleneck = std::min(bottleneck, residual[a]);
        for (uint32_t a : path) {
          residual[a] -= bottleneck;
          residual[network.reverse_[a]] += bottleneck;
        }
        total_flow += bottleneck;
        // Retreat to the tail of the first saturated arc on the path
        size_t k = 0;
        while (residual[path[k]] > 0) ++k;
        path.resize(k);
        u = k == 0 ? s : network.head_[path[k - 1]];
        continue;
      }

      const uint32_t end = network.first_arc_[u + 1];
      uint32_t& a = current_arc[u];
      while (a < end && (residual[a] == 0 || level[network.head_[a]] != level[u] - 1)) ++a;
      if (a < end) {
        path.push_back(a);
        u = network.head_[a];
      } else {
        // Dead end: no augmenting path through u in this phase
        level[u] = -1;
        if (u == s) break;
        const uint32_t back = path.back();
        path.pop_back();
        u = network.head_[network.reverse_[back]];
        ++current_arc[u];
      }
    }
  }
  return total_flow;
}

int64_t push_relabel_max_flow(FlowNetwork& network) {
  network.finalize();
  const uint32_t n = network.numNodes();
  const uint32_t s = network.source;
  const uint32_t t = network.sink;
  if (s == t) return 0;
  std::vector<int32_t>& label = network.level_;
  std::vector<uint32_t>& current_arc = network.current_arc_;
  std::vector<int64_t>& excess = network.excess_;
  std::vector<int64_t>& residual = network.residual_;
  const std::vector<uint32_t>& first_arc = network.first_arc_;
  const std::vector<uint32_t>& head = network.head_;
  const std::vector<uint32_t>& reverse = network.reverse_;
  std::vector<uint32_t>& active = network.path_;  // FIFO queue of active nodes
  const auto max_label = static_cast<int32_t>(n);

  excess.assign(n, 0);
  current_arc.assign(first_arc.begin(), first_arc.end() - 1);
  active.clear();

  // Global relabeling: exact residual distances to t. Nodes that cannot
  // reach t get max_label and are never touched again.
  auto global_relabel = [&] {
    std::vector<uint32_t>& queue = network.queue_;
    label.assign(n, max_label);
    queue.clear();
    label[t] = 0;
    queue.push_back(t);
    for (size_t i = 0; i < queue.size(); ++i) {
      const uint32_t v = queue[i];
      for (uint32_t a = first_arc[v]; a < first_arc[v + 1]; ++a) {
        const uint32_t w = head[a];
        if (label[w] == max_label && w != s && residual[reverse[a]] > 0) {
          label[w] = label[v] + 1;
          queue.push_back(w);
        }
      }
    }
    for (uint32_t u = 0; u < n; ++u) current_arc[u] = first_arc[u];
  };

  global_relabel();
  label[s] = max_label;
  for (uint32_t a = first_arc[s]; a < first_arc[s + 1]; ++a) {
    const uint32_t v = head[a];
    const int64_t d = residual[a];
    if (d == 0 || v == s) continue;
    residual[a] = 0;
    residual[reverse[a]] += d;
    if (excess[v] == 0 && v != t) active.push_back(v);
    excess[v] += d;
  }

  const size_t num_arcs = head.size();
  const size_t global_relabel_threshold = 6 * static_cast<size_t>(n) + num_arcs;
  size_t work_since_relabel = 0;
  size_t queue_head = 0;
  while (queue_head < active.size()) {
    if (work_since_relabel > global_relabel_threshold) {
      global_relabel();
      work_since_relabel = 0;
    }
    const uint32_t u = active[queue_head++];
    // Discharge u until its excess is gone or it provably cannot reach t
    while (excess[u] > 0 && label[u] < max_label) {
      const uint32_t end = first_arc[u + 1];
      uint32_t& a = current_arc[u];
      for (; a < end; ++a) {
        const uint32_t v = head[a];
        if (residual[a] > 0 && label[u] == label[v] + 1) {
          const int64_t d = std::min(excess[u], residual[a]);
          residual[a] -= d;
          residual[reverse[a]] += d;
          excess[u] -= d;
          if (excess[v] == 0 && v != t && v != s) active.push_back(v);
          excess[v] += d;
          if (excess[u] == 0) break;
        }
      }
      if (excess[u] == 0) break;
      // Relabel: every admissible arc is used up
      int32_t new_label = max_label;
      for (uint32_t b = first_arc[u]; b < end; ++b) {
        if (residual[b] > 0) new_label = std::min(new_label, label[head[b]] + 1);
      }
      work_since_relabel += end - first_arc[u] + 12;
      label[u] = new_label;
      a = first_arc[u];
    }
    // Compact the FIFO queue once its consumed prefix dominates
    if (queue_head > 1024 && 2 * queue_head > active.size()) {
      active.erase(active.begin(), active.begin() + queue_head);
      queue_head = 0;
    }
  }
  return excess[t];
}

void min_cut_reachable_from_source(const FlowNetwork& network, std::vector<char>& reachable,
                                   std::vector<uint32_t>& queue) {
  reachable.assign(network.numNodes(), 0);
  queue.clear();
  reachable[network.source] = 1;
  queue.push_back(network.source);
  // After a max preflow, nodes still holding excess are cut off from the
  // sink and belong to the source side as well.
  for (uint32_t u = 0; u < network.numNodes(); ++u) {
    if (network.excess(u) > 0 && u != network.sink && !reachable[u]) {
      reachable[u] = 1;
      queue.push_back(u);
    }
  }
  for (size_t head = 0; head < queue.size(); ++head) {
    const uint32_t u = queue[head];
    for (uint32_t a = network.firstArc(u); a < network.firstArc(u + 1); ++a) {
      const uint32_t v = network.head(a);
      if (network.residual(a) > 0 && !reachable[v]) {
        reachable[v] = 1;
        queue.push_back(v);
      }
    }
  }
}

void min_cut_reaching_sink(const FlowNetwork& network, std::vector<char>& reaches_sink,
                           std::vector<uint32_t>& queue) {
  reaches_sink.assign(network.numNodes(), 0);
  queue.clear();
  reaches_sink[network.sink] = 1;
  queue.push_back(network.sink);
  // Backward search: u reaches v if the residual capacity of arc u -> v,
  // i.e., of the reverse of arc v -> u, is positive.
  for (size_t head = 0; head < queue.size(); ++head) {
    const uint32_t v = queue[head];
    for (uint32_t a = network.firstArc(v); a < network.firstArc(v + 1); ++a) {
      const uint32_t u = network.head(a);
      if (network.residual(network.reverse(a)) > 0 && !reaches_sink[u]) {
        reaches_sink[u] = 1;
        queue.push_back(u);
      }
    }
  }
}

}  // namespace filtering
}  // namespace mt_kahypar
