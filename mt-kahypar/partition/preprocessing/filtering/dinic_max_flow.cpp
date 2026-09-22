#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"

#include <algorithm>
#include <deque>
#include <limits>

namespace mt_kahypar {
namespace filtering {

namespace {
bool bfs_level_graph(FlowNetwork& network, std::vector<int>& level) {
  level.assign(network.numNodes(), -1);
  level[network.source] = 0;
  std::deque<uint32_t> queue{network.source};
  while (!queue.empty()) {
    const uint32_t u = queue.front();
    queue.pop_front();
    for (const auto& arc : network.adj[u]) {
      if (arc.capacity > 0 && level[arc.to] == -1) {
        level[arc.to] = level[u] + 1;
        queue.push_back(arc.to);
      }
    }
  }
  return level[network.sink] != -1;
}

int64_t send_flow(FlowNetwork& network, const std::vector<int>& level, std::vector<size_t>& it,
                   uint32_t u, int64_t pushed) {
  if (u == network.sink) return pushed;
  for (; it[u] < network.adj[u].size(); ++it[u]) {
    auto& arc = network.adj[u][it[u]];
    if (arc.capacity > 0 && level[arc.to] == level[u] + 1) {
      const int64_t bottleneck = send_flow(network, level, it, arc.to, std::min(pushed, arc.capacity));
      if (bottleneck > 0) {
        arc.capacity -= bottleneck;
        network.adj[arc.to][arc.reverse_index].capacity += bottleneck;
        return bottleneck;
      }
    }
  }
  return 0;
}
}  // namespace

int64_t dinic_max_flow(FlowNetwork& network) {
  int64_t total_flow = 0;
  std::vector<int> level;
  while (bfs_level_graph(network, level)) {
    std::vector<size_t> it(network.numNodes(), 0);
    int64_t pushed;
    while ((pushed = send_flow(network, level, it, network.source,
                                std::numeric_limits<int64_t>::max())) > 0) {
      total_flow += pushed;
    }
  }
  return total_flow;
}

std::vector<char> min_cut_reachable_from_source(const FlowNetwork& network) {
  std::vector<char> reachable(network.numNodes(), 0);
  reachable[network.source] = 1;
  std::deque<uint32_t> queue{network.source};
  while (!queue.empty()) {
    const uint32_t u = queue.front();
    queue.pop_front();
    for (const auto& arc : network.adj[u]) {
      if (arc.capacity > 0 && !reachable[arc.to]) {
        reachable[arc.to] = 1;
        queue.push_back(arc.to);
      }
    }
  }
  return reachable;
}

}  // namespace filtering
}  // namespace mt_kahypar
