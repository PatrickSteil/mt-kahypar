#pragma once

#include <cstdint>
#include <vector>

namespace mt_kahypar {
namespace filtering {

// A small local flow network used for one natural-cut subproblem. Nodes are
// local ids [0, numNodes()); `source` and `sink` designate s and t.
struct FlowNetwork {
  struct Arc {
    uint32_t to;
    int64_t capacity;
    uint32_t reverse_index;  // index of the reverse arc in adj[to]
  };

  std::vector<std::vector<Arc>> adj;
  uint32_t source = 0;
  uint32_t sink = 0;

  explicit FlowNetwork(uint32_t num_nodes) : adj(num_nodes) {}

  uint32_t numNodes() const { return static_cast<uint32_t>(adj.size()); }

  // Adds a pair of forward/backward arcs with equal capacity, correctly
  // modeling one undirected edge of capacity `cap` in a directed max-flow
  // formulation (standard technique).
  void add_edge(uint32_t u, uint32_t v, int64_t cap) {
    const uint32_t fwd_index = static_cast<uint32_t>(adj[u].size());
    const uint32_t bwd_index = static_cast<uint32_t>(adj[v].size());
    adj[u].push_back(Arc{v, cap, bwd_index});
    adj[v].push_back(Arc{u, cap, fwd_index});
  }
};

// Runs Dinic's max-flow algorithm and returns the flow value. Mutates
// `network`'s arc capacities in place (to residual capacities).
int64_t dinic_max_flow(FlowNetwork& network);

// Must be called after dinic_max_flow(network) on the same network. Returns,
// per node, whether it is reachable from `source` in the final residual
// graph -- the min-cut separates reachable from unreachable nodes.
std::vector<char> min_cut_reachable_from_source(const FlowNetwork& network);

}  // namespace filtering
}  // namespace mt_kahypar
