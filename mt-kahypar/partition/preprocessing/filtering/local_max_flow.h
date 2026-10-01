#pragma once

#include <cstdint>
#include <vector>

namespace mt_kahypar {
namespace filtering {

// Small flow network for one natural-cut subproblem. Reusable across
// subproblems: reset() keeps the buffers, finalize() builds the CSR arcs.
class FlowNetwork {
 public:
  uint32_t source = 0;
  uint32_t sink = 0;

  FlowNetwork() = default;
  explicit FlowNetwork(uint32_t num_nodes) { reset(num_nodes); }

  void reset(uint32_t num_nodes) {
    num_nodes_ = num_nodes;
    edges_.clear();
    finalized_ = false;
  }

  uint32_t numNodes() const { return num_nodes_; }

  // Undirected edge: two opposite arcs with capacity `cap` each.
  void add_edge(uint32_t u, uint32_t v, int64_t cap) { edges_.push_back(Edge{u, v, cap}); }

  void finalize();

  uint32_t firstArc(uint32_t u) const { return first_arc_[u]; }
  uint32_t head(uint32_t a) const { return head_[a]; }
  int64_t residual(uint32_t a) const { return residual_[a]; }
  uint32_t reverse(uint32_t a) const { return reverse_[a]; }
  int64_t excess(uint32_t u) const { return u < excess_.size() ? excess_[u] : 0; }

 private:
  friend int64_t dinic_max_flow(FlowNetwork& network);
  friend int64_t push_relabel_max_flow(FlowNetwork& network);

  struct Edge {
    uint32_t u;
    uint32_t v;
    int64_t cap;
  };

  uint32_t num_nodes_ = 0;
  bool finalized_ = false;
  std::vector<Edge> edges_;
  std::vector<uint32_t> first_arc_;
  std::vector<uint32_t> head_;
  std::vector<int64_t> residual_;
  std::vector<uint32_t> reverse_;

  std::vector<int32_t> level_;
  std::vector<uint32_t> current_arc_;
  std::vector<uint32_t> queue_;
  std::vector<uint32_t> path_;
  std::vector<int64_t> excess_;
};

int64_t dinic_max_flow(FlowNetwork& network);

// FIFO push-relabel with global relabeling, first phase only (max preflow,
// which suffices for the min cut).
int64_t push_relabel_max_flow(FlowNetwork& network);

// Source side of the min cut closest to the source.
void min_cut_reachable_from_source(const FlowNetwork& network, std::vector<char>& reachable,
                                   std::vector<uint32_t>& queue);

// Vertices that can reach the sink; the complement is the source side of the
// min cut closest to the sink.
void min_cut_reaching_sink(const FlowNetwork& network, std::vector<char>& reaches_sink,
                           std::vector<uint32_t>& queue);

inline std::vector<char> min_cut_reachable_from_source(const FlowNetwork& network) {
  std::vector<char> reachable;
  std::vector<uint32_t> queue;
  min_cut_reachable_from_source(network, reachable, queue);
  return reachable;
}

}  // namespace filtering
}  // namespace mt_kahypar
