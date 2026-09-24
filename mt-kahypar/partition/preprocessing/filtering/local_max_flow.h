#pragma once

#include <cstdint>
#include <vector>

namespace mt_kahypar {
namespace filtering {

// A small local flow network used for one natural-cut subproblem. Nodes are
// local ids [0, numNodes()); `source` and `sink` designate s and t.
//
// Stored in CSR form so that one instance can be reused across many
// subproblems without per-node heap allocations: reset() keeps all buffer
// capacity, add_edge() only appends to a flat edge list, and finalize()
// (called implicitly by the solvers) builds the arc arrays.
class FlowNetwork {
 public:
  uint32_t source = 0;
  uint32_t sink = 0;

  FlowNetwork() = default;
  explicit FlowNetwork(uint32_t num_nodes) { reset(num_nodes); }

  // Clears the network and resizes it to `num_nodes` nodes, keeping capacity.
  void reset(uint32_t num_nodes) {
    num_nodes_ = num_nodes;
    edges_.clear();
    finalized_ = false;
  }

  uint32_t numNodes() const { return num_nodes_; }

  // Adds one undirected edge of capacity `cap`, modeled as a pair of
  // opposite arcs that each have capacity `cap` (standard technique).
  void add_edge(uint32_t u, uint32_t v, int64_t cap) { edges_.push_back(Edge{u, v, cap}); }

  // Builds the CSR arc arrays from the added edges. Idempotent.
  void finalize();

  // CSR access (valid after finalize()). Arcs of node u are
  // [firstArc(u), firstArc(u + 1)); residual(a) is the remaining capacity of
  // arc a, and reverse(a) is the index of its opposite arc.
  uint32_t firstArc(uint32_t u) const { return first_arc_[u]; }
  uint32_t head(uint32_t a) const { return head_[a]; }
  int64_t residual(uint32_t a) const { return residual_[a]; }
  uint32_t reverse(uint32_t a) const { return reverse_[a]; }
  // Excess left at u by push_relabel_max_flow (always 0 after dinic_max_flow).
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
  std::vector<uint32_t> first_arc_;  // size numNodes() + 1
  std::vector<uint32_t> head_;
  std::vector<int64_t> residual_;
  std::vector<uint32_t> reverse_;

  // Solver scratch, kept here so repeated solves do not reallocate.
  std::vector<int32_t> level_;
  std::vector<uint32_t> current_arc_;
  std::vector<uint32_t> queue_;
  std::vector<uint32_t> path_;
  std::vector<int64_t> excess_;
};

// Both solvers yield the same (source-side minimal) min cut; push-relabel is
// several times faster on natural-cut subproblems with non-uniform
// capacities, where Dinic needs many phases with few augmentations each.

// Runs Dinic's max-flow algorithm and returns the flow value. Mutates
// `network`'s arc capacities in place (to residual capacities).
int64_t dinic_max_flow(FlowNetwork& network);

// Computes a maximum preflow with FIFO push-relabel and global relabeling
// (only the first phase: the min cut is fully determined by a max preflow,
// so excess is never returned to the source). Returns the flow value.
// Mutates arc capacities to residual capacities.
int64_t push_relabel_max_flow(FlowNetwork& network);

// Must be called after dinic_max_flow or push_relabel_max_flow on the same
// network. Writes, per node, whether it is on the source side of the min cut
// closest to the source: reachable in the final residual graph from `source`
// (or, after a preflow, from a node holding excess).
void min_cut_reachable_from_source(const FlowNetwork& network, std::vector<char>& reachable,
                                   std::vector<uint32_t>& queue);

// Must be called after dinic_max_flow or push_relabel_max_flow on the same
// network. Writes, per node, whether it can reach `sink` in the final residual
// graph. The complement is the source side of the min cut closest to the sink.
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
