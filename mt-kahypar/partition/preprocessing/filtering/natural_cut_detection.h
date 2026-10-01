#pragma once

#include <random>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/local_max_flow.h"

namespace mt_kahypar {
namespace filtering {

// PUNCH filtering, part 2: natural cuts.

// Which min cut to keep: closest to the core, closest to the ring, or both.
enum class CutSide { Source, Sink, Both };

struct NaturalCutParams {
  NodeWeight U;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  CutSide cut_side = CutSide::Source;
};

// Per-thread buffers. The O(n) arrays are allocated once and only the entries
// touched by the previous call are reset.
struct NaturalCutScratch {
  std::vector<NodeID> bfs_queue;
  std::vector<char> in_tree;
  std::vector<char> in_core;
  std::vector<NodeID> tree_order;
  std::vector<char> in_ring;
  std::vector<NodeID> ring;
  std::vector<NodeID> local_id;
  FlowNetwork network;
  std::vector<char> reachable;
  std::vector<char> reaches_sink;
  std::vector<uint32_t> flow_queue;
};

// Grows a BFS tree from `seed` up to weight alpha * U, contracts the core into
// s and the ring into t, and returns the min s-t cut edges. Marks the core
// vertices as covered.
std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                        const NaturalCutParams& params,
                                        NaturalCutScratch& scratch,
                                        std::vector<char>& covered);

// `coverage` sweeps; each picks uncovered seeds in random order. Returns the
// per-edge flags of all cut edges.
std::vector<char> run_natural_cut_detection_sequential(
    const FilterGraph& graph, const NaturalCutParams& params,
    std::mt19937_64& rng);

std::vector<char> run_natural_cut_detection(const FilterGraph& graph,
                                            const NaturalCutParams& params,
                                            bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
