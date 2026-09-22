// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h
#pragma once

#include <random>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"
#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct NaturalCutParams {
  NodeWeight U;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
};

// Scratch buffers reused across many BFS-growth + local-min-cut calls, to
// avoid per-call heap allocation (design spec section 5). One instance per
// thread when parallelized (Task 15).
struct NaturalCutScratch {
  std::vector<NodeID> bfs_queue;
  std::vector<char> in_tree;
  std::vector<char> in_core;
  std::vector<NodeID> tree_order;
};

// Grows a BFS tree from `seed` until its total vertex weight reaches
// params.U * params.alpha, computes the local min s-t cut between the
// resulting core (contracted to s) and ring (contracted to t), and returns
// the canonical edge ids of `graph` forming that min cut. Every vertex
// visited by the BFS growth (the tree) is marked covered[v] = true; ring
// vertices are not marked (design spec section 5). `covered` must have size
// graph.numNodes().
std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                         const NaturalCutParams& params,
                                         NaturalCutScratch& scratch,
                                         std::vector<char>& covered);

// Runs the full sequential natural-cut detection procedure (design spec
// section 5): for each of params.coverage sweeps, resets per-sweep coverage,
// visits vertices in a freshly shuffled order, and for every not-yet-covered
// vertex runs compute_natural_cut as a new seed, accumulating its cut edges
// into the returned keep set (which persists across all sweeps).
std::vector<char> run_natural_cut_detection_sequential(const FilterGraph& graph,
                                                        const NaturalCutParams& params,
                                                        std::mt19937_64& rng);

}  // namespace filtering
}  // namespace mt_kahypar
