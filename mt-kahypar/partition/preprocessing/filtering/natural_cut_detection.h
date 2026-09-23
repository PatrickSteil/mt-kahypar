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
// Per-thread buffers reused across compute_natural_cut calls. The O(n) arrays
// (in_tree, in_core, in_ring, local_id) are only allocated once per graph size:
// each call resets exactly the entries touched by the previous call (tree_order
// and ring), since a natural-cut subproblem only touches O(alpha * U) vertices
// while there can be hundreds of thousands of subproblems on large graphs.
struct NaturalCutScratch {
  std::vector<NodeID> bfs_queue;
  std::vector<char> in_tree;
  std::vector<char> in_core;
  std::vector<NodeID> tree_order;
  std::vector<char> in_ring;
  std::vector<NodeID> ring;
  std::vector<NodeID> local_id;
  // Used by run_natural_cut_detection as the `covered` output of a single call
  std::vector<char> covered_buffer;
};

// Grows a BFS tree from `seed` until its total vertex weight reaches
// params.U * params.alpha, computes the local min s-t cut between the
// resulting core (contracted to s) and ring (contracted to t), and returns
// the canonical edge ids of `graph` forming that min cut. Only vertices
// that end up in the CORE are marked covered[v] = true -- not every
// vertex visited by the wider BFS growth, and not ring vertices. This
// matches the paper's own stopping rule ("pick uniformly at random among
// vertices that have not yet been part of any core", design spec section
// 5) and is load-bearing for the hard U-invariant: only core membership
// guarantees a vertex's eventual fragment is bounded by U (see the
// in-line comment at the core-selection loop in the .cpp for the full
// argument). `covered` must have size graph.numNodes().
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

// Parallel counterpart to run_natural_cut_detection_sequential, with the
// same semantics (design spec section 5): per sweep, a pre-shuffled vertex
// order is scanned by a TBB parallel_for; each task atomically claims its
// vertex as a new seed via compare-exchange on a per-vertex `covered` flag
// (skipping on failure), then runs compute_natural_cut using thread-local
// scratch buffers. `keep` flags are set with plain relaxed atomic stores
// (a monotonic boolean needs no compare-exchange). A short serial mop-up
// pass handles any vertices left uncovered by races near the scan's end.
// When `verbose` is set, prints a running max-flow-solve count (one Dinic
// solve per compute_natural_cut call) every 5000 solves, plus a per-sweep
// summary, to stderr.
std::vector<char> run_natural_cut_detection(const FilterGraph& graph, const NaturalCutParams& params,
                                             bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
