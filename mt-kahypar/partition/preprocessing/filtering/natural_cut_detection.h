// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h
#pragma once

#include <random>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/local_max_flow.h"

namespace mt_kahypar {
namespace filtering {

// Which minimum cut of a natural-cut subproblem is kept. On unit-capacity
// road networks, a subproblem usually has several distinct min cuts: Source
// keeps the one closest to the core, Sink the one closest to the ring, and
// Both keeps both (more natural cuts, hence more and smaller fragments).
enum class CutSide { Source, Sink, Both };

struct NaturalCutParams {
  NodeWeight U;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  CutSide cut_side = CutSide::Source;
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
  // Local flow network and min-cut buffers, reused across calls
  FlowNetwork network;
  std::vector<char> reachable;
  std::vector<char> reaches_sink;
  std::vector<uint32_t> flow_queue;
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
std::vector<char> run_natural_cut_detection_sequential(
    const FilterGraph& graph, const NaturalCutParams& params,
    std::mt19937_64& rng);

// Parallel counterpart to run_natural_cut_detection_sequential, with the
// same semantics (design spec section 5): per sweep, a pre-shuffled vertex
// order is scanned by a TBB parallel_for; each task atomically claims its
// vertex as a new seed via compare-exchange on a per-vertex `covered` flag
// (skipping on failure), then runs compute_natural_cut using thread-local
// scratch buffers. Core vertices are marked covered as soon as the core is
// fixed, before the flow computation, so concurrent tasks do not pick them
// as redundant seeds. `keep` flags are set with plain relaxed atomic stores
// (a monotonic boolean needs no compare-exchange). When `verbose` is set,
// prints a running max-flow-solve count (one solve per compute_natural_cut
// call) every 5000 solves, plus a per-sweep
// summary, to stderr.
std::vector<char> run_natural_cut_detection(const FilterGraph& graph,
                                            const NaturalCutParams& params,
                                            bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
