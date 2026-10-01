#pragma once

#include <random>
#include <utility>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

// Region-based natural cuts, an alternative to PUNCH's random seeds. The graph
// is split into connected regions (graph Voronoi cells + Lloyd iterations).
// For every pair of adjacent regions (r1, r2), S1 in r1 and S2 in r2 are
// contracted into s and t and a min s-t cut is computed on r1 + r2 plus
// `margin` hops. The subgraph boundary is left open, so neighboring cuts
// overlap and close into fragment boundaries.

// Far: vertices farthest from the r1-r2 border. Center: farthest from the
// border and from the rest of the graph. Ball: a ball of ball_fraction of the
// Far set's weight around its innermost vertex.
enum class TerminalChoice { Far, Center, Ball };

// Spread: initial centers are drawn among the vertices farthest from
// avoid_centers (e.g. the centers of earlier rounds).
enum class CenterChoice { Random, Spread };

struct RegionCutParams {
  NodeID num_regions;
  int lloyd_iterations = 3;
  // Weight of S1 as a fraction of r1 (same for S2).
  double terminal_fraction = 0.3;
  int margin = 1;
  // Contract the outermost margin layer into s or t instead of leaving it open.
  bool pin_ring = false;
  // t additionally contains the layer outside the margin, so every cut fences
  // off the source side as in PUNCH. Both directions are solved. S1 is marked
  // covered if the source side can weigh at most max_enclosed_weight.
  bool closed = false;
  NodeWeight max_enclosed_weight = 0;
  TerminalChoice terminals = TerminalChoice::Far;
  double ball_fraction = 0.3;
  CenterChoice centers = CenterChoice::Random;
  double spread_fraction = 0.3;
  const std::vector<NodeID>* avoid_centers = nullptr;
  CutSide cut_side = CutSide::Source;
  uint64_t seed = 0;
};

struct Regions {
  std::vector<NodeID> region_of;
  std::vector<NodeID> centers;
  NodeID num_regions = 0;
};

// Every region is connected; components without a center get extra centers.
Regions compute_regions(const FilterGraph& graph, NodeID num_regions,
                        int lloyd_iterations, std::mt19937_64& rng,
                        const std::vector<NodeID>* avoid = nullptr, double spread_fraction = 0.3);

// Sorted pairs (r1 < r2) of adjacent regions.
std::vector<std::pair<NodeID, NodeID>> compute_region_edges(
    const FilterGraph& graph, const Regions& regions);

struct RegionCutStats {
  Regions regions;
  std::vector<std::pair<NodeID, NodeID>> region_edges;
  std::vector<int64_t> flow_value;
  std::vector<char> covered;  // closed mode only
};

enum class PairRole : uint8_t {
  Free,
  Source,      // S1
  Sink,        // S2
  Margin,
  RingSource,  // pin_ring
  RingSink     // pin_ring / closed
};

struct PairSubproblem {
  std::vector<NodeID> vertices;
  std::vector<PairRole> role;
  std::vector<char> source_side;
  std::vector<EdgeID> cut_edges;
  int64_t flow = 0;
};

// Single subproblem, for inspection and plotting.
PairSubproblem solve_region_pair(const FilterGraph& graph, const Regions& regions, NodeID r1, NodeID r2,
                                 const RegionCutParams& params);

// Returns the per-edge flags of all cut edges.
std::vector<char> run_region_cut_detection(const FilterGraph& graph,
                                           const RegionCutParams& params,
                                           RegionCutStats* stats = nullptr,
                                           bool verbose = false);

}  // namespace filtering
}  // namespace mt_kahypar
