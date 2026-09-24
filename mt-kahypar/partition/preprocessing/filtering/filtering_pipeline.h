// mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h
#pragma once

#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

// -----------------------------------------------------------------------
// Extension point for a future Mt-KaHyPar coarsening integration: construct
// a FilterGraph (see filter_graph.h) from whatever hypergraph/graph
// representation the caller has, call run_filtering_pipeline, and consume
// FilteringResult.fragment_id as an alternative clustering to feed into
// coarsening -- nothing else in this module needs to be understood. No such
// integration is wired up yet (design spec section 10).
// -----------------------------------------------------------------------
struct FilteringParams {
  NodeWeight U;
  NodeWeight tau = 5;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  CutSide cut_side = CutSide::Source;
  // When set, prints per-phase timing and progress (including a running
  // count of max-flow solves in natural-cut detection) to stderr as the
  // pipeline runs. Off by default so callers (e.g. unit tests) stay silent.
  bool verbose = false;
};

// Runs the full PUNCH filtering pipeline (Parts 1-3) on `graph`. No fragment
// in the result exceeds params.U, guaranteed by alpha <= 1 (design spec
// section 6) -- enforced as a test invariant, not merely documented.
FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params);

}  // namespace filtering
}  // namespace mt_kahypar
