#pragma once

#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

struct FilteringParams {
  NodeWeight U;
  NodeWeight tau = 5;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  CutSide cut_side = CutSide::Source;
  bool verbose = false;
};

// PUNCH filtering (Delling et al., "Graph Partitioning with Natural Cuts").
// With alpha <= 1, no fragment is heavier than U.
FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params);

}  // namespace filtering
}  // namespace mt_kahypar
