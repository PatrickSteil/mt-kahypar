#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

namespace mt_kahypar {
namespace filtering {

FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params) {
  ContractionResult tiny_cut_result =
      run_tiny_cut_detection(graph, TinyCutParams{params.U, params.tau});

  std::vector<char> keep = run_natural_cut_detection(
      tiny_cut_result.graph,
      NaturalCutParams{params.U, params.alpha, params.f, params.coverage});

  return assemble_fragments(tiny_cut_result.graph, keep, tiny_cut_result.mapping);
}

}  // namespace filtering
}  // namespace mt_kahypar
