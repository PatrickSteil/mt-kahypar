#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

#include <chrono>
#include <iostream>

namespace mt_kahypar {
namespace filtering {

FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params) {
  const auto pipeline_start = std::chrono::steady_clock::now();

  if (params.verbose) std::cerr << "[filtering] Part 1: tiny-cut detection\n";
  const auto tiny_cut_start = std::chrono::steady_clock::now();
  ContractionResult tiny_cut_result = run_tiny_cut_detection(
      graph, TinyCutParams{params.U, params.tau}, params.verbose);
  if (params.verbose) {
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tiny_cut_start).count();
    std::cerr << "[filtering] Part 1 done (" << seconds << "s), " << tiny_cut_result.graph.numNodes()
              << " vertices remain\n[filtering] Part 2: natural-cut detection\n";
  }

  const auto natural_cut_start = std::chrono::steady_clock::now();
  std::vector<char> keep = run_natural_cut_detection(
      tiny_cut_result.graph,
      NaturalCutParams{params.U, params.alpha, params.f, params.coverage},
      params.verbose);
  if (params.verbose) {
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - natural_cut_start).count();
    std::cerr << "[filtering] Part 2 done (" << seconds << "s)\n[filtering] Part 3: fragment assembly\n";
  }

  const auto assembly_start = std::chrono::steady_clock::now();
  FilteringResult result = assemble_fragments(tiny_cut_result.graph, keep, tiny_cut_result.mapping);
  if (params.verbose) {
    const double assembly_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - assembly_start).count();
    const double total_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
    std::cerr << "[filtering] Part 3 done (" << assembly_seconds << "s), " << result.fragment_size.size()
              << " fragments\n[filtering] total: " << total_seconds << "s\n";
  }

  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
