#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Reads a METIS graph file using Mt-KaHyPar's own (memory-mapped, parallel)
// METIS parser. Edge weights become flow capacities and vertex weights are
// kept; both default to 1 if the file has none. Parallel edges between the
// same pair of vertices are merged by summing their weights, since
// build_csr_from_edge_list requires a simple graph. Throws on parse errors.
FilterGraph read_metis_graph(const std::string& path);

}  // namespace filtering
}  // namespace mt_kahypar
