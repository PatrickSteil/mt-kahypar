#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Parallel edges are merged by summing their weights.
FilterGraph read_metis_graph(const std::string& path);

void write_metis_graph(const FilterGraph& graph, const std::string& path);

}  // namespace filtering
}  // namespace mt_kahypar
