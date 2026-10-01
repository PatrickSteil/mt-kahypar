#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Reads a DIMACS ".gr" file. Edges get unit weight unless use_arc_weights is
// set (then the minimum of both arc directions is used).
FilterGraph read_dimacs_graph(const std::string& path, bool use_arc_weights = false);

}  // namespace filtering
}  // namespace mt_kahypar
