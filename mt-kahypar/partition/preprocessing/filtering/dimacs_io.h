#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Reads a 9th DIMACS Implementation Challenge ".gr" graph file. Format: a
// header line "p sp <n> <m>" (n nodes, m directed arcs -- each undirected
// edge appears once per direction), followed by "m" arc lines "a <u> <v> <w>"
// (1-indexed). Comment lines start with 'c' and are ignored. By default, every
// edge gets weight 1, i.e., natural cuts minimize the number of cut edges, as
// in the unweighted setting of the PUNCH paper. If `use_arc_weights` is set,
// the arc weight becomes the edge's flow capacity instead; if the two
// directions of an undirected edge disagree, the minimum of the two is used.
// Vertex weights default to 1. Self-loops in the file are dropped. Throws
// std::runtime_error if the file cannot be opened, has no header, or
// references an out-of-range node id.
FilterGraph read_dimacs_graph(const std::string& path, bool use_arc_weights = false);

}  // namespace filtering
}  // namespace mt_kahypar
