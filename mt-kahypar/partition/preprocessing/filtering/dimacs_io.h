#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Reads a 9th DIMACS Implementation Challenge ".gr" graph file. Format: a
// header line "p sp <n> <m>" (n nodes, m directed arcs -- each undirected
// edge appears once per direction), followed by "m" arc lines "a <u> <v> <w>"
// (1-indexed). Comment lines start with 'c' and are ignored. Vertex weights
// default to 1. Self-loops in the file are dropped. Throws std::runtime_error
// if the file cannot be opened, has no header, references an out-of-range
// node id, or if the two directions of an undirected edge disagree on weight.
FilterGraph read_dimacs_graph(const std::string& path);

}  // namespace filtering
}  // namespace mt_kahypar
