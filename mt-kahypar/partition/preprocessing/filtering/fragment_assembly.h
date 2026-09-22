#pragma once

#include <utility>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct FilteringResult {
  std::vector<NodeID> fragment_id;                    // indexed by ORIGINAL vertex id
  std::vector<NodeWeight> fragment_size;               // indexed by fragment id
  std::vector<std::pair<NodeID, NodeID>> kept_edges;   // representative ORIGINAL vertex id pairs
};

// Part 3 (design spec section 6): unions the endpoints of every edge of
// `graph` NOT flagged in `keep`; the resulting connected components are the
// fragments. `graph` is Part 1's output graph, `keep` is Part 2's per-edge
// flag over it (size graph.numEdges()), and `part1_mapping` is Part 1's
// composed mapping from ORIGINAL vertex ids to `graph` vertex ids. Kept
// edges are translated back to a representative original vertex per `graph`
// vertex (its lowest original id), for later plotting against DIMACS
// coordinates.
FilteringResult assemble_fragments(const FilterGraph& graph,
                                    const std::vector<char>& keep,
                                    const std::vector<NodeID>& part1_mapping);

}  // namespace filtering
}  // namespace mt_kahypar
