#pragma once

#include <utility>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct FilteringResult {
  std::vector<NodeID> fragment_id;  // per original vertex
  std::vector<NodeWeight> fragment_size;
  std::vector<std::pair<NodeID, NodeID>> kept_edges;  // as original vertex ids
};

// PUNCH filtering, part 3: fragments are the components of `graph` without the
// kept edges. `part1_mapping` maps original vertices to vertices of `graph`.
FilteringResult assemble_fragments(const FilterGraph& graph,
                                    const std::vector<char>& keep,
                                    const std::vector<NodeID>& part1_mapping);

}  // namespace filtering
}  // namespace mt_kahypar
