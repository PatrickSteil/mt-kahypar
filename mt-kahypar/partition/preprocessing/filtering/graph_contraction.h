#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

struct ContractionResult {
  FilterGraph graph;
  std::vector<NodeID> mapping;  // old node -> new node
};

// Contracts every set of `uf` into one node. Node and parallel edge weights
// are summed, self-loops are dropped.
ContractionResult contract_graph(const FilterGraph& graph, const AtomicUnionFind& uf);

FilterGraph contract_graph_by_mapping(const FilterGraph& graph,
                                      const std::vector<NodeID>& mapping,
                                      size_t num_nodes);

}  // namespace filtering
}  // namespace mt_kahypar
