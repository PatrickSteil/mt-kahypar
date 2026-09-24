#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

struct ContractionResult {
  FilterGraph graph;
  std::vector<NodeID> mapping;  // input node id -> output node id
};

// Contracts `graph` according to the equivalence classes of `uf` (a
// union-find over graph.numNodes() ids). Vertex weights of a class are
// summed; parallel edges created by contraction are merged by summing their
// weights; self-loops created by contraction are dropped.
ContractionResult contract_graph(const FilterGraph& graph, const AtomicUnionFind& uf);

// Contracts `graph` according to `mapping` (input node id -> output node id in
// [0, num_nodes)): output node i gets the summed weight of all input nodes
// mapped to it; parallel edges are merged by summing their weights and
// self-loops are dropped.
FilterGraph contract_graph_by_mapping(const FilterGraph& graph,
                                      const std::vector<NodeID>& mapping,
                                      size_t num_nodes);

}  // namespace filtering
}  // namespace mt_kahypar
