#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

using namespace mt_kahypar::filtering;

namespace {
FilterGraph make_path_graph() {
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 10}, {1, 2, 20}, {2, 3, 30}};
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(GraphContractionTest, MergesParallelEdgesAndSumsWeights) {
  FilterGraph graph = make_path_graph();
  AtomicUnionFind uf(4);
  uf.unite(0, 1);
  uf.unite(2, 3);

  ContractionResult result = contract_graph(graph, uf);
  ASSERT_EQ(result.graph.numNodes(), 2u);
  ASSERT_EQ(result.graph.numEdges(), 1u);
  EXPECT_EQ(result.graph.edge_weight[0], 20);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[2], result.mapping[3]);
  EXPECT_NE(result.mapping[0], result.mapping[2]);
  EXPECT_EQ(result.graph.node_weight[result.mapping[0]], 2u);
  EXPECT_EQ(result.graph.node_weight[result.mapping[2]], 2u);
}

TEST(GraphContractionTest, MergesTriangleIntoParallelEdgeThenDropsSelfLoop) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}, {0, 2, 3}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  AtomicUnionFind uf(3);
  uf.unite(1, 2);

  ContractionResult result = contract_graph(graph, uf);
  ASSERT_EQ(result.graph.numNodes(), 2u);
  ASSERT_EQ(result.graph.numEdges(), 1u);
  EXPECT_EQ(result.graph.edge_weight[0], 8);  // 5 (0-1) + 3 (0-2) merged
}

TEST(GraphContractionTest, ContractingEverythingLeavesNoEdges) {
  FilterGraph graph = make_path_graph();
  AtomicUnionFind uf(4);
  uf.unite(0, 1);
  uf.unite(1, 2);
  uf.unite(2, 3);

  ContractionResult result = contract_graph(graph, uf);
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.graph.numEdges(), 0u);
  EXPECT_EQ(result.graph.node_weight[0], 4u);
}

TEST(GraphContractionTest, ContractsByMappingKeepingOutputIds) {
  // Square 0-1-2-3-0 plus diagonal 0-2; classes {0, 3} -> 1 and {1, 2} -> 0,
  // so output ids follow the mapping, not the order of first occurrence.
  std::vector<NodeWeight> weights = {1, 2, 3, 4};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 5}, {2, 3, 1}, {3, 0, 7}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  FilterGraph contracted = contract_graph_by_mapping(graph, {1, 0, 0, 1}, 2);
  ASSERT_EQ(contracted.numNodes(), 2u);
  ASSERT_EQ(contracted.numEdges(), 1u);
  EXPECT_EQ(contracted.node_weight[0], 5u);
  EXPECT_EQ(contracted.node_weight[1], 5u);
  // Edges 0-1, 2-3 and 0-2 cross the classes; 1-2 and 3-0 become self-loops
  EXPECT_EQ(contracted.edge_weight[0], 3);
}
