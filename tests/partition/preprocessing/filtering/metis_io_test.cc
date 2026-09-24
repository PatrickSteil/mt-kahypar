#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>

#include "mt-kahypar/partition/preprocessing/filtering/metis_io.h"

using namespace mt_kahypar::filtering;

namespace {
void write_file(const std::string& path, const std::string& content) {
  std::ofstream out(path);
  out << content;
}

EdgeWeight weight_of(const FilterGraph& graph, NodeID a, NodeID b) {
  std::vector<std::pair<NodeID, NodeID>> endpoints = compute_edge_endpoints(graph);
  for (EdgeID e = 0; e < graph.numEdges(); ++e) {
    const auto& [u, v] = endpoints[e];
    if ((u == a && v == b) || (u == b && v == a)) return graph.edge_weight[e];
  }
  ADD_FAILURE() << "edge not found";
  return -1;
}
}  // namespace

TEST(MetisIoTest, ParsesUnweightedGraph) {
  const std::string path = "metis_io_test_unweighted.graph";
  // Triangle 1-2-3 plus pendant vertex 4 attached to 3
  write_file(path, "% comment\n4 4\n2 3\n1 3\n1 2 4\n3\n");
  FilterGraph graph = read_metis_graph(path);
  std::remove(path.c_str());

  ASSERT_EQ(graph.numNodes(), 4u);
  ASSERT_EQ(graph.numEdges(), 4u);
  EXPECT_EQ(graph.degree(2), 3u);
  EXPECT_EQ(graph.degree(3), 1u);
  for (NodeID v = 0; v < 4; ++v) EXPECT_EQ(graph.node_weight[v], 1u);
  EXPECT_EQ(weight_of(graph, 0, 1), 1);
  EXPECT_EQ(weight_of(graph, 2, 3), 1);
}

TEST(MetisIoTest, ParsesEdgeAndVertexWeights) {
  const std::string path = "metis_io_test_weighted.graph";
  // fmt 11: vertex weights, then (neighbor, edge weight) pairs
  write_file(path, "3 2 11\n5 2 7\n6 1 7 3 4\n7 2 4\n");
  FilterGraph graph = read_metis_graph(path);
  std::remove(path.c_str());

  ASSERT_EQ(graph.numNodes(), 3u);
  ASSERT_EQ(graph.numEdges(), 2u);
  EXPECT_EQ(graph.node_weight[0], 5u);
  EXPECT_EQ(graph.node_weight[1], 6u);
  EXPECT_EQ(graph.node_weight[2], 7u);
  EXPECT_EQ(weight_of(graph, 0, 1), 7);
  EXPECT_EQ(weight_of(graph, 1, 2), 4);
}

TEST(MetisIoTest, WrittenGraphReadsBackIdentically) {
  std::vector<NodeWeight> weights = {3, 1, 4, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 9}, {0, 3, 2}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  const std::string path = "metis_io_test_roundtrip.graph";
  write_metis_graph(graph, path);
  FilterGraph read_back = read_metis_graph(path);
  std::remove(path.c_str());

  ASSERT_EQ(read_back.numNodes(), 4u);
  ASSERT_EQ(read_back.numEdges(), 3u);
  for (NodeID v = 0; v < 4; ++v) EXPECT_EQ(read_back.node_weight[v], weights[v]);
  EXPECT_EQ(weight_of(read_back, 0, 1), 5);
  EXPECT_EQ(weight_of(read_back, 1, 2), 9);
  EXPECT_EQ(weight_of(read_back, 0, 3), 2);
}
