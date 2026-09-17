// tests/partition/preprocessing/filtering/filter_graph_test.cc
#include <gtest/gtest.h>
#include <algorithm>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

using namespace mt_kahypar::filtering;

TEST(FilterGraphTest, BuildsTriangle) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}, {0, 2, 3}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ASSERT_EQ(graph.numNodes(), 3u);
  ASSERT_EQ(graph.numEdges(), 3u);
  EXPECT_EQ(graph.degree(0), 2u);
  EXPECT_EQ(graph.degree(1), 2u);
  EXPECT_EQ(graph.degree(2), 2u);

  std::vector<std::pair<NodeID, EdgeWeight>> incident_to_0;
  for (EdgeID pos = graph.node_begin[0]; pos < graph.node_begin[1]; ++pos) {
    incident_to_0.emplace_back(graph.adj[pos], graph.edge_weight[graph.adj_edge[pos]]);
  }
  std::sort(incident_to_0.begin(), incident_to_0.end());
  EXPECT_EQ(incident_to_0[0], std::make_pair(NodeID(1), EdgeWeight(5)));
  EXPECT_EQ(incident_to_0[1], std::make_pair(NodeID(2), EdgeWeight(3)));
}

TEST(FilterGraphTest, HandlesEmptyGraph) {
  std::vector<NodeWeight> weights = {1};
  std::vector<EdgeListEntry> edges;
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  ASSERT_EQ(graph.numNodes(), 1u);
  ASSERT_EQ(graph.numEdges(), 0u);
  EXPECT_EQ(graph.degree(0), 0u);
}

TEST(FilterGraphTest, ComputesEdgeEndpoints) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  auto endpoints = compute_edge_endpoints(graph);
  ASSERT_EQ(endpoints.size(), 2u);
  auto has_pair = [&](NodeID a, NodeID b) {
    for (auto& p : endpoints) {
      if ((p.first == a && p.second == b) || (p.first == b && p.second == a)) return true;
    }
    return false;
  };
  EXPECT_TRUE(has_pair(0, 1));
  EXPECT_TRUE(has_pair(1, 2));
}
