#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

using namespace mt_kahypar::filtering;

namespace {
FilterGraph make_two_triangles() {
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(ParallelConnectivityTest, FindsTwoComponents) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  EXPECT_EQ(component[0], component[1]);
  EXPECT_EQ(component[1], component[2]);
  EXPECT_EQ(component[3], component[4]);
  EXPECT_EQ(component[4], component[5]);
  EXPECT_NE(component[0], component[3]);
}

TEST(ParallelConnectivityTest, ExcludingAnEdgeCanSplitAComponent) {
  // Excluding two adjacent edges of a 4-cycle isolates vertex 1.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::vector<char> excluded(4, 0);
  excluded[0] = 1;  // edge (0,1)
  excluded[1] = 1;  // edge (1,2)
  std::vector<NodeID> component = parallel_connected_components_excluding(graph, excluded);
  EXPECT_NE(component[1], component[3]);
}

TEST(ParallelConnectivityTest, SpanningForestCoversAllVerticesWithDefaultRoots) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});

  ASSERT_EQ(forest.roots.size(), 2u);
  ASSERT_EQ(forest.bfs_order.size(), 6u);
  for (NodeID v = 0; v < 6; ++v) EXPECT_NE(forest.parent[v], kInvalidNode);
  for (NodeID root : forest.roots) {
    EXPECT_EQ(forest.parent[root], root);
    EXPECT_EQ(forest.parent_edge[root], kInvalidEdge);
  }
  for (NodeID v = 0; v < 6; ++v) {
    if (forest.parent[v] != v) EXPECT_NE(forest.parent_edge[v], kInvalidEdge);
  }
}

TEST(ParallelConnectivityTest, SpanningForestHonorsExplicitRoots) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  std::vector<NodeID> roots = {2, 5};
  SpanningForest forest = build_spanning_forest(graph, component, roots);
  EXPECT_EQ(forest.parent[2], 2u);
  EXPECT_EQ(forest.parent[5], 5u);
  ASSERT_EQ(forest.bfs_order.size(), 6u);
  for (NodeID v = 0; v < 6; ++v) {
    EXPECT_NE(forest.parent[v], kInvalidNode);
    const NodeID owning_root = (component[v] == component[2]) ? NodeID(2) : NodeID(5);
    NodeID cur = v;
    while (forest.parent[cur] != cur) cur = forest.parent[cur];
    EXPECT_EQ(cur, owning_root);
    if (v != 2 && v != 5) EXPECT_NE(forest.parent_edge[v], kInvalidEdge);
  }
}
