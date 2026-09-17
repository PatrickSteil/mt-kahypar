// tests/partition/preprocessing/filtering/parallel_connectivity_test.cc
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
  // A 4-cycle 0-1-2-3-0 is 2-edge-connected; excluding one edge still leaves
  // it connected (it's a cycle), but excluding two adjacent edges splits it.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::vector<char> excluded(4, 0);
  excluded[0] = 1;  // edge (0,1)
  excluded[1] = 1;  // edge (1,2)
  std::vector<NodeID> component = parallel_connected_components_excluding(graph, excluded);
  EXPECT_NE(component[1], component[3]);  // vertex 1 is now isolated from 3
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
  // Force roots 2 and 5 explicitly (one per component, each a valid member
  // of that component in make_two_triangles' fixed layout: {0,1,2} and
  // {3,4,5}).
  std::vector<NodeID> roots = {2, 5};
  SpanningForest forest = build_spanning_forest(graph, component, roots);
  EXPECT_EQ(forest.parent[2], 2u);
  EXPECT_EQ(forest.parent[5], 5u);
  // The explicit-roots path must be as fully covered as the default-roots
  // path: every vertex reached, every non-root's tree edge in its own
  // component, and each vertex's parent chain leading back to the root
  // that owns its component (not the other root).
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
