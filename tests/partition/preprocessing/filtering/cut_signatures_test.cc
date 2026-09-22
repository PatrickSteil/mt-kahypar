#include <algorithm>
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

using namespace mt_kahypar::filtering;

namespace {
EdgeSignatures signatures_for(const FilterGraph& graph) {
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  return compute_edge_signatures(graph, forest);
}
}  // namespace

TEST(CutSignaturesTest, EveryEdgeOfAPathIsABridge) {
  // 0-1-2-3, a path: every edge is a bridge.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  for (EdgeID e = 0; e < 3; ++e) EXPECT_TRUE(is_bridge[e]);
}

TEST(CutSignaturesTest, NoEdgeOfACycleIsABridge) {
  // A simple cycle 0-1-2-3-0: no bridges (removing any single edge leaves it
  // connected via the rest of the cycle).
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  for (EdgeID e = 0; e < 4; ++e) EXPECT_FALSE(is_bridge[e]);
}

TEST(CutSignaturesTest, BridgeConnectingTwoCyclesIsFound) {
  // Two triangles (0,1,2) and (3,4,5) joined by a single bridge edge (2,3).
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},  // the bridge
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  int bridge_count = 0;
  EdgeID bridge_id = kInvalidEdge;
  for (EdgeID e = 0; e < 7; ++e) {
    if (is_bridge[e]) { ++bridge_count; bridge_id = e; }
  }
  EXPECT_EQ(bridge_count, 1);
  auto endpoints = compute_edge_endpoints(graph);
  EXPECT_TRUE((endpoints[bridge_id] == std::make_pair(NodeID(2), NodeID(3))) ||
              (endpoints[bridge_id] == std::make_pair(NodeID(3), NodeID(2))));
}

namespace {
// Sorts each class and the outer list so equality comparisons are order-independent.
std::vector<std::vector<EdgeID>> normalize(std::vector<std::vector<EdgeID>> classes) {
  for (auto& c : classes) std::sort(c.begin(), c.end());
  std::sort(classes.begin(), classes.end());
  return classes;
}
}  // namespace

TEST(CutSignaturesTest, CycleGraphIsOneClassOfAllEdges) {
  // A 6-cycle: removing any two edges disconnects it, so all 6 edges form a
  // single equivalence class (this is the "Case B" scenario -- a class
  // containing a non-tree edge -- that a naive tree-only bucketing scheme
  // would miss; see design spec section 4.5).
  std::vector<NodeWeight> weights(6, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 5, 1}, {5, 0, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(sigs);
  ASSERT_EQ(classes.size(), 1u);
  EXPECT_EQ(classes[0].size(), 6u);
}

TEST(CutSignaturesTest, ThetaGraphHasTwoIndependentCutPairs) {
  // A "theta graph": nodes 0 and 3 connected by three internally-disjoint
  // paths: 0-1-3, 0-2-3, 0-4-3. Removing both edges of any one path
  // disconnects that path's middle vertex, but does not disconnect the whole
  // graph -- the two edges *incident to a degree-2 middle vertex* (e.g. (0,1)
  // and (1,3)) form a 2-cut isolating {1}. The three middle vertices give
  // three independent 2-cut classes of size 2 each; no larger cross-path
  // pairing is a valid cut.
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 3, 1},
    {0, 2, 1}, {2, 3, 1},
    {0, 4, 1}, {4, 3, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(sigs);
  ASSERT_EQ(classes.size(), 3u);
  for (auto& c : classes) EXPECT_EQ(c.size(), 2u);
}

TEST(CutSignaturesTest, BridgeIsExcludedButEachTriangleIsItsOwnClass) {
  // Two triangles joined by a bridge. An earlier draft of this test wrongly
  // expected NO 2-cut classes at all, reasoning that a triangle is
  // "3-edge-connected internally" -- it is not. A triangle's min edge cut
  // is 2 (any two of its three edges isolate the third vertex, e.g.
  // removing (0,1) and (1,2) here isolates vertex 1 while the rest of the
  // graph -- 0, plus the bridge, plus triangle B -- stays connected via
  // edge (0,2)). By PUNCH's own equivalence relation P ("e,f form a 2-cut,
  // but neither is a 1-cut on its own"), every pair of edges within a
  // single triangle qualifies, since individually none of a triangle's
  // edges is a bridge. So each triangle's 3 edges form their own
  // equivalence class of size 3 -- this matches the general "an n-cycle's
  // edges form one class of size n" pattern already exercised by
  // CycleGraphIsOneClassOfAllEdges above (a triangle is a 3-cycle). Only
  // the bridge (2,3) itself must never appear in any class.
  std::vector<NodeWeight> weights(6, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = normalize(find_two_edge_cut_classes(sigs));
  ASSERT_EQ(classes.size(), 2u);
  EXPECT_EQ(classes[0].size(), 3u);
  EXPECT_EQ(classes[1].size(), 3u);
  const EdgeID bridge_id = 3;  // canonical id of {2,3}, the 4th edge added (index 3)
  for (const auto& cls : classes) {
    EXPECT_EQ(std::find(cls.begin(), cls.end(), bridge_id), cls.end());
  }
}
