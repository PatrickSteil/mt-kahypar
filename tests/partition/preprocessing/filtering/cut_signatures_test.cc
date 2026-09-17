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
