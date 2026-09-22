#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"

using namespace mt_kahypar::filtering;

TEST(DinicMaxFlowTest, ComputesKnownMaxFlowAndMinCut) {
  // Undirected network: s(0)-a(1) cap 3, s(0)-b(2) cap 2, a(1)-t(3) cap 2,
  // b(2)-t(3) cap 3, a(1)-b(2) cap 1. Hand-computed max flow / min cut = 5
  // (e.g. the cut {s} alone already has capacity 3+2=5, and it's achievable:
  // push s-a=3, s-b=2, then a sends 2 to t and 1 to b, b sends 3 to t).
  struct Edge { uint32_t u, v; int64_t cap; };
  std::vector<Edge> edges = {
    {0, 1, 3}, {0, 2, 2}, {1, 3, 2}, {2, 3, 3}, {1, 2, 1}
  };
  FlowNetwork network(4);
  network.source = 0;
  network.sink = 3;
  for (const Edge& e : edges) network.add_edge(e.u, e.v, e.cap);

  const int64_t max_flow = dinic_max_flow(network);
  EXPECT_EQ(max_flow, 5);

  std::vector<char> reachable = min_cut_reachable_from_source(network);
  EXPECT_TRUE(reachable[0]);
  EXPECT_FALSE(reachable[3]);

  int64_t crossing_capacity = 0;
  for (const Edge& e : edges) {
    if (reachable[e.u] != reachable[e.v]) crossing_capacity += e.cap;
  }
  EXPECT_EQ(crossing_capacity, max_flow);
}

TEST(DinicMaxFlowTest, DisconnectedSourceAndSinkGiveZeroFlow) {
  FlowNetwork network(2);
  network.source = 0;
  network.sink = 1;
  EXPECT_EQ(dinic_max_flow(network), 0);
}
