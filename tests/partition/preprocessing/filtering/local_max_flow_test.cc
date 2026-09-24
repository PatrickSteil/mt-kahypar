#include <gtest/gtest.h>

#include <random>

#include "mt-kahypar/partition/preprocessing/filtering/local_max_flow.h"

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

namespace {
struct TestEdge { uint32_t u, v; int64_t cap; };

FlowNetwork make_network(uint32_t n, uint32_t s, uint32_t t, const std::vector<TestEdge>& edges) {
  FlowNetwork network(n);
  network.source = s;
  network.sink = t;
  for (const TestEdge& e : edges) network.add_edge(e.u, e.v, e.cap);
  return network;
}
}  // namespace

TEST(PushRelabelMaxFlowTest, ComputesKnownMaxFlowAndMinCut) {
  std::vector<TestEdge> edges = {{0, 1, 3}, {0, 2, 2}, {1, 3, 2}, {2, 3, 3}, {1, 2, 1}};
  FlowNetwork network = make_network(4, 0, 3, edges);
  EXPECT_EQ(push_relabel_max_flow(network), 5);
  std::vector<char> reachable = min_cut_reachable_from_source(network);
  EXPECT_TRUE(reachable[0]);
  EXPECT_FALSE(reachable[3]);
  int64_t crossing_capacity = 0;
  for (const TestEdge& e : edges) {
    if (reachable[e.u] != reachable[e.v]) crossing_capacity += e.cap;
  }
  EXPECT_EQ(crossing_capacity, 5);
}

TEST(PushRelabelMaxFlowTest, ExcessNodesStayOnSourceSide) {
  // s(0) -10- a(1) -1- t(2): after the preflow, a keeps excess 9 and cannot
  // reach t, so the min cut closest to s is {s, a} | {t} (capacity 1), not {s}.
  FlowNetwork network = make_network(3, 0, 2, {{0, 1, 10}, {1, 2, 1}});
  EXPECT_EQ(push_relabel_max_flow(network), 1);
  std::vector<char> reachable = min_cut_reachable_from_source(network);
  EXPECT_TRUE(reachable[0]);
  EXPECT_TRUE(reachable[1]);
  EXPECT_FALSE(reachable[2]);
}

TEST(PushRelabelMaxFlowTest, MatchesDinicOnRandomGraphs) {
  std::mt19937 rng(7);
  for (int round = 0; round < 200; ++round) {
    const uint32_t n = 2 + rng() % 40;
    std::vector<TestEdge> edges;
    const uint32_t m = rng() % (3 * n);
    for (uint32_t i = 0; i < m; ++i) {
      const uint32_t u = rng() % n, v = rng() % n;
      if (u != v) edges.push_back({u, v, static_cast<int64_t>(1 + rng() % 20)});
    }
    FlowNetwork dinic = make_network(n, 0, 1, edges);
    FlowNetwork push_relabel = make_network(n, 0, 1, edges);
    const int64_t flow = dinic_max_flow(dinic);
    ASSERT_EQ(push_relabel_max_flow(push_relabel), flow) << "round " << round;
    // Both must yield the same source-side minimal cut
    EXPECT_EQ(min_cut_reachable_from_source(dinic), min_cut_reachable_from_source(push_relabel))
        << "round " << round;
  }
}

TEST(PushRelabelMaxFlowTest, NetworkCanBeReusedAcrossSolves) {
  FlowNetwork network;
  for (int round = 0; round < 3; ++round) {
    network.reset(3);
    network.source = 0;
    network.sink = 2;
    network.add_edge(0, 1, 4 + round);
    network.add_edge(1, 2, 7);
    EXPECT_EQ(push_relabel_max_flow(network), 4 + round);
  }
}

TEST(MinCutSidesTest, SourceAndSinkSideCutsDifferOnPathWithTiedCuts) {
  // Path s(0)-a(1)-t(2) with unit capacities: both {s-a} and {a-t} are
  // minimum cuts. The source-side cut is closest to s, the sink-side cut
  // closest to t, so a lies on different sides.
  for (const bool use_push_relabel : {false, true}) {
    FlowNetwork network(3);
    network.source = 0;
    network.sink = 2;
    network.add_edge(0, 1, 1);
    network.add_edge(1, 2, 1);
    EXPECT_EQ(use_push_relabel ? push_relabel_max_flow(network) : dinic_max_flow(network), 1);

    std::vector<char> reachable;
    std::vector<char> reaches_sink;
    std::vector<uint32_t> queue;
    min_cut_reachable_from_source(network, reachable, queue);
    min_cut_reaching_sink(network, reaches_sink, queue);
    EXPECT_TRUE(reachable[0]);
    EXPECT_FALSE(reachable[1]);
    EXPECT_FALSE(reachable[2]);
    EXPECT_FALSE(reaches_sink[0]);
    EXPECT_FALSE(reaches_sink[1]);
    EXPECT_TRUE(reaches_sink[2]);
  }
}
