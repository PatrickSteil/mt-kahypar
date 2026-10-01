#include <algorithm>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

using namespace mt_kahypar::filtering;

TEST(ContractComponentTreeTest, ContractsASmallLeafHangingOffABridge) {
  // Heavy triangle (0,1,2) with a pendant path 2-3-4. With U = 5, {3,4}
  // contracts into one vertex.
  std::vector<NodeWeight> weights = {100, 100, 100, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},
    {3, 4, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{5, 5});
  EXPECT_EQ(result.graph.numNodes(), 4u);
  EXPECT_EQ(result.mapping[3], result.mapping[4]);
  EXPECT_NE(result.mapping[2], result.mapping[3]);
  const NodeWeight contracted_weight = result.graph.node_weight[result.mapping[3]];
  EXPECT_EQ(contracted_weight, 2u);
}

TEST(ContractComponentTreeTest, TauMergeFoldsSmallSubtreeIntoParentSpecifically) {
  // Chain root(1000) - mid(2) - leaf(1), U = 10, tau = 5: {mid, leaf}
  // contracts, but is not merged into the heavy root.
  std::vector<NodeWeight> weights = {1000, 2, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.mapping[1], result.mapping[2]);
  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.graph.numNodes(), 2u);
}

TEST(ContractComponentTreeTest, TauMergeActuallyFusesIntoALightParent) {
  // root(2) - leaf(2), U = 10: everything contracts.
  std::vector<NodeWeight> weights = {2, 2};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
}

TEST(ContractComponentTreeTest, TauMergePreventsCascadingOverflow) {
  // R(1000) - P(1) - {a(4), b(4)}, U = 6, tau = 5: a and b may each be
  // merged into P, but not both (1 + 4 + 4 > U).
  std::vector<NodeWeight> weights = {1000, 1, 4, 4};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {1, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{6, 5});

  // Only merged groups are bounded by U.
  std::vector<size_t> group_size(result.graph.numNodes(), 0);
  for (size_t v = 0; v < graph.numNodes(); ++v) ++group_size[result.mapping[v]];
  for (size_t v = 0; v < result.graph.numNodes(); ++v) {
    if (group_size[v] > 1) {
      EXPECT_LE(result.graph.node_weight[v], 6u) << "merged group " << v << " exceeds U";
    }
  }
  EXPECT_NE(result.mapping[0], result.mapping[1]);
  const bool a_fused_with_p = (result.mapping[1] == result.mapping[2]);
  const bool b_fused_with_p = (result.mapping[1] == result.mapping[3]);
  EXPECT_TRUE(a_fused_with_p != b_fused_with_p);
}

TEST(ContractComponentTreeTest, LargeGraphIsLeftUntouchedWhenNoSubtreeFits) {
  // Heavy triangle without bridges: nothing contracts.
  std::vector<NodeWeight> weights = {100, 100, 100};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.graph.numNodes(), 3u);
  EXPECT_EQ(result.graph.numEdges(), 3u);
}

TEST(ContractDegree2ChainsTest, ContractsAPathBetweenTwoAnchors) {
  // A(deg1) - B(deg2) - C(deg2) - D(deg1). B and C should merge.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 3u);
  EXPECT_EQ(result.mapping[1], result.mapping[2]);
  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_NE(result.mapping[2], result.mapping[3]);
}

TEST(ContractDegree2ChainsTest, LeavesChainUntouchedWhenTooBig) {
  std::vector<NodeWeight> weights = {1, 100, 100, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 4u);
  EXPECT_NE(result.mapping[1], result.mapping[2]);
}

TEST(ContractDegree2ChainsTest, ContractsAPureCycleWithNoAnchor) {
  // A 5-cycle: every vertex has degree 2, so there's no anchor. The whole
  // component is one chain.
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.graph.numEdges(), 0u);
}

namespace {
// Double bridge: leaf 1 - hub 0 = hub 2 - leaf 3, where hubs 0 and 2 are
// joined by two parallel edges (the 2-cut class).
FilterGraph make_double_bridge(NodeWeight w0, NodeWeight w1, NodeWeight w2, NodeWeight w3) {
  std::vector<NodeWeight> weights = {w0, w1, w2, w3};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1},
    {2, 3, 1},
    {0, 2, 1},
    {0, 2, 1}
  };
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(ContractTwoEdgeCutsTest, ContractsBothSmallSidesOfADoubleBridge) {
  FilterGraph graph = make_double_bridge(1, 1, 1, 1);
  ContractionResult result = contract_two_edge_cuts(graph, 5);

  EXPECT_EQ(result.graph.numNodes(), 2u);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[2], result.mapping[3]);
  EXPECT_NE(result.mapping[0], result.mapping[2]);
}

TEST(ContractTwoEdgeCutsTest, LeavesHeavySideUncontracted) {
  // {0,1} is too heavy for U = 5, {2,3} contracts.
  FilterGraph graph = make_double_bridge(50, 50, 1, 1);
  ContractionResult result = contract_two_edge_cuts(graph, 5);

  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[2], result.mapping[3]);
}

TEST(ContractTwoEdgeCutsTest, ManyIndependentClassesProcessInParallelCorrectly) {
  // Heavy center 0 with num_units satellites (hub, leaf), each attached to the
  // center by its own pair of parallel edges: num_units independent classes.
  const int num_units = 50;
  std::vector<NodeWeight> weights = {100};
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i < num_units; ++i) {
    const NodeID hub = static_cast<NodeID>(1 + 2 * i);
    const NodeID leaf = static_cast<NodeID>(2 + 2 * i);
    weights.push_back(1);
    weights.push_back(1);
    edges.push_back({hub, leaf, 1});
    edges.push_back({NodeID(0), hub, 1});
    edges.push_back({NodeID(0), hub, 1});
  }
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_two_edge_cuts(graph, 5);

  for (int i = 0; i < num_units; ++i) {
    const NodeID hub = static_cast<NodeID>(1 + 2 * i);
    const NodeID leaf = static_cast<NodeID>(2 + 2 * i);
    EXPECT_EQ(result.mapping[hub], result.mapping[leaf]) << "unit " << i;
    EXPECT_NE(result.mapping[hub], result.mapping[0]) << "unit " << i << " vs center";
  }
  EXPECT_EQ(result.graph.numNodes(), static_cast<size_t>(num_units + 1));
}

TEST(BoundedComponentsExcludingTest, SplitsAnIsolatedTriangleIntoThreeSingletons) {
  // Removing all edges of a triangle leaves 3 singletons.
  std::vector<NodeWeight> weights = {5, 7, 9};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  std::vector<EdgeID> cls = {0, 1, 2};
  std::vector<NodeID> seeds = {0, 1, 1, 2, 0, 2};

  BoundedComponentsResult result =
      compute_bounded_components_excluding(graph, cls, seeds, /*total_graph_weight=*/21, /*U=*/100);

  ASSERT_EQ(result.small_component_members.size(), 2u);
  NodeWeight explicit_weight = 0;
  std::unordered_set<NodeID> seen;
  for (size_t k = 0; k < result.small_component_members.size(); ++k) {
    ASSERT_EQ(result.small_component_members[k].size(), 1u);
    const NodeID v = result.small_component_members[k][0];
    EXPECT_TRUE(seen.insert(v).second) << "vertex " << v << " reported twice";
    EXPECT_EQ(result.small_component_weight[k], weights[v]);
    explicit_weight += result.small_component_weight[k];
  }
  EXPECT_EQ(result.leftover_weight, 21 - explicit_weight);
  // With U = 100, the leftover component is materialized.
  ASSERT_EQ(result.leftover_members.size(), 1u);
  EXPECT_TRUE(seen.insert(result.leftover_members[0]).second);
}

namespace {
FilterGraph make_random_connected_graph(size_t n, int extra_edges, std::mt19937_64& rng,
                                         std::vector<NodeWeight>& weights_out) {
  weights_out.resize(n);
  std::uniform_int_distribution<int> weight_dist(1, 5);
  for (auto& w : weights_out) w = weight_dist(rng);

  std::vector<NodeID> order(n);
  for (size_t i = 0; i < n; ++i) order[i] = static_cast<NodeID>(i);
  std::shuffle(order.begin(), order.end(), rng);

  std::vector<EdgeListEntry> edges;
  std::set<std::pair<NodeID, NodeID>> present;
  auto add_edge = [&](NodeID u, NodeID v) {
    if (u == v) return;
    auto key = std::minmax(u, v);
    if (!present.insert(key).second) return;
    std::uniform_int_distribution<int> ew(1, 5);
    edges.push_back({u, v, ew(rng)});
  };
  for (size_t i = 1; i < n; ++i) {
    std::uniform_int_distribution<size_t> pick(0, i - 1);
    add_edge(order[pick(rng)], order[i]);
  }
  std::uniform_int_distribution<size_t> vpick(0, n - 1);
  for (int e = 0; e < extra_edges; ++e) add_edge(static_cast<NodeID>(vpick(rng)), static_cast<NodeID>(vpick(rng)));

  return build_csr_from_edge_list(edges, weights_out);
}
}  // namespace

TEST(BoundedComponentsExcludingTest, MatchesFullScanOnRandomGraphs) {
  // Compare against a full scan on random graphs.
  std::mt19937_64 rng(42);
  const int trials = 2000;
  for (int trial = 0; trial < trials; ++trial) {
    std::uniform_int_distribution<size_t> n_dist(4, 30);
    const size_t n = n_dist(rng);
    std::vector<NodeWeight> weights;
    FilterGraph graph = make_random_connected_graph(n, 6, rng, weights);
    if (graph.numEdges() < 2) continue;

    std::uniform_int_distribution<size_t> num_excl_dist(1, std::min<size_t>(3, graph.numEdges()));
    const size_t num_excl = num_excl_dist(rng);
    std::unordered_set<EdgeID> chosen;
    std::uniform_int_distribution<size_t> edge_dist(0, graph.numEdges() - 1);
    while (chosen.size() < num_excl) chosen.insert(static_cast<EdgeID>(edge_dist(rng)));
    std::vector<EdgeID> excluded(chosen.begin(), chosen.end());

    std::vector<std::pair<NodeID, NodeID>> endpoints = compute_edge_endpoints(graph);
    std::vector<NodeID> seeds;
    for (EdgeID e : excluded) {
      seeds.push_back(endpoints[e].first);
      seeds.push_back(endpoints[e].second);
    }

    NodeWeight total_weight = 0;
    for (size_t v = 0; v < n; ++v) total_weight += weights[v];
    std::uniform_int_distribution<NodeWeight> u_dist(1, total_weight);
    const NodeWeight U = u_dist(rng);

    std::vector<char> excl_bitmap(graph.numEdges(), 0);
    for (EdgeID e : excluded) excl_bitmap[e] = 1;
    std::vector<NodeID> comp = parallel_connected_components_excluding(graph, excl_bitmap);
    std::unordered_map<NodeID, NodeWeight> truth_weight;
    for (size_t v = 0; v < n; ++v) truth_weight[comp[v]] += weights[v];
    std::vector<char> truth_qualifies(n, 0);
    for (size_t v = 0; v < n; ++v) truth_qualifies[v] = (truth_weight[comp[v]] <= U);

    BoundedComponentsResult result =
        compute_bounded_components_excluding(graph, excluded, seeds, total_weight, U);
    std::vector<char> new_qualifies(n, 0);
    std::vector<int> new_group_of(n, -1);
    int next_group = 0;
    for (size_t k = 0; k < result.small_component_members.size(); ++k) {
      if (result.small_component_weight[k] <= U) {
        for (NodeID v : result.small_component_members[k]) {
          new_qualifies[v] = 1;
          new_group_of[v] = next_group;
        }
        ++next_group;
      }
    }
    if (result.leftover_weight <= U) {
      ASSERT_TRUE(!result.leftover_members.empty() || result.leftover_weight == 0)
          << "trial " << trial << ": leftover qualifies but wasn't materialized";
      for (NodeID v : result.leftover_members) {
        new_qualifies[v] = 1;
        new_group_of[v] = next_group;
      }
      ++next_group;
    }

    for (size_t v = 0; v < n; ++v) {
      EXPECT_EQ(static_cast<bool>(truth_qualifies[v]), static_cast<bool>(new_qualifies[v]))
          << "trial " << trial << " n=" << n << " vertex " << v << " U=" << U;
    }
    for (size_t v = 0; v < n; ++v) {
      for (size_t w = v + 1; w < n; ++w) {
        if (truth_qualifies[v] && truth_qualifies[w] && comp[v] == comp[w]) {
          EXPECT_EQ(new_group_of[v], new_group_of[w])
              << "trial " << trial << " n=" << n << " v=" << v << " w=" << w << " U=" << U;
        }
      }
    }
  }
}

TEST(RunTinyCutDetectionTest, ComposesMappingAcrossAllThreePasses) {
  // Light path 0..5 attached to a heavy triangle (6,7,8) via bridge (5,6).
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1, 100, 100, 100};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 5, 1},
    {5, 6, 1},
    {6, 7, 1}, {7, 8, 1}, {6, 8, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = run_tiny_cut_detection(graph, TinyCutParams{6, 5});

  ASSERT_EQ(result.mapping.size(), 9u);
  for (NodeID v = 0; v < 9; ++v) EXPECT_LT(result.mapping[v], result.graph.numNodes());

  NodeWeight total_original = 0, total_output = 0;
  for (NodeWeight w : weights) total_original += w;
  for (size_t v = 0; v < result.graph.numNodes(); ++v) total_output += result.graph.node_weight[v];
  EXPECT_EQ(total_original, total_output);

  EXPECT_NE(result.mapping[6], result.mapping[0]);
}
