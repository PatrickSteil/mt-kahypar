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
  // A big triangle (0,1,2) with a small pendant vertex 3 attached via a
  // bridge (2,3), and an even smaller pendant 4 attached to 3 via another
  // bridge. Weights: main triangle vertices weight 100 each (heavy = root),
  // vertices 3 and 4 weight 1 each. With U = 5, the whole {3,4} subtree
  // (total weight 2) should contract into a single vertex, still attached to
  // the triangle via the bridge (2,3) -- so the result has 4 nodes (the
  // triangle's 3 untouched vertices, plus one contracted {3,4} vertex) and
  // the bridge survives as an edge.
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
  // A 3-block chain: root(weight 1000, its subtree alone is far too big to
  // ever be chosen as a whole) -- mid(weight 2) -- leaf(weight 1). With
  // U = 10, tau = 5: mid's own subtree {mid, leaf} has weight 3 <= U, so it
  // gets chosen and contracted as a unit; since 3 <= tau AND 3 + root's own
  // weight... wait, root's own weight (1000) alone already exceeds U, so no
  // merge into root should happen here either. This deliberately isolates
  // the "chosen subtree, but tau-merge condition fails because the parent is
  // heavy" case, distinct from the plain "whole graph fits" case.
  std::vector<NodeWeight> weights = {1000, 2, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.mapping[1], result.mapping[2]);   // mid+leaf contract together
  EXPECT_NE(result.mapping[0], result.mapping[1]);   // but root stays separate
  EXPECT_EQ(result.graph.numNodes(), 2u);
}

TEST(ContractComponentTreeTest, TauMergeActuallyFusesIntoALightParent) {
  // root(weight 2, light enough itself) -- leaf(weight 2). Root's own
  // subtree (root+leaf, weight 4) already fits under U = 10 directly, so
  // this alone contracts everything without needing tau-merge -- included
  // as a baseline. The interesting case is TauMergePreventsCascadingOverflow
  // below, which is the one that actually isolates the tau-merge-into-a-
  // not-otherwise-chosen-parent path.
  std::vector<NodeWeight> weights = {2, 2};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
}

TEST(ContractComponentTreeTest, TauMergePreventsCascadingOverflow) {
  // R(weight 1000, unambiguously the heaviest block so it's picked as root)
  // -- P(weight 1, a light non-root parent) -- {child_a(weight 4),
  // child_b(weight 4)} (P's two children). With U = 6, tau = 5: P's own
  // subtree (1+4+4=9) is too big to be chosen as a whole (9 > 6), so P
  // itself is never marked chosen. child_a's subtree (weight 4) <= U and
  // <= tau, and 4 + P's own weight (1) = 5 <= U, so child_a tau-merges into
  // P. child_b's subtree also individually satisfies 4 + 1 = 5 <= U -- but
  // P has ALREADY absorbed child_a's 4 units of extra weight, so the TRUE
  // combined result of also merging child_b (1 + 4 + 4 = 9) would exceed U.
  // The implementation must track this and skip child_b's tau-merge, leaving
  // child_b's own already-contracted subtree (just itself) standing on its
  // own rather than fusing it into the now-full P group. This is the exact
  // bug caught during this plan's pre-flight review (see the Global
  // Constraints note on this task) -- if the fix regresses, this test will
  // fail by observing a merged group heavier than U.
  std::vector<NodeWeight> weights = {1000, 1, 4, 4};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {1, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{6, 5});

  // Only groups actually formed by merging >1 original vertex are bound by
  // U; an untouched singleton (like R here, deliberately weight 1000 so it
  // is picked as root and legitimately stays separate -- see
  // EXPECT_NE(result.mapping[0], result.mapping[1]) below, and the same
  // pattern in TauMergeFoldsSmallSubtreeIntoParentSpecifically and
  // LargeGraphIsLeftUntouchedWhenNoSubtreeFits) can exceed U on its own,
  // since Part 1 never touches a block it doesn't choose to contract.
  std::vector<size_t> group_size(result.graph.numNodes(), 0);
  for (size_t v = 0; v < graph.numNodes(); ++v) ++group_size[result.mapping[v]];
  for (size_t v = 0; v < result.graph.numNodes(); ++v) {
    if (group_size[v] > 1) {
      EXPECT_LE(result.graph.node_weight[v], 6u) << "merged group " << v << " exceeds U";
    }
  }
  EXPECT_NE(result.mapping[0], result.mapping[1]);  // R never merges with anything
  // Exactly one of P's two children fused with P; the other stands alone.
  const bool a_fused_with_p = (result.mapping[1] == result.mapping[2]);
  const bool b_fused_with_p = (result.mapping[1] == result.mapping[3]);
  EXPECT_TRUE(a_fused_with_p != b_fused_with_p);
}

TEST(ContractComponentTreeTest, LargeGraphIsLeftUntouchedWhenNoSubtreeFits) {
  // A single triangle where every vertex is heavy: no subtree (other than
  // the whole graph, which isn't a proper subtree hanging off a bridge --
  // there are no bridges at all here) can be <= U, so nothing contracts.
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
// A "double bridge": hub 0 with a pendant leaf 1 (bridge (0,1)), hub 2 with a
// pendant leaf 3 (bridge (2,3)), and hub 0 connected to hub 2 by TWO parallel
// edges. The two edges between hub 0 and hub 2 are literal parallel edges
// between the same vertex pair, so their fundamental-cycle relationship is
// direct (one is the spanning tree edge, the other's fundamental cycle
// covers exactly that one edge and nothing else) -- no ambient cycle pulls
// in the bridges (0,1)/(2,3), which are correctly excluded as 1-cuts. The
// genuine 2-cut class is exactly the 2 parallel edges, and removing them
// isolates {0,1} from {2,3}, both genuine multi-vertex components.
FilterGraph make_double_bridge(NodeWeight w0, NodeWeight w1, NodeWeight w2, NodeWeight w3) {
  std::vector<NodeWeight> weights = {w0, w1, w2, w3};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1},  // bridge: leaf 1 off hub 0
    {2, 3, 1},  // bridge: leaf 3 off hub 2
    {0, 2, 1},  // connector 1 of the parallel pair
    {0, 2, 1}   // connector 2 of the parallel pair
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
  // {0,1} together weigh 100, too heavy to contract at U = 5, but {2,3}
  // (weight 2) still qualifies.
  FilterGraph graph = make_double_bridge(50, 50, 1, 1);
  ContractionResult result = contract_two_edge_cuts(graph, 5);

  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[2], result.mapping[3]);
}

TEST(ContractTwoEdgeCutsTest, ManyIndependentClassesProcessInParallelCorrectly) {
  // A "star of double-bridges": a heavy central hub C (vertex 0, weight
  // 100 -- always too heavy to contract with anything) with num_units
  // satellite (hub_i, leaf_i) pairs, each connected to C via its OWN pair
  // of parallel edges (its own independent 2-edge-cut class) plus a
  // hub_i-leaf_i bridge. Deliberately NOT a chain of double-bridges: a
  // chain would make each cut isolate a GROWING prefix of the chain
  // (nesting classes' effects together), which stops being independent
  // once the prefix exceeds U -- exactly the kind of test-design mistake
  // this plan has hit before (see the design comment on make_double_bridge
  // above). A star keeps every class's "small side" to exactly one
  // {hub_i, leaf_i} pair (weight 2 <= U = 5), regardless of how many other
  // units exist, so this genuinely exercises `num_units` independent
  // classes processed in parallel by one contract_two_edge_cuts call --
  // regression-testing the fix for Task 19's real-world-validation finding
  // that this loop needed to be parallelized.
  const int num_units = 50;
  std::vector<NodeWeight> weights = {100};  // vertex 0 = C
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i < num_units; ++i) {
    const NodeID hub = static_cast<NodeID>(1 + 2 * i);
    const NodeID leaf = static_cast<NodeID>(2 + 2 * i);
    weights.push_back(1);  // hub
    weights.push_back(1);  // leaf
    edges.push_back({hub, leaf, 1});      // bridge
    edges.push_back({NodeID(0), hub, 1}); // parallel edge 1 to C
    edges.push_back({NodeID(0), hub, 1}); // parallel edge 2 to C
  }
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_two_edge_cuts(graph, 5);

  for (int i = 0; i < num_units; ++i) {
    const NodeID hub = static_cast<NodeID>(1 + 2 * i);
    const NodeID leaf = static_cast<NodeID>(2 + 2 * i);
    EXPECT_EQ(result.mapping[hub], result.mapping[leaf]) << "unit " << i;
    EXPECT_NE(result.mapping[hub], result.mapping[0]) << "unit " << i << " vs center";
  }
  // Every unit contracts to its own vertex, plus C stays on its own.
  EXPECT_EQ(result.graph.numNodes(), static_cast<size_t>(num_units + 1));
}

TEST(BoundedComponentsExcludingTest, SplitsAnIsolatedTriangleIntoThreeSingletons) {
  // A bare triangle: all 3 edges form one equivalence class (its min cut is
  // 2, not 3 -- see BridgeIsExcludedButEachTriangleIsItsOwnClass in
  // cut_signatures_test.cc), and removing all 3 must split it into exactly
  // 3 singleton components (the structural guarantee this function relies
  // on, independently verified via a 500-trial networkx brute-force check
  // before this was implemented).
  std::vector<NodeWeight> weights = {5, 7, 9};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  std::vector<EdgeID> cls = {0, 1, 2};
  std::vector<NodeID> seeds = {0, 1, 1, 2, 0, 2};

  BoundedComponentsResult result =
      compute_bounded_components_excluding(graph, cls, seeds, /*total_graph_weight=*/21, /*U=*/100);

  // Exactly 2 explicit singletons plus one inferred leftover singleton.
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
  // The leftover is small enough to qualify here (U=100), so it must be
  // materialized rather than left empty.
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
  // The primary safety net for this function: rather than trusting hand-
  // picked small cases alone (this exact class of graph-traversal logic has
  // produced subtle bugs elsewhere in this module), compare against the
  // already-correct full-scan approach (parallel_connected_components_
  // excluding + a manual weight tally) across many random graphs and random
  // small excluded-edge sets, checking both the qualifying-vertex set AND
  // the grouping agree.
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

    // Ground truth via the full-scan approach.
    std::vector<char> excl_bitmap(graph.numEdges(), 0);
    for (EdgeID e : excluded) excl_bitmap[e] = 1;
    std::vector<NodeID> comp = parallel_connected_components_excluding(graph, excl_bitmap);
    std::unordered_map<NodeID, NodeWeight> truth_weight;
    for (size_t v = 0; v < n; ++v) truth_weight[comp[v]] += weights[v];
    std::vector<char> truth_qualifies(n, 0);
    for (size_t v = 0; v < n; ++v) truth_qualifies[v] = (truth_weight[comp[v]] <= U);

    // New bounded-traversal method.
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
  // A path of 6 light vertices (0..5) hanging off a heavy triangle (6,7,8)
  // via a bridge (5,6). Pass 1 should contract the light path's bridge-
  // isolated subtree if small enough; whatever survives should still respect
  // the invariant that every original vertex maps to some final vertex, and
  // the heavy triangle should never merge with anything (it's the largest
  // component and stays intact since U is too small to touch it).
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

  // The heavy triangle vertices must all still be distinguishable from the
  // light path (they're far too heavy, individually, to be absorbed).
  EXPECT_NE(result.mapping[6], result.mapping[0]);
}
