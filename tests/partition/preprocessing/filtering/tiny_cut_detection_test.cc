#include <gtest/gtest.h>

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
