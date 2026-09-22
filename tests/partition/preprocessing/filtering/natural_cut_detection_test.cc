#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

using namespace mt_kahypar::filtering;

TEST(ComputeNaturalCutTest, FindsTheOnlyPossibleCutOnAUnitPath) {
  // Path 0-1-2-3-4, unit weights and unit edge weights. With U = 3,
  // alpha = 1.0: BFS from seed 0 grows the tree until size >= 3, i.e.
  // tree = {0,1,2}. With f = 10: core-size threshold = 3/10 = 0.3, so core =
  // {0} only (0's own weight of 1 already exceeds 0.3). ring = neighbors of
  // the tree outside it = {3}. The local network is a unit-capacity path
  // s(0) - 1 - 2 - t(3), whose unique min cut (value 1) is the edge nearest
  // the source once Dinic's saturates it -- edge (0,1).
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  NaturalCutParams params{3, 1.0, 10.0, 2};
  NaturalCutScratch scratch;
  std::vector<char> covered(5, 0);

  std::vector<EdgeID> cut = compute_natural_cut(graph, 0, params, scratch, covered);

  ASSERT_EQ(cut.size(), 1u);
  EXPECT_EQ(cut[0], 0u);  // canonical edge id of (0,1), the first edge added

  EXPECT_TRUE(covered[0]);
  EXPECT_TRUE(covered[1]);
  EXPECT_TRUE(covered[2]);
  EXPECT_FALSE(covered[3]);  // ring, not tree -- not marked covered
  EXPECT_FALSE(covered[4]);  // never visited
}

TEST(ComputeNaturalCutTest, HandlesMultiVertexCoreRingAndBranchingCut) {
  // Edges (insertion order = canonical edge id): e0=(0,1), e1=(1,2), e2=(1,3),
  // e3=(2,4), e4=(3,4), e5=(3,5). Node weights: 0,1 = 2 each; 2,3 = 1 each;
  // 4,5 = 1 each (ring weights are never summed into tree_size/target
  // thresholds, so their exact value doesn't matter here).
  //
  // With U=6, alpha=1.0: target_tree_size=6. BFS from seed 0 (adjacency
  // built in edge-insertion order) visits 1 (via e0), then, from 1's
  // adjacency [0,2,3], visits 2 (via e1) and 3 (via e2) -- tree={0,1,2,3},
  // tree_size=2+2+1+1=6, stopping exactly at target. Vertices 4,5 are never
  // added to the tree (reached only from 2/3, which are only visited once
  // the tree already hit its target).
  //
  // With f=2: target_core_size=6/2=3. Core-selection walks tree_order
  // [0,1,2,3]: admit 0 (running=2, 2<3 continue), admit 1 (running=4, 4>=3
  // break) -- core={0,1}, weight 4.
  //
  // ring = neighbors of tree{0,1,2,3} outside it: vertex 2's only non-tree
  // neighbor is 4; vertex 3's non-tree neighbors are 4 (again -- reached
  // from TWO different tree-interior vertices, exercising the in_ring
  // dedup) and 5. ring={4,5}.
  //
  // Local network: core{0,1}->s, ring{4,5}->t, tree-interior {2,3}->local
  // ids 2,3. Edges (all added from the tree-side loop; the ring-side loop
  // is provably unreachable here since every neighbor of a ring vertex with
  // a valid local_id is necessarily a tree vertex, already handled):
  //   s-local2 (e1, cap 100), s-local3 (e2, cap 3),
  //   local2-t (e3, cap 50), local3-t (e4, cap 5), local3-t (e5, cap 5).
  // Checking all 4 subsets of {local2,local3} reachable from s gives cut
  // values {s}=103, {s,local2}=53, {s,local3}=110, {s,local2,local3}=60 --
  // a UNIQUE minimum of 53 at {s,local2}, crossed by e2 (s-local3, cap 3)
  // and e3 (local2-t, cap 50). Dinic's must find this exact cut regardless
  // of augmenting-path order, since it's the unique global minimum.
  std::vector<NodeWeight> weights = {2, 2, 1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 100}, {1, 3, 3}, {2, 4, 50}, {3, 4, 5}, {3, 5, 5}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  NaturalCutParams params{6, 1.0, 2.0, 2};
  NaturalCutScratch scratch;
  std::vector<char> covered(6, 0);

  std::vector<EdgeID> cut = compute_natural_cut(graph, 0, params, scratch, covered);

  std::sort(cut.begin(), cut.end());
  ASSERT_EQ(cut.size(), 2u);
  EXPECT_EQ(cut[0], 2u);  // e2 = (1,3)
  EXPECT_EQ(cut[1], 3u);  // e3 = (2,4)

  EXPECT_TRUE(covered[0]);
  EXPECT_TRUE(covered[1]);
  EXPECT_TRUE(covered[2]);
  EXPECT_TRUE(covered[3]);
  EXPECT_FALSE(covered[4]);  // ring, not tree
  EXPECT_FALSE(covered[5]);  // ring, not tree
}

TEST(RunNaturalCutDetectionSequentialTest, WholeSmallGraphAbsorbedYieldsNoCuts) {
  // A single triangle, U so large the whole component fits inside one BFS
  // tree: there is no "outside" to cut against, so no edges get kept.
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::mt19937_64 rng(42);
  std::vector<char> keep = run_natural_cut_detection_sequential(graph, NaturalCutParams{1000, 1.0, 10.0, 2}, rng);
  ASSERT_EQ(keep.size(), 3u);
  for (char k : keep) EXPECT_EQ(k, 0);
}

TEST(RunNaturalCutDetectionSequentialTest, LongPathKeepsSomeEdges) {
  // A path much longer than U forces multiple natural cuts along its length.
  const int len = 40;
  std::vector<NodeWeight> weights(len, 1);
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i + 1 < len; ++i) edges.push_back({static_cast<NodeID>(i), static_cast<NodeID>(i + 1), 1});
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::mt19937_64 rng(7);
  std::vector<char> keep = run_natural_cut_detection_sequential(graph, NaturalCutParams{5, 1.0, 10.0, 2}, rng);
  ASSERT_EQ(keep.size(), static_cast<size_t>(len - 1));
  bool any_kept = false;
  for (char k : keep) any_kept = any_kept || k;
  EXPECT_TRUE(any_kept);
}

TEST(RunNaturalCutDetectionParallelTest, WholeSmallGraphAbsorbedYieldsNoCuts) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::vector<char> keep = run_natural_cut_detection(graph, NaturalCutParams{1000, 1.0, 10.0, 2});
  ASSERT_EQ(keep.size(), 3u);
  for (char k : keep) EXPECT_EQ(k, 0);
}

TEST(RunNaturalCutDetectionParallelTest, LongPathKeepsSomeEdgesAndNeverCrashes) {
  const int len = 200;
  std::vector<NodeWeight> weights(len, 1);
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i + 1 < len; ++i) edges.push_back({static_cast<NodeID>(i), static_cast<NodeID>(i + 1), 1});
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  for (int trial = 0; trial < 5; ++trial) {
    std::vector<char> keep = run_natural_cut_detection(graph, NaturalCutParams{5, 1.0, 10.0, 2});
    ASSERT_EQ(keep.size(), static_cast<size_t>(len - 1));
    bool any_kept = false;
    for (char k : keep) any_kept = any_kept || k;
    EXPECT_TRUE(any_kept);
  }
}
