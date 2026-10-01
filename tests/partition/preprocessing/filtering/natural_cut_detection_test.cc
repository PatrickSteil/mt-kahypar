#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include <tbb/global_control.h>

#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

using namespace mt_kahypar::filtering;

TEST(ComputeNaturalCutTest, FindsTheOnlyPossibleCutOnAUnitPath) {
  // Path 0-1-2-3-4, U = 3: tree {0,1,2}, core {0}, ring {3}. Min cut (0,1).
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  NaturalCutParams params{3, 1.0, 10.0, 2};
  NaturalCutScratch scratch;
  std::vector<char> covered(5, 0);

  std::vector<EdgeID> cut = compute_natural_cut(graph, 0, params, scratch, covered);

  ASSERT_EQ(cut.size(), 1u);
  EXPECT_EQ(cut[0], 0u);

  // Only core vertices are marked covered.
  EXPECT_TRUE(covered[0]);
  EXPECT_FALSE(covered[1]);
  EXPECT_FALSE(covered[2]);
  EXPECT_FALSE(covered[3]);
  EXPECT_FALSE(covered[4]);
}

TEST(ComputeNaturalCutTest, HandlesMultiVertexCoreRingAndBranchingCut) {
  // U = 6, f = 2: tree {0,1,2,3}, core {0,1}, ring {4,5}.
  // The unique min cut (value 53) consists of e2 = (1,3) and e3 = (2,4).
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
  EXPECT_EQ(cut[0], 2u);
  EXPECT_EQ(cut[1], 3u);

  EXPECT_TRUE(covered[0]);
  EXPECT_TRUE(covered[1]);
  EXPECT_FALSE(covered[2]);
  EXPECT_FALSE(covered[3]);
  EXPECT_FALSE(covered[4]);
  EXPECT_FALSE(covered[5]);
}

TEST(RunNaturalCutDetectionSequentialTest, WholeSmallGraphAbsorbedYieldsNoCuts) {
  // The whole graph fits into one tree: nothing to cut.
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
  // Limit the threads to force contention on the shared flags.
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, 8);

  const int len = 500;
  std::vector<NodeWeight> weights(len, 1);
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i + 1 < len; ++i) edges.push_back({static_cast<NodeID>(i), static_cast<NodeID>(i + 1), 1});
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  for (int trial = 0; trial < 10; ++trial) {
    std::vector<char> keep = run_natural_cut_detection(graph, NaturalCutParams{5, 1.0, 10.0, 2});
    ASSERT_EQ(keep.size(), static_cast<size_t>(len - 1));
    bool any_kept = false;
    for (char k : keep) any_kept = any_kept || k;
    EXPECT_TRUE(any_kept);
  }
}
