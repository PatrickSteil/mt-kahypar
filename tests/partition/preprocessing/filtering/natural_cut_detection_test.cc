#include <gtest/gtest.h>

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
