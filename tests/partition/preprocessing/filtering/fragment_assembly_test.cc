#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"

using namespace mt_kahypar::filtering;

TEST(AssembleFragmentsTest, NonKeptEdgesMergeIntoOneFragment) {
  // Part-1-output graph: a path 0-1-2-3. No edge is kept -> the whole path
  // is one fragment. part1_mapping is the identity (as if Part 1 did
  // nothing), i.e. 5 original vertices map onto a 4-vertex Part-1 graph with
  // one collision (originals 3 and 4 both map to Part-1 vertex 3).
  std::vector<NodeWeight> weights = {1, 1, 1, 2};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  std::vector<char> keep(3, 0);
  std::vector<NodeID> part1_mapping = {0, 1, 2, 3, 3};

  FilteringResult result = assemble_fragments(graph, keep, part1_mapping);

  ASSERT_EQ(result.fragment_id.size(), 5u);
  const NodeID frag = result.fragment_id[0];
  for (NodeID orig = 0; orig < 5; ++orig) EXPECT_EQ(result.fragment_id[orig], frag);
  ASSERT_EQ(result.fragment_size.size(), 1u);
  EXPECT_EQ(result.fragment_size[frag], 5u);  // 1+1+1+2 (original weights)
  EXPECT_TRUE(result.kept_edges.empty());
}

TEST(AssembleFragmentsTest, KeptEdgeSeparatesTwoFragments) {
  // Path 0-1-2-3, edge (1,2) [canonical id 1] is kept -> two fragments:
  // {0,1} and {2,3}.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  std::vector<char> keep = {0, 1, 0};
  std::vector<NodeID> part1_mapping = {0, 1, 2, 3};

  FilteringResult result = assemble_fragments(graph, keep, part1_mapping);

  EXPECT_EQ(result.fragment_id[0], result.fragment_id[1]);
  EXPECT_EQ(result.fragment_id[2], result.fragment_id[3]);
  EXPECT_NE(result.fragment_id[0], result.fragment_id[2]);
  ASSERT_EQ(result.fragment_size.size(), 2u);
  ASSERT_EQ(result.kept_edges.size(), 1u);
  EXPECT_TRUE((result.kept_edges[0] == std::make_pair(NodeID(1), NodeID(2))) ||
              (result.kept_edges[0] == std::make_pair(NodeID(2), NodeID(1))));
}
