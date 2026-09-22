#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

using namespace mt_kahypar::filtering;

namespace {
// A simple random connected graph generator: a random spanning tree (so
// connectivity is guaranteed) plus extra random edges, unit vertex weights,
// random small integer edge weights.
FilterGraph make_random_connected_graph(size_t n, size_t extra_edges, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<NodeWeight> weights(n, 1);
  std::vector<EdgeListEntry> edges;
  std::vector<NodeID> shuffled(n);
  for (size_t i = 0; i < n; ++i) shuffled[i] = static_cast<NodeID>(i);
  std::shuffle(shuffled.begin(), shuffled.end(), rng);
  for (size_t i = 1; i < n; ++i) {
    std::uniform_int_distribution<size_t> pick_parent(0, i - 1);
    const NodeID parent = shuffled[pick_parent(rng)];
    edges.push_back({parent, shuffled[i], static_cast<EdgeWeight>(1 + rng() % 10)});
  }
  std::uniform_int_distribution<size_t> pick_node(0, n - 1);
  for (size_t i = 0; i < extra_edges; ++i) {
    const NodeID u = static_cast<NodeID>(pick_node(rng));
    const NodeID v = static_cast<NodeID>(pick_node(rng));
    if (u != v) edges.push_back({u, v, static_cast<EdgeWeight>(1 + rng() % 10)});
  }
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(FilteringPipelineTest, NoFragmentEverExceedsU) {
  for (uint64_t seed = 0; seed < 10; ++seed) {
    FilterGraph graph = make_random_connected_graph(200, 100, seed);
    FilteringParams params{20, 5, 1.0, 10.0, 2};
    FilteringResult result = run_filtering_pipeline(graph, params);

    for (NodeWeight size : result.fragment_size) {
      EXPECT_LE(size, params.U) << "seed=" << seed;
    }

    NodeWeight total_original = static_cast<NodeWeight>(graph.numNodes());  // unit weights
    NodeWeight total_fragments = 0;
    for (NodeWeight size : result.fragment_size) total_fragments += size;
    EXPECT_EQ(total_original, total_fragments) << "seed=" << seed;
  }
}

TEST(FilteringPipelineTest, EveryOriginalVertexGetsAFragment) {
  FilterGraph graph = make_random_connected_graph(50, 20, 123);
  FilteringResult result = run_filtering_pipeline(graph, FilteringParams{10, 5, 1.0, 10.0, 2});
  ASSERT_EQ(result.fragment_id.size(), graph.numNodes());
  for (NodeID id : result.fragment_id) EXPECT_LT(id, result.fragment_size.size());
}
