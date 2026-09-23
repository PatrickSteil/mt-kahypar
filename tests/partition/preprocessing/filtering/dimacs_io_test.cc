#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "mt-kahypar/partition/preprocessing/filtering/dimacs_io.h"

using namespace mt_kahypar::filtering;

namespace {
void write_file(const std::string& path, const std::string& content) {
  std::ofstream out(path);
  out << content;
}
}  // namespace

TEST(DimacsIoTest, ParsesSimpleGraph) {
  const std::string path = "dimacs_io_test_simple.gr";
  write_file(path,
      "c comment line\n"
      "p sp 3 6\n"
      "a 1 2 4\n"
      "a 2 1 4\n"
      "a 2 3 9\n"
      "a 3 2 9\n"
      "a 1 3 2\n"
      "a 3 1 2\n");

  FilterGraph graph = read_dimacs_graph(path);
  std::remove(path.c_str());

  ASSERT_EQ(graph.numNodes(), 3u);
  ASSERT_EQ(graph.numEdges(), 3u);
  for (NodeID v = 0; v < 3; ++v) {
    EXPECT_EQ(graph.degree(v), 2u);
    EXPECT_EQ(graph.node_weight[v], 1u);
  }
  for (EdgeID e = 0; e < graph.numEdges(); ++e) {
    EXPECT_EQ(graph.edge_weight[e], 1);
  }
}

TEST(DimacsIoTest, IgnoresEdgeWeightsEvenWhenInconsistent) {
  const std::string path = "dimacs_io_test_weights.gr";
  write_file(path, "p sp 2 2\na 1 2 4\na 2 1 5\n");
  FilterGraph graph = read_dimacs_graph(path);
  std::remove(path.c_str());

  ASSERT_EQ(graph.numEdges(), 1u);
  EXPECT_EQ(graph.edge_weight[0], 1);
}

TEST(DimacsIoTest, ThrowsOnMissingFile) {
  EXPECT_THROW(read_dimacs_graph("dimacs_io_test_does_not_exist.gr"), std::runtime_error);
}
