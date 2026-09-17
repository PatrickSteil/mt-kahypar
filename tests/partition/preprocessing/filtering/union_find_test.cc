// tests/partition/preprocessing/filtering/union_find_test.cc
#include <gtest/gtest.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

using namespace mt_kahypar::filtering;

TEST(AtomicUnionFindTest, SingletonsStartDisjoint) {
  AtomicUnionFind uf(5);
  for (uint32_t i = 0; i < 5; ++i) {
    EXPECT_EQ(uf.find(i), i);
    EXPECT_EQ(uf.setSize(i), 1u);
  }
}

TEST(AtomicUnionFindTest, UnitesTransitively) {
  AtomicUnionFind uf(5);
  EXPECT_TRUE(uf.unite(0, 1));
  EXPECT_TRUE(uf.unite(1, 2));
  EXPECT_EQ(uf.find(0), uf.find(2));
  EXPECT_EQ(uf.setSize(0), 3u);
  EXPECT_FALSE(uf.unite(0, 2));
  EXPECT_NE(uf.find(0), uf.find(3));
}

TEST(AtomicUnionFindTest, ConcurrentUnionsFormOneSet) {
  const uint32_t n = 1000;
  AtomicUnionFind uf(n);
  tbb::parallel_for(tbb::blocked_range<uint32_t>(0, n - 1),
    [&](const tbb::blocked_range<uint32_t>& range) {
      for (uint32_t i = range.begin(); i < range.end(); ++i) uf.unite(i, i + 1);
    });
  const uint32_t root = uf.find(0);
  for (uint32_t i = 1; i < n; ++i) EXPECT_EQ(uf.find(i), root);
  EXPECT_EQ(uf.setSize(0), n);
}
