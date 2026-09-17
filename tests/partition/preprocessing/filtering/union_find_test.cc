// tests/partition/preprocessing/filtering/union_find_test.cc
#include <gtest/gtest.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

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
  // Deliberately not asserting uf.setSize(0) here: setSize() is a
  // best-effort heuristic that can be permanently undercounted under
  // concurrent structural changes (see the class's doc comment) -- it is
  // never used for correctness-critical logic in this module, only
  // set-membership (find()/unite()), which the assertions above do check.
}

TEST(AtomicUnionFindTest, ConcurrentSwappedArgumentOrderNeverFormsACycle) {
  // Regression test for a race where two concurrent unite() calls resolving
  // to the same two already-formed, equal-size roots -- but with the
  // arguments passed in opposite order -- could both succeed in opposite
  // attach directions, forming a 2-cycle. A racing find() could then use
  // that cycle to fully un-merge two already-united elements. This is
  // exactly the failure mode a naive size-based (rather than id-based)
  // attach-orientation rule allows; repeated many times since it is a
  // timing-dependent race that a single trial has no guarantee of hitting.
  for (int trial = 0; trial < 200; ++trial) {
    AtomicUnionFind uf(4);
    uf.unite(0, 1);  // component A = {0, 1}, size 2
    uf.unite(2, 3);  // component B = {2, 3}, size 2 (tied with A)
    tbb::parallel_invoke(
      [&] { uf.unite(0, 2); },
      [&] { uf.unite(3, 1); });  // same two roots, swapped argument order
    // Regardless of scheduling, all four elements must end up in one set.
    const uint32_t root = uf.find(0);
    EXPECT_EQ(uf.find(1), root) << "trial " << trial;
    EXPECT_EQ(uf.find(2), root) << "trial " << trial;
    EXPECT_EQ(uf.find(3), root) << "trial " << trial;
  }
}
