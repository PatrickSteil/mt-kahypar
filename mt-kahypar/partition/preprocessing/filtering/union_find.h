#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace filtering {

// Lock-free union-find over a fixed universe of ids [0, n). Thread-safe for
// concurrent find()/unite() calls; the universe size is fixed at construction.
// Union by id (not by size -- see below) with path halving. `setSize` is a
// best-effort heuristic ONLY: under concurrent structural changes its
// undercount can be permanent, not merely transient (a root can be demoted
// by an unrelated concurrent unite() between being read and having a
// sibling's size folded into it). It must never be used for
// correctness-critical decisions -- only for coarse balancing/reporting.
// Nothing in this module relies on setSize() for correctness.
class AtomicUnionFind {
 public:
  explicit AtomicUnionFind(size_t n) : _parent(n), _size(n) {
    for (size_t i = 0; i < n; ++i) {
      _parent[i] = static_cast<uint32_t>(i);
      _size[i] = 1;
    }
  }

  uint32_t find(uint32_t x) const {
    uint32_t p = _parent[x].load(std::memory_order_relaxed);
    while (p != x) {
      const uint32_t gp = _parent[p].load(std::memory_order_relaxed);
      _parent[x].compare_exchange_weak(p, gp, std::memory_order_relaxed);
      x = p;
      p = _parent[x].load(std::memory_order_relaxed);
    }
    return x;
  }

  // Unions the sets containing a and b. Returns true if they were disjoint.
  bool unite(uint32_t a, uint32_t b) {
    for (;;) {
      uint32_t ra = find(a);
      uint32_t rb = find(b);
      if (ra == rb) return false;
      // Deterministic orientation: always attach the larger-id root under
      // the smaller-id root. This MUST NOT depend on the caller's argument
      // order, nor on a racy read of `_size` -- two concurrent unite() calls
      // that resolve to the same two pre-existing roots (e.g. unite(x,y) and
      // unite(y,x) racing on already-formed components, or two calls whose
      // relative size reads flip due to a third thread's concurrent update)
      // must always agree on which side attaches to which. Root ids don't
      // change while a node is still a root, so comparing them is race-free
      // in a way comparing `_size` is not. Getting this wrong lets both
      // CASes below succeed in opposite directions, forming a 2-cycle that a
      // racing find() can then use to fully un-merge two already-united
      // elements -- this is exactly the bug this comment exists to prevent
      // a future edit from reintroducing.
      if (ra > rb) std::swap(ra, rb);
      uint32_t expected = rb;
      if (_parent[rb].compare_exchange_strong(expected, ra, std::memory_order_relaxed)) {
        _size[ra].fetch_add(_size[rb].load(std::memory_order_relaxed), std::memory_order_relaxed);
        return true;
      }
      // Lost the race to another thread unioning rb; retry from scratch.
    }
  }

  size_t setSize(uint32_t x) const {
    return _size[find(x)].load(std::memory_order_relaxed);
  }

 private:
  mutable std::vector<parallel::IntegralAtomicWrapper<uint32_t>> _parent;
  std::vector<parallel::IntegralAtomicWrapper<uint64_t>> _size;
};

}  // namespace filtering
}  // namespace mt_kahypar
