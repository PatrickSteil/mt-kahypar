#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace filtering {

// Lock-free union-find (union by id, path halving). setSize() is only
// approximate under concurrent unite() calls.
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

  // Returns true if a and b were in different sets.
  bool unite(uint32_t a, uint32_t b) {
    for (;;) {
      uint32_t ra = find(a);
      uint32_t rb = find(b);
      if (ra == rb) return false;
      // Orient by id, not by size: concurrent unite() calls on the same roots
      // must agree on the direction, otherwise they can create a cycle.
      if (ra > rb) std::swap(ra, rb);
      uint32_t expected = rb;
      if (_parent[rb].compare_exchange_strong(expected, ra, std::memory_order_relaxed)) {
        _size[ra].fetch_add(_size[rb].load(std::memory_order_relaxed), std::memory_order_relaxed);
        return true;
      }
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
