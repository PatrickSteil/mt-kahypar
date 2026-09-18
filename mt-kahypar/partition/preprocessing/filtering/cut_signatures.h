#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

struct EdgeSignatures {
  std::vector<unsigned __int128> signature;  // indexed by canonical edge id
  std::vector<char> is_tree_edge;             // indexed by canonical edge id
};

// Assigns a random 128-bit label to every non-tree edge of `forest`, and for
// every tree edge computes the XOR of labels of non-tree edges whose
// fundamental cycle covers it. A non-tree edge's own signature is its own
// label. See design spec section 4.4/4.5 (Pritchard & Thurimella, ACM TALG
// 7(4), 2011) for the correctness argument. Monte Carlo with negligible
// one-sided error at 128 bits; callers verify candidate results locally
// (Task 7) rather than relying on the labels alone.
EdgeSignatures compute_edge_signatures(const FilterGraph& graph, const SpanningForest& forest);

// A tree edge with signature 0 is a bridge (its fundamental-cycle coverage
// set is empty); non-tree edges are never bridges, since the spanning tree
// alone already connects the graph without them.
std::vector<char> compute_bridges(const FilterGraph& graph, const EdgeSignatures& sigs);

// Finds all 2-edge-cut equivalence classes of `graph` (design spec section
// 4.4): buckets every non-bridge edge by its signature (tree edges by their
// aggregated label, non-tree edges by their own label -- see
// compute_edge_signatures), then verifies each candidate class of size >= 2
// by removing every edge in the class at once and checking that two of its
// endpoints land in different connected components. This is a sufficient
// check (not merely a one-pair spot-check): signature equality already
// implies identical fundamental-cycle coverage sets for every edge in the
// bucket (up to an astronomically rare 128-bit collision), which by the
// cographic-matroid argument in the design spec means EVERY pair within a
// genuine bucket forms a valid 2-cut -- so confirming any single
// representative pair is disconnected is enough to confirm the whole
// bucket. Verified classes are returned; unverified (collision) candidates
// are dropped. Note this deliberately does NOT special-case "the class
// happens to form a simple cycle" or similar shape-based heuristics: a
// bucket like a whole triangle's 3 edges (which all share one signature,
// since the triangle's one chord's fundamental cycle covers both of its
// tree edges) is a genuine, correct class -- every pair of a triangle's
// edges is a valid 2-cut, since a triangle's own minimum edge cut is 2, not
// 3 (see BridgeIsExcludedButEachTriangleIsItsOwnClass in the test file,
// which replaced an earlier, incorrect "no classes at all" expectation for
// exactly this shape).
std::vector<std::vector<EdgeID>> find_two_edge_cut_classes(const FilterGraph& graph,
                                                            const EdgeSignatures& sigs);

}  // namespace filtering
}  // namespace mt_kahypar
