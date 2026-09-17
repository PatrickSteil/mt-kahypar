// mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct SpanningForest {
  std::vector<NodeID> parent;       // BFS-tree parent; roots are their own parent
  std::vector<EdgeID> parent_edge;  // canonical edge id to parent; kInvalidEdge for roots
  std::vector<NodeID> bfs_order;    // visitation order, non-decreasing depth within each root
  std::vector<NodeID> roots;
};

// Computes connected components of `graph` in parallel via TBB parallel_for
// unioning every adjacency pair through an atomic union-find. Component ids
// are union-find representatives (small vertex ids), not a dense range. See
// design spec section 4.2 for why this simple approach (rather than
// Shiloach-Vishkin pointer jumping) is adequate here.
std::vector<NodeID> parallel_connected_components(const FilterGraph& graph);

// Same as above, but treats every edge e with excluded_edge[e] == true as
// absent. excluded_edge must have size graph.numEdges().
std::vector<NodeID> parallel_connected_components_excluding(
    const FilterGraph& graph, const std::vector<char>& excluded_edge);

// Builds a BFS spanning forest given a graph and its component labels.
// `component` MUST be the output of parallel_connected_components(graph) on
// this same graph -- NOT parallel_connected_components_excluding(graph,
// ...). The BFS below walks `graph`'s real, unmodified adjacency list; it
// has no notion of which edges an "excluding" component computation
// treated as absent. If `component` came from the excluding variant, two
// "components" it labels as separate can still be mutually reachable via
// the excluded edges in the real graph, so their concurrently-running BFS
// instances would walk into the same actual vertices -- a genuine
// unsynchronized data race on `forest.parent`/`forest.parent_edge` (not
// just a logic bug), since the safety of processing components in parallel
// depends entirely on each root's true BFS-reachable set in `graph` being
// exactly its own component, which only holds for the non-excluding
// variant. If a future task needs a spanning forest of an edge-excluded
// view, it needs its own edge-aware BFS, not this function.
//
// If `roots` is empty, one root per distinct component value is chosen
// automatically (its lowest-numbered vertex); otherwise `roots` must
// contain exactly one vertex per distinct component value (each within
// [0, graph.numNodes())) -- violating this reproduces the same class of
// race described above, since it's equivalent to supplying two "roots"
// inside one real connected component. Debug builds assert this
// precondition. BFS runs sequentially per component, but components are
// processed in parallel via TBB (safe under the precondition above: each
// root's BFS only ever reaches vertices in its own component, which are
// disjoint vertex sets across components, never touching another
// in-flight BFS's output indices).
SpanningForest build_spanning_forest(const FilterGraph& graph,
                                      const std::vector<NodeID>& component,
                                      const std::vector<NodeID>& roots);

}  // namespace filtering
}  // namespace mt_kahypar
