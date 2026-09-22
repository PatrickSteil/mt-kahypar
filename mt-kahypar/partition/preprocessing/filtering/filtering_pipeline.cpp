#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

namespace mt_kahypar {
namespace filtering {

namespace {

// Natural-cut detection (Part 2) bounds only the small "reachable" core side
// of each seed's local min cut (design spec section 5) -- that side has
// weight <= alpha*U <= U by construction of the BFS growth. Everything else
// a seed's BFS visits (its own tree beyond the cut, and every vertex that
// happens to be swept up as *someone else's* tree/ring member before ever
// getting a turn as a seed itself) is left for other seeds to resolve. On a
// graph with many small "natural" cuts near cores rather than near a full
// tree of weight ~U (e.g. a sparse, near-tree-like graph), most vertices
// never end up on the reachable side of any single seed's cut in one sweep,
// so they default to being unioned together into one unbounded "leftover"
// component that can exceed U -- a real bug found by
// filtering_invariant_test.cc's NoFragmentEverExceedsU, not merely a
// theoretical risk.
//
// Fixed here, at the orchestrator level, by iterative refinement: after a
// sweep, any resulting fragment that still exceeds U is re-analyzed as its
// own fresh induced subgraph, giving every one of its vertices a new chance
// to be its own seed. This is guaranteed to make progress each round: a
// fragment's total weight exceeding U (== alpha*U here) means no single BFS
// tree can absorb it whole, so its ring is always non-empty and at least one
// new edge is found and marked kept, strictly increasing the total number of
// kept edges (bounded by the graph's edge count) -- so the loop terminates.
void refine_oversized_component(const FilterGraph& g, const std::vector<NodeID>& vertices,
                                 const NaturalCutParams& params, std::vector<char>& keep) {
  if (vertices.size() <= 1) return;  // a single vertex cannot be split further

  std::vector<NodeID> local_of(g.numNodes(), kInvalidNode);
  for (size_t i = 0; i < vertices.size(); ++i) {
    local_of[vertices[i]] = static_cast<NodeID>(i);
  }

  std::vector<NodeWeight> local_weights(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i) {
    local_weights[i] = g.node_weight[vertices[i]];
  }

  std::vector<EdgeListEntry> local_edges;
  std::vector<EdgeID> local_edge_to_global;
  for (NodeID v : vertices) {
    for (EdgeID pos = g.node_begin[v]; pos < g.node_begin[v + 1]; ++pos) {
      const NodeID u = g.adj[pos];
      if (u <= v) continue;  // process each undirected edge once
      if (local_of[u] == kInvalidNode) continue;  // outside this component
      const EdgeID global_edge = g.adj_edge[pos];
      local_edges.push_back({local_of[v], local_of[u], g.edge_weight[global_edge]});
      local_edge_to_global.push_back(global_edge);
    }
  }
  if (local_edges.empty()) return;  // no internal edges left to cut (shouldn't happen for a
                                     // connected component of size > 1, but guard anyway)

  FilterGraph local_graph = build_csr_from_edge_list(local_edges, local_weights);
  std::vector<char> local_keep = run_natural_cut_detection(local_graph, params);

  // build_csr_from_edge_list assigns canonical edge ids in insertion order,
  // so local_keep[i] corresponds directly to local_edge_to_global[i].
  for (size_t i = 0; i < local_keep.size(); ++i) {
    if (local_keep[i]) keep[local_edge_to_global[i]] = 1;
  }
}

}  // namespace

FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params) {
  ContractionResult tiny_cut_result =
      run_tiny_cut_detection(graph, TinyCutParams{params.U, params.tau});
  const FilterGraph& g = tiny_cut_result.graph;

  const NaturalCutParams natural_params{params.U, params.alpha, params.f, params.coverage};
  std::vector<char> keep = run_natural_cut_detection(g, natural_params);

  // One representative original vertex per g-vertex, used below to recover
  // each g-vertex's current fragment id from a FilteringResult (which is
  // indexed by original vertex id).
  std::vector<NodeID> rep_orig(g.numNodes(), kInvalidNode);
  for (size_t orig = 0; orig < tiny_cut_result.mapping.size(); ++orig) {
    NodeID& rep = rep_orig[tiny_cut_result.mapping[orig]];
    if (rep == kInvalidNode) rep = static_cast<NodeID>(orig);
  }

  const size_t max_iterations = g.numEdges() + 1;
  for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
    FilteringResult probe = assemble_fragments(g, keep, tiny_cut_result.mapping);

    bool any_oversized = false;
    for (NodeWeight size : probe.fragment_size) {
      if (size > params.U) { any_oversized = true; break; }
    }
    if (!any_oversized) break;

    std::vector<std::vector<NodeID>> members(probe.fragment_size.size());
    for (NodeID gv = 0; gv < static_cast<NodeID>(g.numNodes()); ++gv) {
      members[probe.fragment_id[rep_orig[gv]]].push_back(gv);
    }
    for (size_t fid = 0; fid < members.size(); ++fid) {
      if (probe.fragment_size[fid] > params.U) {
        refine_oversized_component(g, members[fid], natural_params, keep);
      }
    }
  }

  return assemble_fragments(g, keep, tiny_cut_result.mapping);
}

}  // namespace filtering
}  // namespace mt_kahypar
