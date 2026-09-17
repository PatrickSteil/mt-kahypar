# PUNCH-style graph filtering module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement PUNCH's filtering phase (tiny-cut + natural-cut detection) as a standalone, parallel module operating on plain CSR weighted graphs, with its own CLI and tests, validated against small DIMACS road instances.

**Architecture:** A new `mt_kahypar::filtering` namespace under `mt-kahypar/partition/preprocessing/filtering/`, built bottom-up: generic primitives (CSR graph, union-find, contraction, connectivity, cut signatures) first, then Part 1 (tiny-cut), then Part 2 (natural-cut), then Part 3 (assembly) + CLI + validation. Every task is TDD: failing test, minimal implementation, passing test, commit.

**Tech Stack:** C++17, TBB (oneTBB, already a hard build dependency), GoogleTest/GoogleMock (already vendored via `KAHYPAR_ENABLE_TESTING`), `ds::Array` for persistent CSR storage (matching existing codebase convention), plain `std::vector`/`std::unordered_map` for transient working structures.

**Spec:** `docs/superpowers/specs/2026-09-17-punch-filtering-design.md`

## Global Constraints

- Module headers must never include `ds::StaticGraph`, `ds::StaticHypergraph`, `ds::Graph`, `HypernodeID`, `Context`, or anything from `partition/` outside `filtering/`. Only TBB, the standard library, and `parallel/atomic_wrapper.h`, `parallel/stl/thread_locals.h`, `parallel/stl/scalable_vector.h`, `datastructures/array.h` are allowed as Mt-KaHyPar dependencies.
- `NodeID`/`EdgeID` are `uint32_t`, `NodeWeight` is `uint64_t`, `EdgeWeight` is `int64_t` (spec section 3).
- Default parameters: `tau = 5`, `alpha = 1.0`, `f = 10`, `coverage (C) = 2` (spec sections 4.2, 5).
- Build registration follows the existing pattern: sources are added to the `PreprocessingSources` list in `mt-kahypar/partition/preprocessing/CMakeLists.txt` (no new nested `CMakeLists.txt`, matching how `community_detection/` is wired in) via `target_sources(MtKaHyPar-Sources INTERFACE ...)` and `target_sources(MtKaHyPar-ToolsSources INTERFACE ...)`. Tests are added to `tests/partition/preprocessing/CMakeLists.txt` via `target_sources(mtkahypar_tests PRIVATE ...)`.
- **Correction from the spec during planning:** spec section 4.2 ("build spanning tree, contract subtree ≤ U top-down") actually requires the **bridge tree** (nodes = 2-edge-connected blocks, edges = bridges), not a spanning tree of the whole graph — a plain spanning tree of a connected graph has no substructure to contract. Task 8 below implements the corrected version. Pass 1 and pass 3 (2-edge-cuts) both need the same spanning-tree-plus-random-XOR-label machinery, factored into a shared `cut_signatures` module (Task 6).
- **Correction found during pre-flight review (before any task was dispatched):** Task 8's original draft unconditionally merged every vertex of every block into one vertex per block, and translated tau-merges through a separate block-level union-find — both wrong. The first would collapse blocks that were never chosen for contraction (destroying structure later passes need); the second could silently reference an uninitialized representative and crash. Worse, a naive per-merge tau-merge check (`child_subtree + parent's own weight <= U`) allows two *separate* small children of the same parent to each individually pass the check while their *combined* weight, once both are merged into that parent, exceeds `U` — a real violation of the hard `U`-invariant that Task 17's `filtering_invariant_test.cc` exists specifically to catch. Task 8's implementation below is the corrected version: block representatives are computed unconditionally, vertex-level unions happen directly (no block-level union-find detour), and a `tau_merged_extra_weight` accumulator per parent block prevents the cascading overflow. Task 8's test `TauMergePreventsCascadingOverflow` exercises exactly this case.
- **Deferred from the spec:** the "two-components-at-a-time" traversal optimization for processing 2-edge-cut classes (spec section 4.4 step 8) is not implemented in this pass. Each class is processed via a direct connected-components-excluding-edges call instead. This is simpler and correct; flagged as a follow-up optimization if profiling on real DIMACS instances shows it's needed (real road networks are expected to produce a small number of 2-cut classes per contraction stage).

---

## File Structure

```
mt-kahypar/partition/preprocessing/filtering/
    filter_graph.h / .cpp        # FilterGraph CSR struct + builder + edge-endpoint helper
    dimacs_io.h / .cpp           # DIMACS .gr reader
    union_find.h                 # AtomicUnionFind
    graph_contraction.h / .cpp   # contract_graph()
    parallel_connectivity.h / .cpp # connected components (+ excluding edges), spanning forest
    cut_signatures.h / .cpp      # shared PT XOR-signature machinery: bridges + 2-cut classes
    tiny_cut_detection.h / .cpp  # Part 1: passes 1-3 + orchestrator
    dinic_max_flow.h / .cpp      # local s-t min-cut solver
    natural_cut_detection.h / .cpp # Part 2: BFS core/ring, sequential + parallel seed loops
    fragment_assembly.h / .cpp   # Part 3: final union-find + FilteringResult
    filtering_pipeline.h / .cpp  # orchestrates Parts 1-3, documented extension point

tools/
    punch_filter.cc              # standalone CLI

tests/partition/preprocessing/filtering/
    filter_graph_test.cc
    dimacs_io_test.cc
    union_find_test.cc
    graph_contraction_test.cc
    parallel_connectivity_test.cc
    cut_signatures_test.cc
    tiny_cut_detection_test.cc
    dinic_max_flow_test.cc
    natural_cut_detection_test.cc
    fragment_assembly_test.cc
    filtering_invariant_test.cc
```

---

### Task 1: Module scaffolding — FilterGraph + build wiring

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/filter_graph.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/filter_graph.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/filter_graph_test.cc`

**Interfaces:**
- Produces: `NodeID`, `EdgeID`, `NodeWeight`, `EdgeWeight`, `kInvalidNode`, `kInvalidEdge`, `FilterGraph{node_begin, adj, adj_edge, node_weight, edge_weight, numNodes(), numEdges(), degree(v)}`, `EdgeListEntry{u,v,weight}`, `build_csr_from_edge_list(edges, node_weights) -> FilterGraph`, `compute_edge_endpoints(graph) -> std::vector<std::pair<NodeID,NodeID>>` (indexed by canonical edge id).

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/filter_graph_test.cc
#include <gtest/gtest.h>
#include <algorithm>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

using namespace mt_kahypar::filtering;

TEST(FilterGraphTest, BuildsTriangle) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}, {0, 2, 3}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ASSERT_EQ(graph.numNodes(), 3u);
  ASSERT_EQ(graph.numEdges(), 3u);
  EXPECT_EQ(graph.degree(0), 2u);
  EXPECT_EQ(graph.degree(1), 2u);
  EXPECT_EQ(graph.degree(2), 2u);

  std::vector<std::pair<NodeID, EdgeWeight>> incident_to_0;
  for (EdgeID pos = graph.node_begin[0]; pos < graph.node_begin[1]; ++pos) {
    incident_to_0.emplace_back(graph.adj[pos], graph.edge_weight[graph.adj_edge[pos]]);
  }
  std::sort(incident_to_0.begin(), incident_to_0.end());
  EXPECT_EQ(incident_to_0[0], std::make_pair(NodeID(1), EdgeWeight(5)));
  EXPECT_EQ(incident_to_0[1], std::make_pair(NodeID(2), EdgeWeight(3)));
}

TEST(FilterGraphTest, HandlesEmptyGraph) {
  std::vector<NodeWeight> weights = {1};
  std::vector<EdgeListEntry> edges;
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  ASSERT_EQ(graph.numNodes(), 1u);
  ASSERT_EQ(graph.numEdges(), 0u);
  EXPECT_EQ(graph.degree(0), 0u);
}

TEST(FilterGraphTest, ComputesEdgeEndpoints) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  auto endpoints = compute_edge_endpoints(graph);
  ASSERT_EQ(endpoints.size(), 2u);
  auto has_pair = [&](NodeID a, NodeID b) {
    for (auto& p : endpoints) {
      if ((p.first == a && p.second == b) || (p.first == b && p.second == a)) return true;
    }
    return false;
  };
  EXPECT_TRUE(has_pair(0, 1));
  EXPECT_TRUE(has_pair(1, 2));
}
```

- [ ] **Step 2: Create the filtering test directory and register it**

Edit `tests/partition/preprocessing/CMakeLists.txt`, adding to the existing `target_sources(mtkahypar_tests PRIVATE ...)` block for this directory:

```cmake
target_sources(mtkahypar_tests PRIVATE
        louvain_test.cc
        filtering/filter_graph_test.cc
        )
```

- [ ] **Step 3: Run test to verify it fails (build error: no such header)**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `filter_graph.h: No such file or directory`

- [ ] **Step 4: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/filter_graph.h
#pragma once

#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "mt-kahypar/datastructures/array.h"

namespace mt_kahypar {
namespace filtering {

using NodeID = uint32_t;
using EdgeID = uint32_t;
using NodeWeight = uint64_t;
using EdgeWeight = int64_t;

constexpr NodeID kInvalidNode = std::numeric_limits<NodeID>::max();
constexpr EdgeID kInvalidEdge = std::numeric_limits<EdgeID>::max();

// Plain CSR representation of an undirected weighted graph, independent of
// Mt-KaHyPar's hypergraph/partitioning types (see design spec section 3).
struct FilterGraph {
  ds::Array<EdgeID> node_begin;       // size numNodes() + 1
  ds::Array<NodeID> adj;              // size 2 * numEdges(), neighbor per half-edge
  ds::Array<EdgeID> adj_edge;         // size 2 * numEdges(), half-edge -> canonical edge id
  ds::Array<NodeWeight> node_weight;  // size numNodes()
  ds::Array<EdgeWeight> edge_weight;  // size numEdges()

  size_t numNodes() const { return node_weight.size(); }
  size_t numEdges() const { return edge_weight.size(); }
  size_t degree(NodeID v) const { return node_begin[v + 1] - node_begin[v]; }
};

// One undirected edge, used while building a FilterGraph from an edge list.
struct EdgeListEntry {
  NodeID u;
  NodeID v;
  EdgeWeight weight;
};

// Builds a FilterGraph in CSR form. node_weights.size() determines the number
// of nodes; every edge must reference node ids in [0, node_weights.size()).
// Callers must not pass self-loops (u == v) or duplicate edges between the
// same pair; graph_contraction.cpp is responsible for merging/dropping those
// when a contraction would create them.
FilterGraph build_csr_from_edge_list(const std::vector<EdgeListEntry>& edges,
                                      const std::vector<NodeWeight>& node_weights);

// Recovers the (u, v) endpoint pair for every canonical edge id, by scanning
// the adjacency structure once (O(numNodes() + numEdges())).
std::vector<std::pair<NodeID, NodeID>> compute_edge_endpoints(const FilterGraph& graph);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/filter_graph.cpp
#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

FilterGraph build_csr_from_edge_list(const std::vector<EdgeListEntry>& edges,
                                      const std::vector<NodeWeight>& node_weights) {
  const size_t n = node_weights.size();
  const size_t m = edges.size();

  FilterGraph graph;
  graph.node_begin.resize(n + 1, 0);
  graph.adj.resize(2 * m, kInvalidNode);
  graph.adj_edge.resize(2 * m, kInvalidEdge);
  graph.node_weight.resize(n, 0);
  graph.edge_weight.resize(m, 0);

  for (size_t v = 0; v < n; ++v) graph.node_weight[v] = node_weights[v];
  for (size_t e = 0; e < m; ++e) graph.edge_weight[e] = edges[e].weight;

  std::vector<size_t> degree(n, 0);
  for (const EdgeListEntry& e : edges) {
    ++degree[e.u];
    ++degree[e.v];
  }

  EdgeID offset = 0;
  for (size_t v = 0; v < n; ++v) {
    graph.node_begin[v] = offset;
    offset += static_cast<EdgeID>(degree[v]);
  }
  graph.node_begin[n] = offset;

  std::vector<EdgeID> cursor(n);
  for (size_t v = 0; v < n; ++v) cursor[v] = graph.node_begin[v];

  for (EdgeID e = 0; e < static_cast<EdgeID>(m); ++e) {
    const NodeID u = edges[e].u;
    const NodeID v = edges[e].v;
    graph.adj[cursor[u]] = v;
    graph.adj_edge[cursor[u]] = e;
    ++cursor[u];
    graph.adj[cursor[v]] = u;
    graph.adj_edge[cursor[v]] = e;
    ++cursor[v];
  }

  return graph;
}

std::vector<std::pair<NodeID, NodeID>> compute_edge_endpoints(const FilterGraph& graph) {
  std::vector<std::pair<NodeID, NodeID>> endpoints(
      graph.numEdges(), {kInvalidNode, kInvalidNode});
  for (size_t v = 0; v < graph.numNodes(); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const EdgeID e = graph.adj_edge[pos];
      if (endpoints[e].first == kInvalidNode) {
        endpoints[e].first = static_cast<NodeID>(v);
      } else {
        endpoints[e].second = static_cast<NodeID>(v);
      }
    }
  }
  return endpoints;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 6: Register the source file**

Edit `mt-kahypar/partition/preprocessing/CMakeLists.txt`:

```cmake
set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp
        filtering/filter_graph.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${PreprocessingSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${PreprocessingSources})
```

(Note: this also adds `MtKaHyPar-ToolsSources` registration for `PreprocessingSources`, which the file didn't have before — needed so `tools/punch_filter.cc`, added in Task 18, can link against this module.)

- [ ] **Step 7: Build and run the test**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=FilterGraphTest.*`
Expected: PASS (3 tests)

- [ ] **Step 8: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/filter_graph.h \
        mt-kahypar/partition/preprocessing/filtering/filter_graph.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/filter_graph_test.cc
git commit -m "feat: add FilterGraph CSR struct for the PUNCH filtering module

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: DIMACS `.gr` reader

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/dimacs_io.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/dimacs_io.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/dimacs_io_test.cc`

**Interfaces:**
- Consumes: `FilterGraph`, `EdgeListEntry`, `build_csr_from_edge_list` (Task 1).
- Produces: `read_dimacs_graph(path) -> FilterGraph`.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/dimacs_io_test.cc
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
}

TEST(DimacsIoTest, ThrowsOnInconsistentWeights) {
  const std::string path = "dimacs_io_test_bad.gr";
  write_file(path, "p sp 2 2\na 1 2 4\na 2 1 5\n");
  EXPECT_THROW(read_dimacs_graph(path), std::runtime_error);
  std::remove(path.c_str());
}

TEST(DimacsIoTest, ThrowsOnMissingFile) {
  EXPECT_THROW(read_dimacs_graph("dimacs_io_test_does_not_exist.gr"), std::runtime_error);
}
```

- [ ] **Step 2: Register the test file**

Edit `tests/partition/preprocessing/CMakeLists.txt`:

```cmake
target_sources(mtkahypar_tests PRIVATE
        louvain_test.cc
        filtering/filter_graph_test.cc
        filtering/dimacs_io_test.cc
        )
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `dimacs_io.h: No such file or directory`

- [ ] **Step 4: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/dimacs_io.h
#pragma once

#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

// Reads a 9th DIMACS Implementation Challenge ".gr" graph file. Format: a
// header line "p sp <n> <m>" (n nodes, m directed arcs -- each undirected
// edge appears once per direction), followed by "m" arc lines "a <u> <v> <w>"
// (1-indexed). Comment lines start with 'c' and are ignored. Vertex weights
// default to 1. Self-loops in the file are dropped. Throws std::runtime_error
// if the file cannot be opened, has no header, references an out-of-range
// node id, or if the two directions of an undirected edge disagree on weight.
FilterGraph read_dimacs_graph(const std::string& path);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/dimacs_io.cpp
#include "mt-kahypar/partition/preprocessing/filtering/dimacs_io.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace mt_kahypar {
namespace filtering {

namespace {
uint64_t undirected_key(NodeID u, NodeID v) {
  if (u > v) std::swap(u, v);
  return (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
}
}  // namespace

FilterGraph read_dimacs_graph(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("could not open DIMACS graph file: " + path);

  size_t n = 0;
  size_t num_arcs = 0;
  bool header_seen = false;
  std::unordered_map<uint64_t, EdgeListEntry> edges_by_key;

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == 'c') continue;
    std::istringstream iss(line);
    char tag;
    iss >> tag;
    if (tag == 'p') {
      std::string problem_type;
      iss >> problem_type >> n >> num_arcs;
      header_seen = true;
    } else if (tag == 'a') {
      if (!header_seen) {
        throw std::runtime_error("DIMACS file has an arc line before the 'p' header: " + path);
      }
      size_t u1 = 0, v1 = 0;
      long long w = 0;
      iss >> u1 >> v1 >> w;
      if (u1 == 0 || v1 == 0 || u1 > n || v1 > n) {
        throw std::runtime_error("DIMACS arc references an out-of-range node id: " + path);
      }
      const NodeID u = static_cast<NodeID>(u1 - 1);
      const NodeID v = static_cast<NodeID>(v1 - 1);
      if (u == v) continue;
      const uint64_t key = undirected_key(u, v);
      auto it = edges_by_key.find(key);
      if (it == edges_by_key.end()) {
        edges_by_key.emplace(key, EdgeListEntry{u, v, static_cast<EdgeWeight>(w)});
      } else if (it->second.weight != static_cast<EdgeWeight>(w)) {
        throw std::runtime_error(
            "DIMACS file has inconsistent weights for the two directions of an undirected edge: " + path);
      }
    }
  }

  if (!header_seen) throw std::runtime_error("DIMACS file has no 'p' header: " + path);

  std::vector<EdgeListEntry> edges;
  edges.reserve(edges_by_key.size());
  for (auto& [key, edge] : edges_by_key) edges.push_back(edge);

  std::vector<NodeWeight> node_weights(n, 1);
  return build_csr_from_edge_list(edges, node_weights);
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 6: Register the source file**

Edit `mt-kahypar/partition/preprocessing/CMakeLists.txt`, adding `filtering/dimacs_io.cpp` to `PreprocessingSources`.

- [ ] **Step 7: Build and run the test**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=DimacsIoTest.*`
Expected: PASS (3 tests)

- [ ] **Step 8: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/dimacs_io.h \
        mt-kahypar/partition/preprocessing/filtering/dimacs_io.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/dimacs_io_test.cc
git commit -m "feat: add DIMACS .gr reader for the filtering module

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: AtomicUnionFind

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/union_find.h` (header-only)
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/union_find_test.cc`

**Interfaces:**
- Produces: `AtomicUnionFind(n)`, `.find(x) -> uint32_t`, `.unite(a,b) -> bool`, `.setSize(x) -> size_t`.

- [ ] **Step 1: Write the failing test**

```cpp
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
```

- [ ] **Step 2: Register the test file**

Edit `tests/partition/preprocessing/CMakeLists.txt`, add `filtering/union_find_test.cc`.

- [ ] **Step 3: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `union_find.h: No such file or directory`

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/union_find.h
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
```

- [ ] **Step 5: Build and run the test**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=AtomicUnionFindTest.*`
Expected: PASS (4 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/union_find.h \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/union_find_test.cc
git commit -m "feat: add lock-free AtomicUnionFind for the filtering module

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: `contract_graph`

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/graph_contraction.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/graph_contraction.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/graph_contraction_test.cc`

**Interfaces:**
- Consumes: `FilterGraph`, `build_csr_from_edge_list` (Task 1), `AtomicUnionFind` (Task 3).
- Produces: `ContractionResult{graph, mapping}`, `contract_graph(graph, uf) -> ContractionResult`. `mapping[v]` (indexed by *input*-graph node id) gives the node id in the *output* graph.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/graph_contraction_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

using namespace mt_kahypar::filtering;

namespace {
FilterGraph make_path_graph() {
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 10}, {1, 2, 20}, {2, 3, 30}};
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(GraphContractionTest, MergesParallelEdgesAndSumsWeights) {
  FilterGraph graph = make_path_graph();
  AtomicUnionFind uf(4);
  uf.unite(0, 1);
  uf.unite(2, 3);

  ContractionResult result = contract_graph(graph, uf);
  ASSERT_EQ(result.graph.numNodes(), 2u);
  ASSERT_EQ(result.graph.numEdges(), 1u);
  EXPECT_EQ(result.graph.edge_weight[0], 20);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[2], result.mapping[3]);
  EXPECT_NE(result.mapping[0], result.mapping[2]);
  EXPECT_EQ(result.graph.node_weight[result.mapping[0]], 2u);
  EXPECT_EQ(result.graph.node_weight[result.mapping[2]], 2u);
}

TEST(GraphContractionTest, MergesTriangleIntoParallelEdgeThenDropsSelfLoop) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 5}, {1, 2, 7}, {0, 2, 3}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  AtomicUnionFind uf(3);
  uf.unite(1, 2);

  ContractionResult result = contract_graph(graph, uf);
  ASSERT_EQ(result.graph.numNodes(), 2u);
  ASSERT_EQ(result.graph.numEdges(), 1u);
  EXPECT_EQ(result.graph.edge_weight[0], 8);  // 5 (0-1) + 3 (0-2) merged
}

TEST(GraphContractionTest, ContractingEverythingLeavesNoEdges) {
  FilterGraph graph = make_path_graph();
  AtomicUnionFind uf(4);
  uf.unite(0, 1);
  uf.unite(1, 2);
  uf.unite(2, 3);

  ContractionResult result = contract_graph(graph, uf);
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.graph.numEdges(), 0u);
  EXPECT_EQ(result.graph.node_weight[0], 4u);
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/graph_contraction_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `graph_contraction.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/graph_contraction.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

struct ContractionResult {
  FilterGraph graph;
  std::vector<NodeID> mapping;  // input node id -> output node id
};

// Contracts `graph` according to the equivalence classes of `uf` (a
// union-find over graph.numNodes() ids). Vertex weights of a class are
// summed; parallel edges created by contraction are merged by summing their
// weights; self-loops created by contraction are dropped.
ContractionResult contract_graph(const FilterGraph& graph, const AtomicUnionFind& uf);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/graph_contraction.cpp
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

#include <map>

namespace mt_kahypar {
namespace filtering {

ContractionResult contract_graph(const FilterGraph& graph, const AtomicUnionFind& uf) {
  const size_t n = graph.numNodes();

  std::vector<NodeID> mapping(n, kInvalidNode);
  NodeID next_id = 0;
  for (size_t v = 0; v < n; ++v) {
    const uint32_t rep = uf.find(static_cast<uint32_t>(v));
    if (mapping[rep] == kInvalidNode) mapping[rep] = next_id++;
  }
  for (size_t v = 0; v < n; ++v) {
    mapping[v] = mapping[uf.find(static_cast<uint32_t>(v))];
  }

  std::vector<NodeWeight> node_weights(next_id, 0);
  for (size_t v = 0; v < n; ++v) node_weights[mapping[v]] += graph.node_weight[v];

  // Merge parallel edges (sum weights), drop self-loops created by contraction.
  // std::map keeps this simple/correct; can become a parallel hash
  // aggregation later if profiling shows it matters (tiny-cut detection is a
  // cheap stage of the pipeline overall -- see design spec section 4).
  std::map<std::pair<NodeID, NodeID>, EdgeWeight> merged_edges;
  for (NodeID v = 0; v < static_cast<NodeID>(n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;  // visit every undirected edge exactly once
      const NodeID mu = mapping[v];
      const NodeID mv = mapping[other];
      if (mu == mv) continue;  // self-loop created by contraction
      const std::pair<NodeID, NodeID> key = mu < mv ? std::make_pair(mu, mv) : std::make_pair(mv, mu);
      merged_edges[key] += graph.edge_weight[graph.adj_edge[pos]];
    }
  }

  std::vector<EdgeListEntry> edges;
  edges.reserve(merged_edges.size());
  for (auto& [endpoints, weight] : merged_edges) {
    edges.push_back(EdgeListEntry{endpoints.first, endpoints.second, weight});
  }

  ContractionResult result;
  result.graph = build_csr_from_edge_list(edges, node_weights);
  result.mapping = std::move(mapping);
  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/graph_contraction.cpp` to `PreprocessingSources` in `mt-kahypar/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=GraphContractionTest.*`
Expected: PASS (3 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/graph_contraction.h \
        mt-kahypar/partition/preprocessing/filtering/graph_contraction.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/graph_contraction_test.cc
git commit -m "feat: add shared contract_graph primitive for the filtering module

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Parallel connectivity + spanning forest

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/parallel_connectivity_test.cc`

**Interfaces:**
- Consumes: `FilterGraph`, `AtomicUnionFind` (Task 3).
- Produces:
  - `parallel_connected_components(graph) -> std::vector<NodeID>` (component id per vertex; ids are union-find representatives, i.e. small vertex ids, not a dense range).
  - `parallel_connected_components_excluding(graph, excluded_edge /* size numEdges(), true = excluded */) -> std::vector<NodeID>` (same, but treats excluded edges as absent).
  - `SpanningForest{parent, parent_edge, bfs_order, roots}`.
  - `build_spanning_forest(graph, component, roots) -> SpanningForest`. If `roots` is empty, one root per distinct value in `component` is chosen automatically (its lowest-numbered vertex); otherwise `roots` must contain exactly one vertex per distinct component value, and those are used as BFS roots.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/parallel_connectivity_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

using namespace mt_kahypar::filtering;

namespace {
FilterGraph make_two_triangles() {
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(ParallelConnectivityTest, FindsTwoComponents) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  EXPECT_EQ(component[0], component[1]);
  EXPECT_EQ(component[1], component[2]);
  EXPECT_EQ(component[3], component[4]);
  EXPECT_EQ(component[4], component[5]);
  EXPECT_NE(component[0], component[3]);
}

TEST(ParallelConnectivityTest, ExcludingAnEdgeCanSplitAComponent) {
  // A 4-cycle 0-1-2-3-0 is 2-edge-connected; excluding one edge still leaves
  // it connected (it's a cycle), but excluding two adjacent edges splits it.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::vector<char> excluded(4, 0);
  excluded[0] = 1;  // edge (0,1)
  excluded[1] = 1;  // edge (1,2)
  std::vector<NodeID> component = parallel_connected_components_excluding(graph, excluded);
  EXPECT_NE(component[1], component[3]);  // vertex 1 is now isolated from 3
}

TEST(ParallelConnectivityTest, SpanningForestCoversAllVerticesWithDefaultRoots) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});

  ASSERT_EQ(forest.roots.size(), 2u);
  ASSERT_EQ(forest.bfs_order.size(), 6u);
  for (NodeID v = 0; v < 6; ++v) EXPECT_NE(forest.parent[v], kInvalidNode);
  for (NodeID root : forest.roots) {
    EXPECT_EQ(forest.parent[root], root);
    EXPECT_EQ(forest.parent_edge[root], kInvalidEdge);
  }
  for (NodeID v = 0; v < 6; ++v) {
    if (forest.parent[v] != v) EXPECT_NE(forest.parent_edge[v], kInvalidEdge);
  }
}

TEST(ParallelConnectivityTest, SpanningForestHonorsExplicitRoots) {
  FilterGraph graph = make_two_triangles();
  std::vector<NodeID> component = parallel_connected_components(graph);
  // Force roots 2 and 5 explicitly (one per component, each a valid member
  // of that component in make_two_triangles' fixed layout: {0,1,2} and
  // {3,4,5}).
  std::vector<NodeID> roots = {2, 5};
  SpanningForest forest = build_spanning_forest(graph, component, roots);
  EXPECT_EQ(forest.parent[2], 2u);
  EXPECT_EQ(forest.parent[5], 5u);
  // The explicit-roots path must be as fully covered as the default-roots
  // path: every vertex reached, every non-root's tree edge in its own
  // component, and each vertex's parent chain leading back to the root
  // that owns its component (not the other root).
  ASSERT_EQ(forest.bfs_order.size(), 6u);
  for (NodeID v = 0; v < 6; ++v) {
    EXPECT_NE(forest.parent[v], kInvalidNode);
    const NodeID owning_root = (component[v] == component[2]) ? NodeID(2) : NodeID(5);
    NodeID cur = v;
    while (forest.parent[cur] != cur) cur = forest.parent[cur];
    EXPECT_EQ(cur, owning_root);
    if (v != 2 && v != 5) EXPECT_NE(forest.parent_edge[v], kInvalidEdge);
  }
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/parallel_connectivity_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `parallel_connectivity.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
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
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.cpp
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

#include <cassert>
#include <deque>
#include <unordered_map>

#include <tbb/parallel_for.h>

#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

namespace {
std::vector<NodeID> components_from_union_find(const FilterGraph& graph, AtomicUnionFind& uf) {
  const size_t n = graph.numNodes();
  std::vector<NodeID> component(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    component[v] = uf.find(static_cast<uint32_t>(v));
  });
  return component;
}
}  // namespace

std::vector<NodeID> parallel_connected_components(const FilterGraph& graph) {
  const size_t n = graph.numNodes();
  AtomicUnionFind uf(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      uf.unite(static_cast<uint32_t>(v), static_cast<uint32_t>(graph.adj[pos]));
    }
  });
  return components_from_union_find(graph, uf);
}

std::vector<NodeID> parallel_connected_components_excluding(
    const FilterGraph& graph, const std::vector<char>& excluded_edge) {
  const size_t n = graph.numNodes();
  AtomicUnionFind uf(n);
  tbb::parallel_for(size_t(0), n, [&](size_t v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      if (excluded_edge[graph.adj_edge[pos]]) continue;
      uf.unite(static_cast<uint32_t>(v), static_cast<uint32_t>(graph.adj[pos]));
    }
  });
  return components_from_union_find(graph, uf);
}

SpanningForest build_spanning_forest(const FilterGraph& graph,
                                      const std::vector<NodeID>& component,
                                      const std::vector<NodeID>& roots_in) {
  const size_t n = graph.numNodes();

  std::vector<NodeID> roots = roots_in;
  if (roots.empty()) {
    std::unordered_map<NodeID, bool> seen;
    for (size_t v = 0; v < n; ++v) {
      if (seen.emplace(component[v], true).second) roots.push_back(static_cast<NodeID>(v));
    }
  } else {
    // Debug-only precondition check: no two supplied roots may belong to
    // the same real component -- that reproduces the same class of race as
    // passing an excluding-variant component here (see the header comment).
    assert([&] {
      std::vector<char> seen_component(n, 0);
      for (NodeID root : roots_in) {
        if (root >= n) return false;
        if (seen_component[component[root]]) return false;
        seen_component[component[root]] = 1;
      }
      return true;
    }() && "roots must contain at most one entry per distinct component value, each a valid vertex id");
  }

  SpanningForest forest;
  forest.parent.assign(n, kInvalidNode);
  forest.parent_edge.assign(n, kInvalidEdge);
  forest.roots = roots;
  std::vector<std::vector<NodeID>> order_per_root(roots.size());

  // Each root's true BFS-reachable set in `graph` is exactly its own
  // component (guaranteed only when `component` came from
  // parallel_connected_components on this same graph -- see the header
  // comment), so concurrently-running BFS instances never touch the same
  // vertex, hence never the same index of forest.parent/forest.parent_edge
  // -- no data race despite plain (non-atomic) vectors.
  tbb::parallel_for(size_t(0), roots.size(), [&](size_t i) {
    const NodeID root = roots[i];
    std::vector<NodeID>& order = order_per_root[i];
    std::deque<NodeID> queue;
    forest.parent[root] = root;
    queue.push_back(root);
    order.push_back(root);
    while (!queue.empty()) {
      const NodeID u = queue.front();
      queue.pop_front();
      for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
        const NodeID v = graph.adj[pos];
        if (forest.parent[v] == kInvalidNode) {
          forest.parent[v] = u;
          forest.parent_edge[v] = graph.adj_edge[pos];
          queue.push_back(v);
          order.push_back(v);
        }
      }
    }
  });

  forest.bfs_order.reserve(n);
  for (const std::vector<NodeID>& order : order_per_root) {
    forest.bfs_order.insert(forest.bfs_order.end(), order.begin(), order.end());
  }
  return forest;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/parallel_connectivity.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=ParallelConnectivityTest.*`
Expected: PASS (4 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h \
        mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/parallel_connectivity_test.cc
git commit -m "feat: add parallel connectivity and spanning forest primitives

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

**Checkpoint after Task 5:** all shared primitives (graph, IO, union-find, contraction, connectivity) are in place and tested. Tasks 6-11 build Part 1 (tiny-cut detection) on top of them.

---

### Task 6: Cut signatures (Pritchard-Thurimella core) — bridges

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/cut_signatures.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/cut_signatures.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/cut_signatures_test.cc`

This is the shared primitive behind both Part 1's bridge detection and Part 1's 2-edge-cut detection (Task 7): a spanning tree plus random-label XOR aggregation (Pritchard & Thurimella, ACM TALG 7(4), 2011 — design spec section 4.4/4.5, correctness verified there against brute force on 300 random graphs plus the cycle-graph edge case).

**Interfaces:**
- Consumes: `FilterGraph`, `SpanningForest` (Task 5).
- Produces:
  - `EdgeSignatures{signature (vector<unsigned __int128>, size numEdges()), is_tree_edge (vector<char>, size numEdges())}`.
  - `compute_edge_signatures(graph, forest) -> EdgeSignatures`. For a tree edge, `signature[e]` is the XOR of random labels of every non-tree edge whose fundamental cycle covers it (0 if none — i.e. a bridge). For a non-tree edge, `signature[e]` is its own random label.
  - `compute_bridges(graph, sigs) -> std::vector<char>` (size numEdges(), true for bridges). A bridge is a tree edge with `signature == 0`; non-tree edges are never bridges.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/cut_signatures_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

using namespace mt_kahypar::filtering;

namespace {
EdgeSignatures signatures_for(const FilterGraph& graph) {
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  return compute_edge_signatures(graph, forest);
}
}  // namespace

TEST(CutSignaturesTest, EveryEdgeOfAPathIsABridge) {
  // 0-1-2-3, a path: every edge is a bridge.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  for (EdgeID e = 0; e < 3; ++e) EXPECT_TRUE(is_bridge[e]);
}

TEST(CutSignaturesTest, NoEdgeOfACycleIsABridge) {
  // A simple cycle 0-1-2-3-0: no bridges (removing any single edge leaves it
  // connected via the rest of the cycle).
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  for (EdgeID e = 0; e < 4; ++e) EXPECT_FALSE(is_bridge[e]);
}

TEST(CutSignaturesTest, BridgeConnectingTwoCyclesIsFound) {
  // Two triangles (0,1,2) and (3,4,5) joined by a single bridge edge (2,3).
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},  // the bridge
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  EdgeSignatures sigs = signatures_for(graph);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);
  int bridge_count = 0;
  EdgeID bridge_id = kInvalidEdge;
  for (EdgeID e = 0; e < 7; ++e) {
    if (is_bridge[e]) { ++bridge_count; bridge_id = e; }
  }
  EXPECT_EQ(bridge_count, 1);
  auto endpoints = compute_edge_endpoints(graph);
  EXPECT_TRUE((endpoints[bridge_id] == std::make_pair(NodeID(2), NodeID(3))) ||
              (endpoints[bridge_id] == std::make_pair(NodeID(3), NodeID(2))));
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/cut_signatures_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `cut_signatures.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/cut_signatures.h
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

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/cut_signatures.cpp
#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"

#include <random>

namespace mt_kahypar {
namespace filtering {

namespace {
unsigned __int128 random_128() {
  static thread_local std::mt19937_64 rng(std::random_device{}());
  const uint64_t hi = rng();
  const uint64_t lo = rng();
  return (static_cast<unsigned __int128>(hi) << 64) | lo;
}
}  // namespace

EdgeSignatures compute_edge_signatures(const FilterGraph& graph, const SpanningForest& forest) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();

  EdgeSignatures sigs;
  sigs.signature.assign(m, 0);
  sigs.is_tree_edge.assign(m, 0);

  for (size_t v = 0; v < n; ++v) {
    if (forest.parent_edge[v] != kInvalidEdge) sigs.is_tree_edge[forest.parent_edge[v]] = 1;
  }

  // diff[v] will hold the XOR of labels of non-tree edges whose fundamental
  // cycle's tree-path endpoint is exactly v (i.e. one endpoint of the
  // non-tree edge); accumulating diff bottom-up along the tree turns this
  // into, for each tree edge (parent(v), v), the XOR over all non-tree edges
  // whose path covers that edge.
  std::vector<unsigned __int128> diff(n, 0);

  for (size_t v = 0; v < n; ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const EdgeID e = graph.adj_edge[pos];
      const NodeID other = graph.adj[pos];
      if (sigs.is_tree_edge[e]) continue;
      if (other <= static_cast<NodeID>(v)) continue;  // process each non-tree edge once
      const unsigned __int128 label = random_128();
      sigs.signature[e] = label;  // non-tree edge's own signature is its label
      diff[v] ^= label;
      diff[other] ^= label;
    }
  }

  // Accumulate diff into per-tree-edge signatures via reverse BFS order:
  // BFS order is non-decreasing in depth, so processing in reverse guarantees
  // every child is folded into its parent before the parent itself is used.
  for (auto it = forest.bfs_order.rbegin(); it != forest.bfs_order.rend(); ++it) {
    const NodeID v = *it;
    const NodeID parent = forest.parent[v];
    if (parent == v) continue;  // root
    sigs.signature[forest.parent_edge[v]] = diff[v];
    diff[parent] ^= diff[v];
  }

  return sigs;
}

std::vector<char> compute_bridges(const FilterGraph& graph, const EdgeSignatures& sigs) {
  std::vector<char> is_bridge(graph.numEdges(), 0);
  for (EdgeID e = 0; e < static_cast<EdgeID>(graph.numEdges()); ++e) {
    is_bridge[e] = sigs.is_tree_edge[e] && sigs.signature[e] == 0;
  }
  return is_bridge;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/cut_signatures.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=CutSignaturesTest.*`
Expected: PASS (3 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/cut_signatures.h \
        mt-kahypar/partition/preprocessing/filtering/cut_signatures.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/cut_signatures_test.cc
git commit -m "feat: add Pritchard-Thurimella edge signatures and bridge detection

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 7: 2-edge-cut equivalence classes

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/cut_signatures.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/cut_signatures.cpp`
- Test: `tests/partition/preprocessing/filtering/cut_signatures_test.cc` (append)

**Interfaces:**
- Consumes: `EdgeSignatures`, `compute_edge_endpoints`, `parallel_connected_components_excluding` (for the local disconnection check).
- Produces: `find_two_edge_cut_classes(graph, sigs) -> std::vector<std::vector<EdgeID>>` — each inner vector is a verified 2-edge-cut equivalence class (size >= 2), excluding bridges. Matches PUNCH's equivalence relation `P` (spec section 4.4).

- [ ] **Step 1: Write the failing tests (appended to `cut_signatures_test.cc`)**

```cpp
namespace {
// Sorts each class and the outer list so equality comparisons are order-independent.
std::vector<std::vector<EdgeID>> normalize(std::vector<std::vector<EdgeID>> classes) {
  for (auto& c : classes) std::sort(c.begin(), c.end());
  std::sort(classes.begin(), classes.end());
  return classes;
}
}  // namespace

TEST(CutSignaturesTest, CycleGraphIsOneClassOfAllEdges) {
  // A 6-cycle: removing any two edges disconnects it, so all 6 edges form a
  // single equivalence class (this is the "Case B" scenario -- a class
  // containing a non-tree edge -- that a naive tree-only bucketing scheme
  // would miss; see design spec section 4.5).
  std::vector<NodeWeight> weights(6, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 5, 1}, {5, 0, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(graph, sigs);
  ASSERT_EQ(classes.size(), 1u);
  EXPECT_EQ(classes[0].size(), 6u);
}

TEST(CutSignaturesTest, ThetaGraphHasTwoIndependentCutPairs) {
  // A "theta graph": nodes 0 and 3 connected by three internally-disjoint
  // paths: 0-1-3, 0-2-3, 0-4-3. Removing both edges of any one path
  // disconnects that path's middle vertex, but does not disconnect the whole
  // graph -- the two edges *incident to a degree-2 middle vertex* (e.g. (0,1)
  // and (1,3)) form a 2-cut isolating {1}. The three middle vertices give
  // three independent 2-cut classes of size 2 each; no larger cross-path
  // pairing is a valid cut.
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 3, 1},
    {0, 2, 1}, {2, 3, 1},
    {0, 4, 1}, {4, 3, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(graph, sigs);
  ASSERT_EQ(classes.size(), 3u);
  for (auto& c : classes) EXPECT_EQ(c.size(), 2u);
}

TEST(CutSignaturesTest, BridgesAreExcludedFromCutClasses) {
  // Two triangles joined by a bridge: the bridge must not appear in any
  // 2-cut class (it's a 1-cut on its own), and the two triangles' edges are
  // 3-edge-connected internally, so no 2-cut classes exist at all.
  std::vector<NodeWeight> weights(6, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},
    {3, 4, 1}, {4, 5, 1}, {3, 5, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);
  EdgeSignatures sigs = signatures_for(graph);

  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(graph, sigs);
  EXPECT_TRUE(classes.empty());
}
```

Add `#include <algorithm>` to the top of `cut_signatures_test.cc`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `find_two_edge_cut_classes` is not declared

- [ ] **Step 3: Add the declaration**

Append to `cut_signatures.h` (inside `namespace filtering`):

```cpp
// Finds all 2-edge-cut equivalence classes of `graph` (design spec section
// 4.4): buckets every non-bridge edge by its signature (tree edges by their
// aggregated label, non-tree edges by their own label -- see
// compute_edge_signatures), then verifies each candidate class of size >= 2
// by checking that its first two edges' removal actually disconnects the
// graph locally. Verified classes are returned; unverified (collision)
// candidates are dropped (Monte Carlo false positives are astronomically
// rare at 128 bits, but the check is cheap enough to always do -- see design
// spec section 4.4 step 7).
std::vector<std::vector<EdgeID>> find_two_edge_cut_classes(const FilterGraph& graph,
                                                            const EdgeSignatures& sigs);
```

- [ ] **Step 4: Write the implementation**

Append to `cut_signatures.cpp`:

```cpp
#include <unordered_map>

#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

namespace {
struct SignatureHash {
  size_t operator()(unsigned __int128 x) const {
    const uint64_t hi = static_cast<uint64_t>(x >> 64);
    const uint64_t lo = static_cast<uint64_t>(x);
    return std::hash<uint64_t>()(hi) ^ (std::hash<uint64_t>()(lo) * 0x9e3779b97f4a7c15ULL);
  }
};

// A class is verified if removing all its edges actually disconnects at
// least one pair of its incident vertices from each other. Cheap at road-
// network scale since classes are small (size 2 in the overwhelming common
// case) and this check runs once per candidate class, not per edge.
bool verify_class_disconnects(const FilterGraph& graph, const std::vector<EdgeID>& cls) {
  std::vector<char> excluded(graph.numEdges(), 0);
  for (EdgeID e : cls) excluded[e] = 1;
  std::vector<NodeID> component = parallel_connected_components_excluding(graph, excluded);
  auto endpoints = compute_edge_endpoints(graph);
  const NodeID u = endpoints[cls[0]].first;
  const NodeID v = endpoints[cls[0]].second;
  return component[u] != component[v];
}
}  // namespace

std::vector<std::vector<EdgeID>> find_two_edge_cut_classes(const FilterGraph& graph,
                                                            const EdgeSignatures& sigs) {
  std::unordered_map<unsigned __int128, std::vector<EdgeID>, SignatureHash> buckets;
  for (EdgeID e = 0; e < static_cast<EdgeID>(graph.numEdges()); ++e) {
    const bool is_bridge = sigs.is_tree_edge[e] && sigs.signature[e] == 0;
    if (is_bridge) continue;
    buckets[sigs.signature[e]].push_back(e);
  }

  std::vector<std::vector<EdgeID>> classes;
  for (auto& [signature, edges] : buckets) {
    if (edges.size() < 2) continue;
    if (verify_class_disconnects(graph, edges)) classes.push_back(std::move(edges));
  }
  return classes;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=CutSignaturesTest.*`
Expected: PASS (6 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/cut_signatures.h \
        mt-kahypar/partition/preprocessing/filtering/cut_signatures.cpp \
        tests/partition/preprocessing/filtering/cut_signatures_test.cc
git commit -m "feat: find 2-edge-cut equivalence classes via signature bucketing

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 8: Pass 1 — component (bridge) tree contraction

This is the most involved task in Part 1: build the bridge tree (nodes = 2-edge-connected blocks, edges = bridges — see the "Correction from the spec" note in Global Constraints), root it at its heaviest block, and contract any subtree with total weight <= U top-down, with the tau-merge-into-parent extension.

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc`

**Interfaces:**
- Consumes: `FilterGraph`, `ContractionResult`, `contract_graph` (Task 4); `parallel_connected_components`, `parallel_connected_components_excluding`, `build_spanning_forest`, `SpanningForest` (Task 5); `compute_edge_signatures`, `compute_bridges` (Task 6); `compute_edge_endpoints`, `build_csr_from_edge_list` (Task 1); `AtomicUnionFind` (Task 3).
- Produces: `TinyCutParams{U, tau=5}`, `contract_component_tree(graph, params) -> ContractionResult`.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

using namespace mt_kahypar::filtering;

TEST(ContractComponentTreeTest, ContractsASmallLeafHangingOffABridge) {
  // A big triangle (0,1,2) with a small pendant vertex 3 attached via a
  // bridge (2,3), and an even smaller pendant 4 attached to 3 via another
  // bridge. Weights: main triangle vertices weight 100 each (heavy = root),
  // vertices 3 and 4 weight 1 each. With U = 5, the whole {3,4} subtree
  // (total weight 2) should contract into a single vertex, still attached to
  // the triangle via the bridge (2,3) -- so the result has 4 nodes (the
  // triangle's 3 untouched vertices, plus one contracted {3,4} vertex) and
  // the bridge survives as an edge.
  std::vector<NodeWeight> weights = {100, 100, 100, 1, 1};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1},
    {3, 4, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{5, 5});
  EXPECT_EQ(result.graph.numNodes(), 4u);
  EXPECT_EQ(result.mapping[3], result.mapping[4]);
  EXPECT_NE(result.mapping[2], result.mapping[3]);
  const NodeWeight contracted_weight = result.graph.node_weight[result.mapping[3]];
  EXPECT_EQ(contracted_weight, 2u);
}

TEST(ContractComponentTreeTest, TauMergeFoldsSmallSubtreeIntoParentSpecifically) {
  // A 3-block chain: root(weight 1000, its subtree alone is far too big to
  // ever be chosen as a whole) -- mid(weight 2) -- leaf(weight 1). With
  // U = 10, tau = 5: mid's own subtree {mid, leaf} has weight 3 <= U, so it
  // gets chosen and contracted as a unit; since 3 <= tau AND 3 + root's own
  // weight... wait, root's own weight (1000) alone already exceeds U, so no
  // merge into root should happen here either. This deliberately isolates
  // the "chosen subtree, but tau-merge condition fails because the parent is
  // heavy" case, distinct from the plain "whole graph fits" case.
  std::vector<NodeWeight> weights = {1000, 2, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.mapping[1], result.mapping[2]);   // mid+leaf contract together
  EXPECT_NE(result.mapping[0], result.mapping[1]);   // but root stays separate
  EXPECT_EQ(result.graph.numNodes(), 2u);
}

TEST(ContractComponentTreeTest, TauMergeActuallyFusesIntoALightParent) {
  // root(weight 2, light enough itself) -- leaf(weight 2). Root's own
  // subtree (root+leaf, weight 4) already fits under U = 10 directly, so
  // this alone contracts everything without needing tau-merge -- included
  // as a baseline. The interesting case is TauMergePreventsCascadingOverflow
  // below, which is the one that actually isolates the tau-merge-into-a-
  // not-otherwise-chosen-parent path.
  std::vector<NodeWeight> weights = {2, 2};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
}

TEST(ContractComponentTreeTest, TauMergePreventsCascadingOverflow) {
  // R(weight 1000, unambiguously the heaviest block so it's picked as root)
  // -- P(weight 1, a light non-root parent) -- {child_a(weight 4),
  // child_b(weight 4)} (P's two children). With U = 6, tau = 5: P's own
  // subtree (1+4+4=9) is too big to be chosen as a whole (9 > 6), so P
  // itself is never marked chosen. child_a's subtree (weight 4) <= U and
  // <= tau, and 4 + P's own weight (1) = 5 <= U, so child_a tau-merges into
  // P. child_b's subtree also individually satisfies 4 + 1 = 5 <= U -- but
  // P has ALREADY absorbed child_a's 4 units of extra weight, so the TRUE
  // combined result of also merging child_b (1 + 4 + 4 = 9) would exceed U.
  // The implementation must track this and skip child_b's tau-merge, leaving
  // child_b's own already-contracted subtree (just itself) standing on its
  // own rather than fusing it into the now-full P group. This is the exact
  // bug caught during this plan's pre-flight review (see the Global
  // Constraints note on this task) -- if the fix regresses, this test will
  // fail by observing a merged group heavier than U.
  std::vector<NodeWeight> weights = {1000, 1, 4, 4};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {1, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{6, 5});

  for (size_t v = 0; v < result.graph.numNodes(); ++v) {
    EXPECT_LE(result.graph.node_weight[v], 6u) << "vertex " << v << " exceeds U";
  }
  EXPECT_NE(result.mapping[0], result.mapping[1]);  // R never merges with anything
  // Exactly one of P's two children fused with P; the other stands alone.
  const bool a_fused_with_p = (result.mapping[1] == result.mapping[2]);
  const bool b_fused_with_p = (result.mapping[1] == result.mapping[3]);
  EXPECT_TRUE(a_fused_with_p != b_fused_with_p);
}

TEST(ContractComponentTreeTest, LargeGraphIsLeftUntouchedWhenNoSubtreeFits) {
  // A single triangle where every vertex is heavy: no subtree (other than
  // the whole graph, which isn't a proper subtree hanging off a bridge --
  // there are no bridges at all here) can be <= U, so nothing contracts.
  std::vector<NodeWeight> weights = {100, 100, 100};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_component_tree(graph, TinyCutParams{10, 5});
  EXPECT_EQ(result.graph.numNodes(), 3u);
  EXPECT_EQ(result.graph.numEdges(), 3u);
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/tiny_cut_detection_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `tiny_cut_detection.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"

namespace mt_kahypar {
namespace filtering {

struct TinyCutParams {
  NodeWeight U;
  NodeWeight tau = 5;
};

// Part 1, pass 1 (design spec section 4.2, corrected per the "Correction
// from the spec" note in this plan's Global Constraints): builds the bridge
// tree (nodes = 2-edge-connected blocks, edges = bridges), roots it at its
// heaviest block, and contracts any subtree with total weight <= U top-down.
// A contracted subtree of weight <= tau is additionally folded into its
// parent block if the merged weight is still <= U.
ContractionResult contract_component_tree(const FilterGraph& graph, const TinyCutParams& params);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

#include <unordered_map>

#include "mt-kahypar/partition/preprocessing/filtering/cut_signatures.h"
#include "mt-kahypar/partition/preprocessing/filtering/parallel_connectivity.h"

namespace mt_kahypar {
namespace filtering {

ContractionResult contract_component_tree(const FilterGraph& graph, const TinyCutParams& params) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();

  // Step 1: bridges, via the shared spanning-tree XOR-signature technique.
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  EdgeSignatures sigs = compute_edge_signatures(graph, forest);
  std::vector<char> is_bridge = compute_bridges(graph, sigs);

  // Step 2: blocks = connected components of the graph with bridges removed.
  std::vector<NodeID> block_rep = parallel_connected_components_excluding(graph, is_bridge);
  std::unordered_map<NodeID, NodeID> dense_block_id;
  std::vector<NodeID> block_of(n);
  for (size_t v = 0; v < n; ++v) {
    auto it = dense_block_id.find(block_rep[v]);
    if (it == dense_block_id.end()) {
      const NodeID id = static_cast<NodeID>(dense_block_id.size());
      dense_block_id.emplace(block_rep[v], id);
      block_of[v] = id;
    } else {
      block_of[v] = it->second;
    }
  }
  const NodeID num_blocks = static_cast<NodeID>(dense_block_id.size());

  std::vector<NodeWeight> block_weight(num_blocks, 0);
  for (size_t v = 0; v < n; ++v) block_weight[block_of[v]] += graph.node_weight[v];

  // Step 3: the block quotient graph. Its edges are exactly the bridges
  // (each bridge connects exactly two distinct blocks); this quotient graph
  // is a forest -- one tree per connected component of the original graph --
  // matching the design spec's "component tree T".
  std::vector<std::pair<NodeID, NodeID>> edge_endpoints = compute_edge_endpoints(graph);
  std::vector<EdgeListEntry> quotient_edges;
  for (EdgeID e = 0; e < static_cast<EdgeID>(m); ++e) {
    if (!is_bridge[e]) continue;
    quotient_edges.push_back(EdgeListEntry{
        block_of[edge_endpoints[e].first], block_of[edge_endpoints[e].second], graph.edge_weight[e]});
  }
  FilterGraph quotient = build_csr_from_edge_list(quotient_edges, block_weight);

  // Step 4: root each quotient-forest component at its heaviest block
  // ("the edge-connected component with maximum size", design spec 4.2).
  std::vector<NodeID> quotient_component = parallel_connected_components(quotient);
  std::unordered_map<NodeID, NodeID> best_root_for_component;
  std::unordered_map<NodeID, NodeWeight> best_weight_for_component;
  for (NodeID b = 0; b < num_blocks; ++b) {
    const NodeID comp = quotient_component[b];
    auto it = best_weight_for_component.find(comp);
    if (it == best_weight_for_component.end() || quotient.node_weight[b] > it->second) {
      best_weight_for_component[comp] = quotient.node_weight[b];
      best_root_for_component[comp] = b;
    }
  }
  std::vector<NodeID> roots;
  roots.reserve(best_root_for_component.size());
  for (auto& [comp, root] : best_root_for_component) roots.push_back(root);
  SpanningForest quotient_forest = build_spanning_forest(quotient, quotient_component, roots);

  // Step 5: subtree weights over the quotient tree, bottom-up via reverse
  // BFS order (BFS order is non-decreasing in depth, so every child is
  // folded into its parent before the parent itself is used).
  std::vector<NodeWeight> subtree_weight(num_blocks);
  for (NodeID b = 0; b < num_blocks; ++b) subtree_weight[b] = quotient.node_weight[b];
  for (auto it = quotient_forest.bfs_order.rbegin(); it != quotient_forest.bfs_order.rend(); ++it) {
    const NodeID b = *it;
    const NodeID parent = quotient_forest.parent[b];
    if (parent != b) subtree_weight[parent] += subtree_weight[b];
  }

  // Step 6: one representative original vertex per block, computed
  // unconditionally (a tau-merge target block may never itself be "chosen"
  // -- see below -- so every block needs a representative available).
  std::vector<NodeID> block_representative_vertex(num_blocks, kInvalidNode);
  for (size_t v = 0; v < n; ++v) {
    NodeID& rep = block_representative_vertex[block_of[v]];
    if (rep == kInvalidNode) rep = static_cast<NodeID>(v);
  }

  // Step 7: top-down selection of subtrees to contract (weight <= U), plus
  // the tau-merge into the parent block, building the vertex-level
  // union-find DIRECTLY (not via an intermediate block-level union-find --
  // see the note below on why that translation is unsound).
  //
  // `is_within_chosen_subtree[b]` marks a block whose own subtree_weight
  // triggered contraction, either directly or by inheriting from an
  // already-chosen ancestor (forward BFS order visits parents before
  // children, so this propagates correctly downward). A tau-merge target
  // (the immediate parent of a chosen small subtree) is unioned into the
  // vertex group directly, WITHOUT setting is_within_chosen_subtree on the
  // parent itself: doing so would incorrectly make every OTHER child of that
  // parent auto-absorb regardless of its own size, which is not what a
  // tau-merge means (it is a narrow, single-subtree-into-its-immediate-
  // parent extension, not "the parent is now fully contracted").
  //
  // `tau_merged_extra_weight[parent]` guards against a genuine correctness
  // bug that a naive per-merge check misses: TWO SEPARATE small subtrees
  // hanging off the SAME parent can each individually satisfy
  // "subtree + parent's own weight <= U", yet their COMBINED weight once
  // both are merged into the same parent can exceed U (e.g. parent weight 1,
  // two children of weight 4 each, U = 6: 4+1 <= 6 passes twice, but
  // 4+4+1 = 9 > 6). Tracking how much extra weight has already been folded
  // into a given parent via prior tau-merges, and including it in the check,
  // prevents this -- since Part 1's own U-cap is part of what keeps Part 2's
  // final fragments within U (see design spec section 6), not merely Part
  // 2's alpha <= 1 guarantee on its own.
  std::vector<char> is_within_chosen_subtree(num_blocks, 0);
  std::vector<NodeWeight> tau_merged_extra_weight(num_blocks, 0);
  AtomicUnionFind vertex_groups(n);
  for (NodeID b : quotient_forest.bfs_order) {
    const NodeID parent = quotient_forest.parent[b];
    if (parent != b && is_within_chosen_subtree[parent]) {
      is_within_chosen_subtree[b] = 1;
      vertex_groups.unite(block_representative_vertex[b], block_representative_vertex[parent]);
      continue;
    }
    if (subtree_weight[b] <= params.U) {
      is_within_chosen_subtree[b] = 1;
      if (parent != b && subtree_weight[b] <= params.tau &&
          subtree_weight[b] + quotient.node_weight[parent] + tau_merged_extra_weight[parent] <= params.U) {
        vertex_groups.unite(block_representative_vertex[b], block_representative_vertex[parent]);
        tau_merged_extra_weight[parent] += subtree_weight[b];
      }
    }
  }

  // Step 8: collapse every vertex of a chosen block into that block's own
  // representative (blocks that were never chosen, and never tau-merge
  // targets, keep their original internal vertices and edges untouched --
  // this is intentional: a block too big to contract must retain its
  // internal structure for pass 2/3 to potentially reduce further).
  for (size_t v = 0; v < n; ++v) {
    const NodeID b = block_of[v];
    if (is_within_chosen_subtree[b]) {
      vertex_groups.unite(block_representative_vertex[b], static_cast<NodeID>(v));
    }
  }

  return contract_graph(graph, vertex_groups);
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/tiny_cut_detection.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=ContractComponentTreeTest.*`
Expected: PASS (5 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc
git commit -m "feat: contract the bridge tree top-down (Part 1 pass 1)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 9: Pass 2 — degree-2 chain contraction

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp`
- Test: `tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc` (append)

**Interfaces:**
- Produces: `contract_degree2_chains(graph, U) -> ContractionResult`. Whole-chain-or-nothing: a maximal chain of degree-2 vertices contracts to one vertex only if its total weight is `<= U`; otherwise it is left untouched (no partial contraction). A connected component consisting entirely of degree-2 vertices (a pure cycle with no anchor) is treated as one chain.

- [ ] **Step 1: Write the failing tests (appended to `tiny_cut_detection_test.cc`)**

```cpp
TEST(ContractDegree2ChainsTest, ContractsAPathBetweenTwoAnchors) {
  // A(deg1) - B(deg2) - C(deg2) - D(deg1). B and C should merge.
  std::vector<NodeWeight> weights = {1, 1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 3u);
  EXPECT_EQ(result.mapping[1], result.mapping[2]);
  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_NE(result.mapping[2], result.mapping[3]);
}

TEST(ContractDegree2ChainsTest, LeavesChainUntouchedWhenTooBig) {
  std::vector<NodeWeight> weights = {1, 100, 100, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 4u);
  EXPECT_NE(result.mapping[1], result.mapping[2]);
}

TEST(ContractDegree2ChainsTest, ContractsAPureCycleWithNoAnchor) {
  // A 5-cycle: every vertex has degree 2, so there's no anchor. The whole
  // component is one chain.
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = contract_degree2_chains(graph, 10);
  EXPECT_EQ(result.graph.numNodes(), 1u);
  EXPECT_EQ(result.graph.numEdges(), 0u);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `contract_degree2_chains` is not declared

- [ ] **Step 3: Add the declaration**

Append to `tiny_cut_detection.h` (inside `namespace filtering`):

```cpp
// Part 1, pass 2 (design spec section 4.3): contracts each maximal chain of
// degree-2 vertices into a single vertex, provided the chain's total weight
// is <= U (whole-chain-or-nothing; a chain that's too big is left as-is). A
// connected component made entirely of degree-2 vertices (a pure cycle, with
// no anchor of different degree) is treated as one chain.
ContractionResult contract_degree2_chains(const FilterGraph& graph, NodeWeight U);
```

- [ ] **Step 4: Write the implementation**

Append to `tiny_cut_detection.cpp`:

```cpp
namespace {
std::vector<NodeID> walk_chain_from(const FilterGraph& graph, const std::vector<char>& is_deg2,
                                     std::vector<char>& visited, NodeID prev, NodeID cur) {
  std::vector<NodeID> chain;
  while (is_deg2[cur] && !visited[cur]) {
    visited[cur] = 1;
    chain.push_back(cur);
    NodeID next = kInvalidNode;
    for (EdgeID pos = graph.node_begin[cur]; pos < graph.node_begin[cur + 1]; ++pos) {
      if (graph.adj[pos] != prev) { next = graph.adj[pos]; break; }
    }
    prev = cur;
    cur = next;
  }
  return chain;
}
}  // namespace

ContractionResult contract_degree2_chains(const FilterGraph& graph, NodeWeight U) {
  const size_t n = graph.numNodes();
  std::vector<char> is_deg2(n);
  for (size_t v = 0; v < n; ++v) is_deg2[v] = (graph.degree(static_cast<NodeID>(v)) == 2);

  std::vector<char> visited(n, 0);
  AtomicUnionFind uf(n);

  auto try_contract_chain = [&](const std::vector<NodeID>& chain) {
    if (chain.empty()) return;
    NodeWeight total = 0;
    for (NodeID v : chain) total += graph.node_weight[v];
    if (total <= U) {
      for (size_t i = 1; i < chain.size(); ++i) uf.unite(chain[0], chain[i]);
    }
  };

  // Chains anchored at both ends by a vertex of degree != 2.
  for (NodeID anchor = 0; anchor < static_cast<NodeID>(n); ++anchor) {
    if (is_deg2[anchor]) continue;
    for (EdgeID pos = graph.node_begin[anchor]; pos < graph.node_begin[anchor + 1]; ++pos) {
      const NodeID cur = graph.adj[pos];
      if (!is_deg2[cur] || visited[cur]) continue;
      try_contract_chain(walk_chain_from(graph, is_deg2, visited, anchor, cur));
    }
  }

  // Pure cycles: components with no anchor at all.
  for (NodeID v = 0; v < static_cast<NodeID>(n); ++v) {
    if (!is_deg2[v] || visited[v]) continue;
    try_contract_chain(walk_chain_from(graph, is_deg2, visited, kInvalidNode, v));
  }

  return contract_graph(graph, uf);
}
```

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=ContractDegree2ChainsTest.*`
Expected: PASS (3 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp \
        tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc
git commit -m "feat: contract degree-2 chains (Part 1 pass 2)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 10: Pass 3 processing — contract 2-edge-cut classes

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp`
- Test: `tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc` (append)

**Interfaces:**
- Consumes: `find_two_edge_cut_classes` (Task 7).
- Produces: `contract_two_edge_cuts(graph, U) -> ContractionResult`. For each equivalence class `S`, computes connected components of `(V, E \ S)` and unions every component with total weight `<= U` (per-class recomputation, no two-at-a-time traversal optimization — see the "Deferred from the spec" note in Global Constraints).

- [ ] **Step 1: Write the failing tests (appended to `tiny_cut_detection_test.cc`)**

```cpp
namespace {
// A "bowtie": two triangles {0,1,2} and {2,3,4} sharing vertex 2. This graph
// has exactly two 2-edge-cut classes: {(1,2),(0,2)} isolating {0,1} from
// {2,3,4}, and {(2,3),(2,4)} isolating {3,4} from {0,1,2}. Edge (0,1) and
// (3,4) are each 3-edge-connected within their own triangle and belong to no
// cut class.
FilterGraph make_bowtie(NodeWeight w0, NodeWeight w1, NodeWeight w2, NodeWeight w3, NodeWeight w4) {
  std::vector<NodeWeight> weights = {w0, w1, w2, w3, w4};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {0, 2, 1},
    {2, 3, 1}, {2, 4, 1}, {3, 4, 1}
  };
  return build_csr_from_edge_list(edges, weights);
}
}  // namespace

TEST(ContractTwoEdgeCutsTest, ContractsBothSmallSidesOfABowtie) {
  FilterGraph graph = make_bowtie(1, 1, 100, 1, 1);
  ContractionResult result = contract_two_edge_cuts(graph, 5);

  EXPECT_EQ(result.graph.numNodes(), 3u);
  EXPECT_EQ(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[3], result.mapping[4]);
  EXPECT_NE(result.mapping[0], result.mapping[2]);
  EXPECT_NE(result.mapping[2], result.mapping[3]);
  EXPECT_NE(result.mapping[0], result.mapping[3]);
}

TEST(ContractTwoEdgeCutsTest, LeavesHeavySideUncontracted) {
  // {0,1} together weigh 100, too heavy to contract at U = 5, but {3,4}
  // (weight 2) still qualifies.
  FilterGraph graph = make_bowtie(50, 50, 1, 1, 1);
  ContractionResult result = contract_two_edge_cuts(graph, 5);

  EXPECT_NE(result.mapping[0], result.mapping[1]);
  EXPECT_EQ(result.mapping[3], result.mapping[4]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `contract_two_edge_cuts` is not declared

- [ ] **Step 3: Add the declaration**

Append to `tiny_cut_detection.h`:

```cpp
// Part 1, pass 3 (design spec section 4.4): finds 2-edge-cut equivalence
// classes and, for each class S, contracts every connected component of
// (V, E \ S) whose total weight is <= U.
ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U);
```

- [ ] **Step 4: Write the implementation**

Append to `tiny_cut_detection.cpp`:

```cpp
ContractionResult contract_two_edge_cuts(const FilterGraph& graph, NodeWeight U) {
  const size_t n = graph.numNodes();
  std::vector<NodeID> component = parallel_connected_components(graph);
  SpanningForest forest = build_spanning_forest(graph, component, {});
  EdgeSignatures sigs = compute_edge_signatures(graph, forest);
  std::vector<std::vector<EdgeID>> classes = find_two_edge_cut_classes(graph, sigs);

  AtomicUnionFind uf(n);
  for (const std::vector<EdgeID>& cls : classes) {
    std::vector<char> excluded(graph.numEdges(), 0);
    for (EdgeID e : cls) excluded[e] = 1;
    std::vector<NodeID> comp = parallel_connected_components_excluding(graph, excluded);

    std::unordered_map<NodeID, NodeWeight> comp_weight;
    for (size_t v = 0; v < n; ++v) comp_weight[comp[v]] += graph.node_weight[v];

    for (size_t v = 0; v < n; ++v) {
      if (comp_weight[comp[v]] <= U) uf.unite(static_cast<NodeID>(v), comp[v]);
    }
  }
  return contract_graph(graph, uf);
}
```

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=ContractTwoEdgeCutsTest.*`
Expected: PASS (2 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp \
        tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc
git commit -m "feat: contract 2-edge-cut equivalence classes (Part 1 pass 3)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 11: Part 1 orchestrator

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp`
- Test: `tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc` (append)

**Interfaces:**
- Produces: `run_tiny_cut_detection(graph, params) -> ContractionResult`, chaining passes 1-3 (Tasks 8-10) and composing their node mappings into one mapping from the *original* input graph's vertices to the final output graph's vertices.

- [ ] **Step 1: Write the failing test (appended to `tiny_cut_detection_test.cc`)**

```cpp
TEST(RunTinyCutDetectionTest, ComposesMappingAcrossAllThreePasses) {
  // A path of 6 light vertices (0..5) hanging off a heavy triangle (6,7,8)
  // via a bridge (5,6). Pass 1 should contract the light path's bridge-
  // isolated subtree if small enough; whatever survives should still respect
  // the invariant that every original vertex maps to some final vertex, and
  // the heavy triangle should never merge with anything (it's the largest
  // component and stays intact since U is too small to touch it).
  std::vector<NodeWeight> weights = {1, 1, 1, 1, 1, 1, 100, 100, 100};
  std::vector<EdgeListEntry> edges = {
    {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 5, 1},
    {5, 6, 1},
    {6, 7, 1}, {7, 8, 1}, {6, 8, 1}
  };
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  ContractionResult result = run_tiny_cut_detection(graph, TinyCutParams{6, 5});

  ASSERT_EQ(result.mapping.size(), 9u);
  for (NodeID v = 0; v < 9; ++v) EXPECT_LT(result.mapping[v], result.graph.numNodes());

  NodeWeight total_original = 0, total_output = 0;
  for (NodeWeight w : weights) total_original += w;
  for (size_t v = 0; v < result.graph.numNodes(); ++v) total_output += result.graph.node_weight[v];
  EXPECT_EQ(total_original, total_output);

  // The heavy triangle vertices must all still be distinguishable from the
  // light path (they're far too heavy, individually, to be absorbed).
  EXPECT_NE(result.mapping[6], result.mapping[0]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `run_tiny_cut_detection` is not declared

- [ ] **Step 3: Add the declaration**

Append to `tiny_cut_detection.h`:

```cpp
// Runs all three tiny-cut passes in sequence (design spec section 4),
// composing their node mappings so the result maps the ORIGINAL input
// graph's vertices directly to the final, smallest output graph's vertices.
ContractionResult run_tiny_cut_detection(const FilterGraph& graph, const TinyCutParams& params);
```

- [ ] **Step 4: Write the implementation**

Append to `tiny_cut_detection.cpp`:

```cpp
ContractionResult run_tiny_cut_detection(const FilterGraph& graph, const TinyCutParams& params) {
  ContractionResult r1 = contract_component_tree(graph, params);
  ContractionResult r2 = contract_degree2_chains(r1.graph, params.U);
  ContractionResult r3 = contract_two_edge_cuts(r2.graph, params.U);

  std::vector<NodeID> composed(graph.numNodes());
  for (size_t v = 0; v < graph.numNodes(); ++v) {
    composed[v] = r3.mapping[r2.mapping[r1.mapping[v]]];
  }

  ContractionResult result;
  result.graph = std::move(r3.graph);
  result.mapping = std::move(composed);
  return result;
}
```

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=RunTinyCutDetectionTest.*`
Expected: PASS (1 test)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.cpp \
        tests/partition/preprocessing/filtering/tiny_cut_detection_test.cc
git commit -m "feat: orchestrate Part 1 tiny-cut detection passes

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

**Checkpoint after Task 11: Part 1 (tiny-cut detection) is complete, tested end-to-end, and CLI-reportable.** This matches the phased-checkpoint pacing: pause here to validate Part 1 in isolation (e.g. via a quick manual CLI smoke test once Task 18 exists, or just trusting the unit test suite) before starting Part 2. Tasks 12-16 build Part 2 (natural-cut detection).

---

### Task 12: Dinic's max-flow / min-cut solver

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/dinic_max_flow_test.cc`

**Interfaces:**
- Produces: `FlowNetwork{adj, source, sink, numNodes(), add_edge(u,v,cap)}`, `dinic_max_flow(network) -> int64_t`, `min_cut_reachable_from_source(network) -> std::vector<char>`.

Design decision (spec section 5): Dinic's (`O(V^2 E)`) is used instead of push-relabel with global relabeling, since this module solves many small, independent local subproblems (bounded by BFS growth to size ~`alpha*U`) rather than one large flow network reused across phases — push-relabel's higher setup cost doesn't pay off at this scale. Swappable later if benchmarking on real DIMACS instances shows otherwise.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/dinic_max_flow_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"

using namespace mt_kahypar::filtering;

TEST(DinicMaxFlowTest, ComputesKnownMaxFlowAndMinCut) {
  // Undirected network: s(0)-a(1) cap 3, s(0)-b(2) cap 2, a(1)-t(3) cap 2,
  // b(2)-t(3) cap 3, a(1)-b(2) cap 1. Hand-computed max flow / min cut = 5
  // (e.g. the cut {s} alone already has capacity 3+2=5, and it's achievable:
  // push s-a=3, s-b=2, then a sends 2 to t and 1 to b, b sends 3 to t).
  struct Edge { uint32_t u, v; int64_t cap; };
  std::vector<Edge> edges = {
    {0, 1, 3}, {0, 2, 2}, {1, 3, 2}, {2, 3, 3}, {1, 2, 1}
  };
  FlowNetwork network(4);
  network.source = 0;
  network.sink = 3;
  for (const Edge& e : edges) network.add_edge(e.u, e.v, e.cap);

  const int64_t max_flow = dinic_max_flow(network);
  EXPECT_EQ(max_flow, 5);

  std::vector<char> reachable = min_cut_reachable_from_source(network);
  EXPECT_TRUE(reachable[0]);
  EXPECT_FALSE(reachable[3]);

  int64_t crossing_capacity = 0;
  for (const Edge& e : edges) {
    if (reachable[e.u] != reachable[e.v]) crossing_capacity += e.cap;
  }
  EXPECT_EQ(crossing_capacity, max_flow);
}

TEST(DinicMaxFlowTest, DisconnectedSourceAndSinkGiveZeroFlow) {
  FlowNetwork network(2);
  network.source = 0;
  network.sink = 1;
  EXPECT_EQ(dinic_max_flow(network), 0);
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/dinic_max_flow_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `dinic_max_flow.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h
#pragma once

#include <cstdint>
#include <vector>

namespace mt_kahypar {
namespace filtering {

// A small local flow network used for one natural-cut subproblem. Nodes are
// local ids [0, numNodes()); `source` and `sink` designate s and t.
struct FlowNetwork {
  struct Arc {
    uint32_t to;
    int64_t capacity;
    uint32_t reverse_index;  // index of the reverse arc in adj[to]
  };

  std::vector<std::vector<Arc>> adj;
  uint32_t source = 0;
  uint32_t sink = 0;

  explicit FlowNetwork(uint32_t num_nodes) : adj(num_nodes) {}

  uint32_t numNodes() const { return static_cast<uint32_t>(adj.size()); }

  // Adds a pair of forward/backward arcs with equal capacity, correctly
  // modeling one undirected edge of capacity `cap` in a directed max-flow
  // formulation (standard technique).
  void add_edge(uint32_t u, uint32_t v, int64_t cap) {
    const uint32_t fwd_index = static_cast<uint32_t>(adj[u].size());
    const uint32_t bwd_index = static_cast<uint32_t>(adj[v].size());
    adj[u].push_back(Arc{v, cap, bwd_index});
    adj[v].push_back(Arc{u, cap, fwd_index});
  }
};

// Runs Dinic's max-flow algorithm and returns the flow value. Mutates
// `network`'s arc capacities in place (to residual capacities).
int64_t dinic_max_flow(FlowNetwork& network);

// Must be called after dinic_max_flow(network) on the same network. Returns,
// per node, whether it is reachable from `source` in the final residual
// graph -- the min-cut separates reachable from unreachable nodes.
std::vector<char> min_cut_reachable_from_source(const FlowNetwork& network);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.cpp
#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"

#include <algorithm>
#include <deque>
#include <limits>

namespace mt_kahypar {
namespace filtering {

namespace {
bool bfs_level_graph(FlowNetwork& network, std::vector<int>& level) {
  level.assign(network.numNodes(), -1);
  level[network.source] = 0;
  std::deque<uint32_t> queue{network.source};
  while (!queue.empty()) {
    const uint32_t u = queue.front();
    queue.pop_front();
    for (const auto& arc : network.adj[u]) {
      if (arc.capacity > 0 && level[arc.to] == -1) {
        level[arc.to] = level[u] + 1;
        queue.push_back(arc.to);
      }
    }
  }
  return level[network.sink] != -1;
}

int64_t send_flow(FlowNetwork& network, const std::vector<int>& level, std::vector<size_t>& it,
                   uint32_t u, int64_t pushed) {
  if (u == network.sink) return pushed;
  for (; it[u] < network.adj[u].size(); ++it[u]) {
    auto& arc = network.adj[u][it[u]];
    if (arc.capacity > 0 && level[arc.to] == level[u] + 1) {
      const int64_t bottleneck = send_flow(network, level, it, arc.to, std::min(pushed, arc.capacity));
      if (bottleneck > 0) {
        arc.capacity -= bottleneck;
        network.adj[arc.to][arc.reverse_index].capacity += bottleneck;
        return bottleneck;
      }
    }
  }
  return 0;
}
}  // namespace

int64_t dinic_max_flow(FlowNetwork& network) {
  int64_t total_flow = 0;
  std::vector<int> level;
  while (bfs_level_graph(network, level)) {
    std::vector<size_t> it(network.numNodes(), 0);
    int64_t pushed;
    while ((pushed = send_flow(network, level, it, network.source,
                                std::numeric_limits<int64_t>::max())) > 0) {
      total_flow += pushed;
    }
  }
  return total_flow;
}

std::vector<char> min_cut_reachable_from_source(const FlowNetwork& network) {
  std::vector<char> reachable(network.numNodes(), 0);
  reachable[network.source] = 1;
  std::deque<uint32_t> queue{network.source};
  while (!queue.empty()) {
    const uint32_t u = queue.front();
    queue.pop_front();
    for (const auto& arc : network.adj[u]) {
      if (arc.capacity > 0 && !reachable[arc.to]) {
        reachable[arc.to] = 1;
        queue.push_back(arc.to);
      }
    }
  }
  return reachable;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/dinic_max_flow.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=DinicMaxFlowTest.*`
Expected: PASS (2 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h \
        mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/dinic_max_flow_test.cc
git commit -m "feat: add Dinic's max-flow solver for local natural-cut subproblems

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 13: BFS core/ring construction + single-seed natural cut

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/natural_cut_detection_test.cc`

**Interfaces:**
- Consumes: `FilterGraph` (Task 1), `FlowNetwork`, `dinic_max_flow`, `min_cut_reachable_from_source` (Task 12).
- Produces: `NaturalCutParams{U, alpha=1.0, f=10.0, coverage=2}`, `NaturalCutScratch{bfs_queue, in_tree, in_core, tree_order}`, `compute_natural_cut(graph, seed, params, scratch, covered) -> std::vector<EdgeID>`.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/natural_cut_detection_test.cc
#include <gtest/gtest.h>

#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

using namespace mt_kahypar::filtering;

TEST(ComputeNaturalCutTest, FindsTheOnlyPossibleCutOnAUnitPath) {
  // Path 0-1-2-3-4, unit weights and unit edge weights. With U = 3,
  // alpha = 1.0: BFS from seed 0 grows the tree until size >= 3, i.e.
  // tree = {0,1,2}. With f = 10: core-size threshold = 3/10 = 0.3, so core =
  // {0} only (0's own weight of 1 already exceeds 0.3). ring = neighbors of
  // the tree outside it = {3}. The local network is a unit-capacity path
  // s(0) - 1 - 2 - t(3), whose unique min cut (value 1) is the edge nearest
  // the source once Dinic's saturates it -- edge (0,1).
  std::vector<NodeWeight> weights(5, 1);
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  NaturalCutParams params{3, 1.0, 10.0, 2};
  NaturalCutScratch scratch;
  std::vector<char> covered(5, 0);

  std::vector<EdgeID> cut = compute_natural_cut(graph, 0, params, scratch, covered);

  ASSERT_EQ(cut.size(), 1u);
  EXPECT_EQ(cut[0], 0u);  // canonical edge id of (0,1), the first edge added

  EXPECT_TRUE(covered[0]);
  EXPECT_TRUE(covered[1]);
  EXPECT_TRUE(covered[2]);
  EXPECT_FALSE(covered[3]);  // ring, not tree -- not marked covered
  EXPECT_FALSE(covered[4]);  // never visited
}
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/natural_cut_detection_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `natural_cut_detection.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h
#pragma once

#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/dinic_max_flow.h"
#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct NaturalCutParams {
  NodeWeight U;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
};

// Scratch buffers reused across many BFS-growth + local-min-cut calls, to
// avoid per-call heap allocation (design spec section 5). One instance per
// thread when parallelized (Task 15).
struct NaturalCutScratch {
  std::vector<NodeID> bfs_queue;
  std::vector<char> in_tree;
  std::vector<char> in_core;
  std::vector<NodeID> tree_order;
};

// Grows a BFS tree from `seed` until its total vertex weight reaches
// params.U * params.alpha, computes the local min s-t cut between the
// resulting core (contracted to s) and ring (contracted to t), and returns
// the canonical edge ids of `graph` forming that min cut. Every vertex
// visited by the BFS growth (the tree) is marked covered[v] = true; ring
// vertices are not marked (design spec section 5). `covered` must have size
// graph.numNodes().
std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                         const NaturalCutParams& params,
                                         NaturalCutScratch& scratch,
                                         std::vector<char>& covered);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

std::vector<EdgeID> compute_natural_cut(const FilterGraph& graph, NodeID seed,
                                         const NaturalCutParams& params,
                                         NaturalCutScratch& scratch,
                                         std::vector<char>& covered) {
  const size_t n = graph.numNodes();
  scratch.in_tree.assign(n, 0);
  scratch.in_core.assign(n, 0);
  scratch.tree_order.clear();
  scratch.bfs_queue.clear();

  const auto target_tree_size = static_cast<NodeWeight>(params.alpha * params.U);
  const auto target_core_size = static_cast<NodeWeight>(params.alpha * params.U / params.f);

  NodeWeight tree_size = 0;
  size_t head = 0;
  scratch.bfs_queue.push_back(seed);
  scratch.in_tree[seed] = 1;
  scratch.tree_order.push_back(seed);
  tree_size += graph.node_weight[seed];
  covered[seed] = 1;

  while (head < scratch.bfs_queue.size() && tree_size < target_tree_size) {
    const NodeID u = scratch.bfs_queue[head++];
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      if (tree_size >= target_tree_size) break;
      const NodeID v = graph.adj[pos];
      if (scratch.in_tree[v]) continue;
      scratch.in_tree[v] = 1;
      scratch.tree_order.push_back(v);
      scratch.bfs_queue.push_back(v);
      tree_size += graph.node_weight[v];
      covered[v] = 1;
    }
  }

  NodeWeight running = 0;
  for (NodeID v : scratch.tree_order) {
    if (running >= target_core_size) break;
    scratch.in_core[v] = 1;
    running += graph.node_weight[v];
  }

  std::vector<NodeID> ring;
  std::vector<char> in_ring(n, 0);
  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (!scratch.in_tree[v] && !in_ring[v]) {
        in_ring[v] = 1;
        ring.push_back(v);
      }
    }
  }

  // local id 0 = s (core), local id 1 = t (ring), 2.. = tree-interior.
  std::vector<NodeID> local_id(n, kInvalidNode);
  uint32_t next_local_id = 2;
  for (NodeID v : scratch.tree_order) {
    local_id[v] = scratch.in_core[v] ? 0 : next_local_id++;
  }
  for (NodeID v : ring) local_id[v] = 1;

  FlowNetwork network(next_local_id);
  network.source = 0;
  network.sink = 1;
  std::vector<std::vector<EdgeID>> local_arc_to_original(next_local_id);

  auto add_local_edge = [&](NodeID lu, NodeID lv, EdgeID original_edge, EdgeWeight weight) {
    if (lu == lv) return;  // both endpoints collapsed to the same local node
    network.add_edge(lu, lv, weight);
    local_arc_to_original[lu].push_back(original_edge);
    local_arc_to_original[lv].push_back(original_edge);
  };

  for (NodeID u : scratch.tree_order) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (local_id[v] == kInvalidNode) continue;
      if (scratch.in_tree[v] && v <= u) continue;  // dedup tree-tree edges
      add_local_edge(local_id[u], local_id[v], graph.adj_edge[pos], graph.edge_weight[graph.adj_edge[pos]]);
    }
  }
  for (NodeID u : ring) {
    for (EdgeID pos = graph.node_begin[u]; pos < graph.node_begin[u + 1]; ++pos) {
      const NodeID v = graph.adj[pos];
      if (local_id[v] == kInvalidNode || local_id[v] == 1) continue;  // outside net, or ring-ring
      if (scratch.in_tree[v]) continue;  // handled from the tree side above
      add_local_edge(local_id[u], local_id[v], graph.adj_edge[pos], graph.edge_weight[graph.adj_edge[pos]]);
    }
  }

  dinic_max_flow(network);
  std::vector<char> reachable = min_cut_reachable_from_source(network);

  std::vector<EdgeID> cut_edges;
  for (uint32_t lu = 0; lu < network.numNodes(); ++lu) {
    if (!reachable[lu]) continue;
    for (size_t i = 0; i < network.adj[lu].size(); ++i) {
      if (!reachable[network.adj[lu][i].to]) cut_edges.push_back(local_arc_to_original[lu][i]);
    }
  }
  return cut_edges;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

Note the `v <= u` dedup check in the tree-tree loop: it relies on tree vertices being compared by their original `NodeID` value, which is a valid (if slightly arbitrary) way to visit each tree-internal edge exactly once from a consistent side; it does not depend on `tree_order` position.

- [ ] **Step 5: Register the source, build, run**

Add `filtering/natural_cut_detection.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=ComputeNaturalCutTest.*`
Expected: PASS (1 test)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/natural_cut_detection_test.cc
git commit -m "feat: compute a single natural cut from BFS core/ring construction

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 14: Sequential seed loop (coverage sweeps)

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp`
- Test: `tests/partition/preprocessing/filtering/natural_cut_detection_test.cc` (append)

A sequential baseline matching PUNCH's own description ("we first pick all centers sequentially, for simplicity") — establishes correct sweep/coverage semantics before parallelizing in Task 15.

**Interfaces:**
- Produces: `run_natural_cut_detection_sequential(graph, params, rng) -> std::vector<char>` (size `numEdges()`, `keep[e] == 1` iff `e` is a kept natural-cut edge). `keep` accumulates across all `params.coverage` sweeps; only the per-sweep `covered` set resets between sweeps.

- [ ] **Step 1: Write the failing tests (appended to `natural_cut_detection_test.cc`)**

```cpp
#include <random>

TEST(RunNaturalCutDetectionSequentialTest, WholeSmallGraphAbsorbedYieldsNoCuts) {
  // A single triangle, U so large the whole component fits inside one BFS
  // tree: there is no "outside" to cut against, so no edges get kept.
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::mt19937_64 rng(42);
  std::vector<char> keep = run_natural_cut_detection_sequential(graph, NaturalCutParams{1000, 1.0, 10.0, 2}, rng);
  ASSERT_EQ(keep.size(), 3u);
  for (char k : keep) EXPECT_EQ(k, 0);
}

TEST(RunNaturalCutDetectionSequentialTest, LongPathKeepsSomeEdges) {
  // A path much longer than U forces multiple natural cuts along its length.
  const int len = 40;
  std::vector<NodeWeight> weights(len, 1);
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i + 1 < len; ++i) edges.push_back({static_cast<NodeID>(i), static_cast<NodeID>(i + 1), 1});
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::mt19937_64 rng(7);
  std::vector<char> keep = run_natural_cut_detection_sequential(graph, NaturalCutParams{5, 1.0, 10.0, 2}, rng);
  ASSERT_EQ(keep.size(), static_cast<size_t>(len - 1));
  bool any_kept = false;
  for (char k : keep) any_kept = any_kept || k;
  EXPECT_TRUE(any_kept);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `run_natural_cut_detection_sequential` is not declared

- [ ] **Step 3: Add the declaration**

Append to `natural_cut_detection.h` (add `#include <random>` to its includes):

```cpp
// Runs the full sequential natural-cut detection procedure (design spec
// section 5): for each of params.coverage sweeps, resets per-sweep coverage,
// visits vertices in a freshly shuffled order, and for every not-yet-covered
// vertex runs compute_natural_cut as a new seed, accumulating its cut edges
// into the returned keep set (which persists across all sweeps).
std::vector<char> run_natural_cut_detection_sequential(const FilterGraph& graph,
                                                        const NaturalCutParams& params,
                                                        std::mt19937_64& rng);
```

- [ ] **Step 4: Write the implementation**

Append to `natural_cut_detection.cpp` (add `#include <algorithm>` to its includes):

```cpp
std::vector<char> run_natural_cut_detection_sequential(const FilterGraph& graph,
                                                        const NaturalCutParams& params,
                                                        std::mt19937_64& rng) {
  const size_t n = graph.numNodes();
  std::vector<char> keep(graph.numEdges(), 0);
  NaturalCutScratch scratch;

  std::vector<NodeID> order(n);
  for (size_t v = 0; v < n; ++v) order[v] = static_cast<NodeID>(v);

  for (int sweep = 0; sweep < params.coverage; ++sweep) {
    std::vector<char> covered(n, 0);
    std::shuffle(order.begin(), order.end(), rng);
    for (NodeID v : order) {
      if (covered[v]) continue;
      std::vector<EdgeID> cut = compute_natural_cut(graph, v, params, scratch, covered);
      for (EdgeID e : cut) keep[e] = 1;
    }
  }
  return keep;
}
```

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=RunNaturalCutDetectionSequentialTest.*`
Expected: PASS (2 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp \
        tests/partition/preprocessing/filtering/natural_cut_detection_test.cc
git commit -m "feat: add sequential coverage-sweep loop for natural-cut detection

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 15: Parallel seed loop

**Files:**
- Modify: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h`
- Modify: `mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp`
- Test: `tests/partition/preprocessing/filtering/natural_cut_detection_test.cc` (append)

This is the new parallel-seed-scheduling design called for by the spec (PUNCH's own scheduling is sequential "for simplicity" — see design spec section 5): atomic CAS-based seed claiming, atomic covered marking, thread-local scratch, relaxed atomic `keep` flags.

**Interfaces:**
- Produces: `run_natural_cut_detection(graph, params) -> std::vector<char>` (size `numEdges()`), the parallel counterpart to Task 14's sequential version, same semantics.

- [ ] **Step 1: Write the failing tests (appended to `natural_cut_detection_test.cc`)**

```cpp
TEST(RunNaturalCutDetectionParallelTest, WholeSmallGraphAbsorbedYieldsNoCuts) {
  std::vector<NodeWeight> weights = {1, 1, 1};
  std::vector<EdgeListEntry> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 1}};
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  std::vector<char> keep = run_natural_cut_detection(graph, NaturalCutParams{1000, 1.0, 10.0, 2});
  ASSERT_EQ(keep.size(), 3u);
  for (char k : keep) EXPECT_EQ(k, 0);
}

TEST(RunNaturalCutDetectionParallelTest, LongPathKeepsSomeEdgesAndNeverCrashes) {
  const int len = 200;
  std::vector<NodeWeight> weights(len, 1);
  std::vector<EdgeListEntry> edges;
  for (int i = 0; i + 1 < len; ++i) edges.push_back({static_cast<NodeID>(i), static_cast<NodeID>(i + 1), 1});
  FilterGraph graph = build_csr_from_edge_list(edges, weights);

  for (int trial = 0; trial < 5; ++trial) {
    std::vector<char> keep = run_natural_cut_detection(graph, NaturalCutParams{5, 1.0, 10.0, 2});
    ASSERT_EQ(keep.size(), static_cast<size_t>(len - 1));
    bool any_kept = false;
    for (char k : keep) any_kept = any_kept || k;
    EXPECT_TRUE(any_kept);
  }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `run_natural_cut_detection` is not declared

- [ ] **Step 3: Add the declaration**

Append to `natural_cut_detection.h`:

```cpp
// Parallel counterpart to run_natural_cut_detection_sequential, with the
// same semantics (design spec section 5): per sweep, a pre-shuffled vertex
// order is scanned by a TBB parallel_for; each task atomically claims its
// vertex as a new seed via compare-exchange on a per-vertex `covered` flag
// (skipping on failure), then runs compute_natural_cut using thread-local
// scratch buffers. `keep` flags are set with plain relaxed atomic stores
// (a monotonic boolean needs no compare-exchange). A short serial mop-up
// pass handles any vertices left uncovered by races near the scan's end.
std::vector<char> run_natural_cut_detection(const FilterGraph& graph, const NaturalCutParams& params);
```

- [ ] **Step 4: Write the implementation**

Append to `natural_cut_detection.cpp` (add includes `<tbb/parallel_for.h>`, `"mt-kahypar/parallel/atomic_wrapper.h"`, `"mt-kahypar/parallel/stl/thread_locals.h"`):

```cpp
std::vector<char> run_natural_cut_detection(const FilterGraph& graph, const NaturalCutParams& params) {
  const size_t n = graph.numNodes();
  const size_t m = graph.numEdges();
  std::vector<parallel::IntegralAtomicWrapper<uint8_t>> keep(m);
  for (size_t e = 0; e < m; ++e) keep[e] = 0;

  tls_enumerable_thread_specific<NaturalCutScratch> scratch;

  std::vector<NodeID> order(n);
  for (size_t v = 0; v < n; ++v) order[v] = static_cast<NodeID>(v);

  for (int sweep = 0; sweep < params.coverage; ++sweep) {
    std::vector<parallel::IntegralAtomicWrapper<uint8_t>> covered(n);
    for (size_t v = 0; v < n; ++v) covered[v] = 0;

    std::mt19937_64 shuffle_rng(0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(sweep));
    std::shuffle(order.begin(), order.end(), shuffle_rng);

    // A "wide" covered marker vector usable from compute_natural_cut, which
    // expects std::vector<char>&. We bridge via a thread-local plain vector
    // that mirrors the atomic one only for the duration of one seed's BFS,
    // then flushes newly-visited vertices back into the atomic array.
    tbb::parallel_for(size_t(0), order.size(), [&](size_t i) {
      const NodeID v = order[i];
      uint8_t expected = 0;
      if (!covered[v].compare_exchange_strong(expected, 1)) return;

      NaturalCutScratch& local_scratch = scratch.local();
      std::vector<char> covered_snapshot(n, 0);
      for (size_t u = 0; u < n; ++u) covered_snapshot[u] = covered[u].load(std::memory_order_relaxed);

      std::vector<EdgeID> cut = compute_natural_cut(graph, v, params, local_scratch, covered_snapshot);

      for (size_t u = 0; u < n; ++u) {
        if (covered_snapshot[u]) covered[u].store(1, std::memory_order_relaxed);
      }
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    });

    // Mop-up pass: handles any vertex left uncovered by races near the end
    // of the parallel scan (covered is monotonic, so this terminates).
    for (NodeID v : order) {
      if (covered[v].load(std::memory_order_relaxed)) continue;
      covered[v].store(1, std::memory_order_relaxed);
      NaturalCutScratch& local_scratch = scratch.local();
      std::vector<char> covered_snapshot(n, 0);
      for (size_t u = 0; u < n; ++u) covered_snapshot[u] = covered[u].load(std::memory_order_relaxed);
      std::vector<EdgeID> cut = compute_natural_cut(graph, v, params, local_scratch, covered_snapshot);
      for (size_t u = 0; u < n; ++u) {
        if (covered_snapshot[u]) covered[u].store(1, std::memory_order_relaxed);
      }
      for (EdgeID e : cut) keep[e].store(1, std::memory_order_relaxed);
    }
  }

  std::vector<char> result(m);
  for (size_t e = 0; e < m; ++e) result[e] = keep[e].load(std::memory_order_relaxed);
  return result;
}
```

**Note on the `covered_snapshot` copy:** `compute_natural_cut` takes `std::vector<char>&`, but the parallel loop's shared state is an atomic array (needed for the CAS-based seed claim). Rather than changing `compute_natural_cut`'s signature (which would force every caller, including the simple sequential Task 14 path, to deal with atomics), each task takes an O(n) snapshot of `covered` into a plain vector, lets `compute_natural_cut` mark it freely, then flushes the newly-true entries back atomically. This is O(n) extra work per seed, which is acceptable for a first correctness-focused parallel implementation; if profiling on real DIMACS instances shows this dominates, the follow-up optimization is to change `compute_natural_cut` to accept an abstract "mark covered" callback instead of a concrete vector reference, avoiding the snapshot copy entirely. Document this as a deferred optimization, matching this plan's other explicit, benchmarking-gated design decisions.

- [ ] **Step 5: Build and run**

Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=RunNaturalCutDetectionParallelTest.*`
Expected: PASS (2 tests, including the 5-trial loop)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h \
        mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.cpp \
        tests/partition/preprocessing/filtering/natural_cut_detection_test.cc
git commit -m "feat: parallelize natural-cut seed scheduling with TBB

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

**Checkpoint after Task 15: Part 2 (natural-cut detection) is complete and tested, sequential and parallel.** Tasks 16-19 build Part 3 (assembly), the CLI, and validation.

---

### Task 16: Fragment assembly

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/fragment_assembly.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/fragment_assembly_test.cc`

**Interfaces:**
- Consumes: `FilterGraph`, `AtomicUnionFind` (Task 3).
- Produces: `FilteringResult{fragment_id (indexed by ORIGINAL vertex id), fragment_size (indexed by fragment id), kept_edges (pairs of representative ORIGINAL vertex ids)}`, `assemble_fragments(graph, keep, part1_mapping) -> FilteringResult`.

`graph` here is Part 1's output graph; `keep` is Part 2's per-edge flag over `graph`; `part1_mapping` is Part 1's composed mapping (original vertex id -> `graph` vertex id, from Task 11). A representative original vertex is chosen per `graph` vertex (its lowest original id) so `kept_edges` can be translated back for later plotting against DIMACS coordinates.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/fragment_assembly_test.cc
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
```

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/fragment_assembly_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `fragment_assembly.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h
#pragma once

#include <utility>
#include <vector>

#include "mt-kahypar/partition/preprocessing/filtering/filter_graph.h"

namespace mt_kahypar {
namespace filtering {

struct FilteringResult {
  std::vector<NodeID> fragment_id;                    // indexed by ORIGINAL vertex id
  std::vector<NodeWeight> fragment_size;               // indexed by fragment id
  std::vector<std::pair<NodeID, NodeID>> kept_edges;   // representative ORIGINAL vertex id pairs
};

// Part 3 (design spec section 6): unions the endpoints of every edge of
// `graph` NOT flagged in `keep`; the resulting connected components are the
// fragments. `graph` is Part 1's output graph, `keep` is Part 2's per-edge
// flag over it (size graph.numEdges()), and `part1_mapping` is Part 1's
// composed mapping from ORIGINAL vertex ids to `graph` vertex ids. Kept
// edges are translated back to a representative original vertex per `graph`
// vertex (its lowest original id), for later plotting against DIMACS
// coordinates.
FilteringResult assemble_fragments(const FilterGraph& graph,
                                    const std::vector<char>& keep,
                                    const std::vector<NodeID>& part1_mapping);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/fragment_assembly.cpp
#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"

#include "mt-kahypar/partition/preprocessing/filtering/union_find.h"

namespace mt_kahypar {
namespace filtering {

FilteringResult assemble_fragments(const FilterGraph& graph,
                                    const std::vector<char>& keep,
                                    const std::vector<NodeID>& part1_mapping) {
  const size_t graph_n = graph.numNodes();
  const size_t orig_n = part1_mapping.size();

  AtomicUnionFind uf(graph_n);
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;
      if (!keep[graph.adj_edge[pos]]) uf.unite(v, other);
    }
  }

  std::vector<NodeID> fragment_of_graph_vertex(graph_n, kInvalidNode);
  NodeID next_fragment = 0;
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    const uint32_t rep = uf.find(v);
    if (fragment_of_graph_vertex[rep] == kInvalidNode) fragment_of_graph_vertex[rep] = next_fragment++;
    fragment_of_graph_vertex[v] = fragment_of_graph_vertex[rep];
  }

  FilteringResult result;
  result.fragment_id.resize(orig_n);
  for (size_t orig = 0; orig < orig_n; ++orig) {
    result.fragment_id[orig] = fragment_of_graph_vertex[part1_mapping[orig]];
  }

  // Fragment sizes are the sum of Part-1-output vertex weights per fragment,
  // which already equals the sum of ORIGINAL vertex weights, since Part 1's
  // contraction sums weights transitively (see graph_contraction.cpp) --
  // original weights aren't needed directly here.
  result.fragment_size.assign(next_fragment, 0);
  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    result.fragment_size[fragment_of_graph_vertex[v]] += graph.node_weight[v];
  }

  std::vector<NodeID> representative(graph_n, kInvalidNode);
  for (size_t orig = 0; orig < orig_n; ++orig) {
    NodeID& rep = representative[part1_mapping[orig]];
    if (rep == kInvalidNode) rep = static_cast<NodeID>(orig);
  }

  for (NodeID v = 0; v < static_cast<NodeID>(graph_n); ++v) {
    for (EdgeID pos = graph.node_begin[v]; pos < graph.node_begin[v + 1]; ++pos) {
      const NodeID other = graph.adj[pos];
      if (other <= v) continue;
      if (keep[graph.adj_edge[pos]]) {
        result.kept_edges.emplace_back(representative[v], representative[other]);
      }
    }
  }

  return result;
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/fragment_assembly.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=AssembleFragmentsTest.*`
Expected: PASS (2 tests)

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h \
        mt-kahypar/partition/preprocessing/filtering/fragment_assembly.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/fragment_assembly_test.cc
git commit -m "feat: assemble final fragments from kept natural-cut edges

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 17: Pipeline orchestrator + `U`-invariant

**Files:**
- Create: `mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h`
- Create: `mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.cpp`
- Modify: `mt-kahypar/partition/preprocessing/CMakeLists.txt`
- Modify: `tests/partition/preprocessing/CMakeLists.txt`
- Test: `tests/partition/preprocessing/filtering/filtering_invariant_test.cc`

**Interfaces:**
- Consumes: `run_tiny_cut_detection` (Task 11), `run_natural_cut_detection` (Task 15), `assemble_fragments` (Task 16).
- Produces: `FilteringParams{U, tau=5, alpha=1.0, f=10, coverage=2}`, `run_filtering_pipeline(graph, params) -> FilteringResult`.

This header is the documented extension point (design spec section 2/10) for a future Mt-KaHyPar coarsening integration: a caller only needs to construct a `FilterGraph`, call `run_filtering_pipeline`, and consume `FilteringResult.fragment_id` — nothing else in this module needs to be understood.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/partition/preprocessing/filtering/filtering_invariant_test.cc
#include <gtest/gtest.h>

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
```

Add `#include <algorithm>` to the test file for `std::shuffle`.

- [ ] **Step 2: Register the test file and run to verify failure**

Add `filtering/filtering_invariant_test.cc` to `tests/partition/preprocessing/CMakeLists.txt`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc)`
Expected: FAIL — `filtering_pipeline.h: No such file or directory`

- [ ] **Step 3: Write the header**

```cpp
// mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h
#pragma once

#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"
#include "mt-kahypar/partition/preprocessing/filtering/natural_cut_detection.h"
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

namespace mt_kahypar {
namespace filtering {

// -----------------------------------------------------------------------
// Extension point for a future Mt-KaHyPar coarsening integration: construct
// a FilterGraph (see filter_graph.h) from whatever hypergraph/graph
// representation the caller has, call run_filtering_pipeline, and consume
// FilteringResult.fragment_id as an alternative clustering to feed into
// coarsening -- nothing else in this module needs to be understood. No such
// integration is wired up yet (design spec section 10).
// -----------------------------------------------------------------------
struct FilteringParams {
  NodeWeight U;
  NodeWeight tau = 5;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
};

// Runs the full PUNCH filtering pipeline (Parts 1-3) on `graph`. No fragment
// in the result exceeds params.U, guaranteed by alpha <= 1 (design spec
// section 6) -- enforced as a test invariant, not merely documented.
FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params);

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 4: Write the implementation**

```cpp
// mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.cpp
#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

namespace mt_kahypar {
namespace filtering {

FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params) {
  ContractionResult tiny_cut_result =
      run_tiny_cut_detection(graph, TinyCutParams{params.U, params.tau});

  std::vector<char> keep = run_natural_cut_detection(
      tiny_cut_result.graph,
      NaturalCutParams{params.U, params.alpha, params.f, params.coverage});

  return assemble_fragments(tiny_cut_result.graph, keep, tiny_cut_result.mapping);
}

}  // namespace filtering
}  // namespace mt_kahypar
```

- [ ] **Step 5: Register the source, build, run**

Add `filtering/filtering_pipeline.cpp` to `PreprocessingSources`.
Run: `cmake --build build --target mtkahypar_tests -j$(nproc) && ./build/tests/mtkahypar_tests --gtest_filter=FilteringPipelineTest.*`
Expected: PASS (2 tests)

If `NoFragmentEverExceedsU` fails for some seed, that is a real bug to fix before moving on (this is the plan's hard correctness invariant from the spec) — do not weaken the test.

- [ ] **Step 6: Commit**

```bash
git add mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h \
        mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.cpp \
        mt-kahypar/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/CMakeLists.txt \
        tests/partition/preprocessing/filtering/filtering_invariant_test.cc
git commit -m "feat: orchestrate the full PUNCH filtering pipeline

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 18: CLI tool

**Files:**
- Create: `tools/punch_filter.cc`
- Modify: `tools/CMakeLists.txt`

**Interfaces:**
- Consumes: `read_dimacs_graph` (Task 2), `run_filtering_pipeline`, `FilteringParams`, `FilteringResult` (Task 17).

- [ ] **Step 1: Write the CLI**

```cpp
// tools/punch_filter.cc
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/dimacs_io.h"
#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"

using namespace mt_kahypar::filtering;

namespace {
void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " --graph <dimacs.gr> --U <size> "
            << "[--alpha 1.0] [--f 10] [--coverage 2] [--tau 5] "
            << "[--dump-kept-edges <path>]\n";
}
}  // namespace

int main(int argc, char** argv) {
  std::string graph_path;
  std::string dump_path;
  NodeWeight U = 0;
  NodeWeight tau = 5;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  bool has_U = false;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) { print_usage(argv[0]); std::exit(1); }
      return argv[++i];
    };
    if (arg == "--graph") graph_path = next();
    else if (arg == "--U") { U = std::stoull(next()); has_U = true; }
    else if (arg == "--alpha") alpha = std::stod(next());
    else if (arg == "--f") f = std::stod(next());
    else if (arg == "--coverage") coverage = std::stoi(next());
    else if (arg == "--tau") tau = std::stoull(next());
    else if (arg == "--dump-kept-edges") dump_path = next();
    else { print_usage(argv[0]); return 1; }
  }

  if (graph_path.empty() || !has_U) {
    print_usage(argv[0]);
    return 1;
  }

  FilterGraph graph = read_dimacs_graph(graph_path);
  std::cout << "Loaded graph: |V| = " << graph.numNodes() << ", |E| = " << graph.numEdges() << "\n";

  const auto start = std::chrono::steady_clock::now();
  FilteringResult result = run_filtering_pipeline(graph, FilteringParams{U, tau, alpha, f, coverage});
  const auto end = std::chrono::steady_clock::now();
  const double seconds = std::chrono::duration<double>(end - start).count();

  std::cout << "Filtering complete in " << seconds << "s\n";
  std::cout << "|V'| (fragments) = " << result.fragment_size.size() << "\n";

  NodeWeight min_size = result.fragment_size.empty() ? 0 : result.fragment_size[0];
  NodeWeight max_size = 0;
  size_t near_cap = 0;
  for (NodeWeight size : result.fragment_size) {
    min_size = std::min(min_size, size);
    max_size = std::max(max_size, size);
    if (size > static_cast<NodeWeight>(0.9 * static_cast<double>(U))) ++near_cap;
  }
  std::cout << "fragment size: min=" << min_size << " max=" << max_size
            << " count>0.9*U=" << near_cap << "\n";

  if (!dump_path.empty()) {
    std::ofstream out(dump_path);
    for (const auto& [u, v] : result.kept_edges) out << u << " " << v << "\n";
    std::cout << "Wrote " << result.kept_edges.size() << " kept edges to " << dump_path << "\n";
  }

  return 0;
}
```

- [ ] **Step 2: Register the executable**

Append to `tools/CMakeLists.txt` (in the "Tools with more dependencies" section):

```cmake
add_executable(PunchFilter punch_filter.cc)
target_link_libraries(PunchFilter MtKaHyPar-BuildTools)
```

- [ ] **Step 3: Build**

Run: `cmake --build build --target PunchFilter -j$(nproc)`
Expected: builds successfully

- [ ] **Step 4: Manual smoke test**

```bash
cat > /tmp/punch_smoke.gr << 'EOF'
p sp 6 12
a 1 2 1
a 2 1 1
a 2 3 1
a 3 2 1
a 3 1 1
a 1 3 1
a 3 4 1
a 4 3 1
a 4 5 1
a 5 4 1
a 5 6 1
a 6 5 1
EOF
./build/tools/PunchFilter --graph /tmp/punch_smoke.gr --U 3 --dump-kept-edges /tmp/punch_smoke_kept.txt
```

Expected output: a line reporting `|V| = 6, |E| = 6`, a completion time, a `|V'|` fragment count, a fragment size summary line, and a "Wrote N kept edges" line. Confirm the process exits with code 0 (`echo $?`) and `/tmp/punch_smoke_kept.txt` was created (may be empty if no edges were kept, which is fine for this small a graph).

- [ ] **Step 5: Commit**

```bash
git add tools/punch_filter.cc tools/CMakeLists.txt
git commit -m "feat: add punch_filter CLI for the filtering module

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 19: Validation against small DIMACS road instances

Not a unit test (network-dependent) — a manual validation step recorded here per the design spec's validation plan (section 9), to run once Tasks 1-18 are complete and committed.

**Files:** none created; this task only runs the CLI built in Task 18 against downloaded data and records results in a follow-up note (e.g. a comment on the relevant tracking issue, or directly reported back in conversation — this plan does not mandate a specific place, since it's a one-time validation run, not a repo artifact).

- [ ] **Step 1: Download small DIMACS road instances**

```bash
mkdir -p /tmp/dimacs
for region in NY BAY COL; do
  curl -fSL "http://www.diag.uniroma1.it/challenge9/data/USA-road-d/USA-road-d.${region}.gr.gz" \
    -o "/tmp/dimacs/USA-road-d.${region}.gr.gz"
  gunzip -k "/tmp/dimacs/USA-road-d.${region}.gr.gz"
done
```

If this URL has moved, search for "9th DIMACS Implementation Challenge shortest paths download" to find the current location (confirmed live during this plan's design phase at `http://www.diag.uniroma1.it/challenge9/data/USA-road-d/`).

- [ ] **Step 2: Run the CLI at a few `U` values per instance**

```bash
for region in NY BAY COL; do
  for U in 1000 10000 100000; do
    echo "=== ${region}, U=${U} ==="
    ./build/tools/PunchFilter --graph "/tmp/dimacs/USA-road-d.${region}.gr" --U "$U"
  done
done
```

- [ ] **Step 3: Sanity-check the results**

For each run, confirm:
- `|V'|` (fragment count) decreases as `U` increases, and is meaningfully smaller than `|V|` at every tested `U` (order-of-magnitude comparison against PUNCH's Table 1 — not an exact match, since we don't have their code).
- `fragment size: max=...` never exceeds the `U` passed on that run (the hard invariant from design spec section 6 — if this is ever violated on real data despite `filtering_invariant_test.cc` passing on synthetic graphs, that indicates a real-world edge case the synthetic tests didn't cover, and must be investigated before considering this module validated).
- Runtime is reasonable (PUNCH reports tiny-cut detection as a flat ~25-30s regardless of scale; these instances are far smaller than Europe/USA, so a much shorter runtime is expected — a multi-minute runtime on the NY-scale instance would be a signal to look at the `contract_two_edge_cuts`/`graph_contraction` `O(m log m)` merge step and the natural-cut parallel loop's `O(n)`-per-seed snapshot copy, both flagged as deferred optimizations above).

- [ ] **Step 4: Record findings and decide on the full-scale follow-up**

Report the `|V| -> |V'|` numbers per instance/`U` back in conversation. Full Europe/USA reproduction (Table 1 numbers, Alps/border visual check against Figures 7-8) is an explicit follow-up per design spec section 9, given the multi-GB download and runtime involved — decide with the user whether to proceed with it now or treat this module as validated at this scale.

---

## Plan self-review

**Spec coverage:** every section of `docs/superpowers/specs/2026-09-17-punch-filtering-design.md` maps to a task — section 2 (layout/build/extension point) -> Tasks 1, 17; section 3 (data model) -> Tasks 1, 2; section 4 (Part 1) -> Tasks 6-11 (with the bridge-tree correction recorded in Global Constraints); section 5 (Part 2) -> Tasks 12-15; section 6 (Part 3) -> Task 16; section 7 (CLI) -> Task 18; section 8 (testing) -> woven into every task's TDD steps; section 9 (validation) -> Task 19; section 10 (deferred scope) -> left untouched, as intended.

**Placeholder scan:** no "TBD"/"TODO" remain; the one previously-dead code fragment (in Task 16's draft) was removed during this review rather than left for the implementer to clean up.

**Type consistency:** `NodeID`/`EdgeID`/`NodeWeight`/`EdgeWeight` are defined once (Task 1) and used identically throughout; `ContractionResult{graph, mapping}` (Task 4) is threaded unchanged through Tasks 8-11; `FilteringResult{fragment_id, fragment_size, kept_edges}` (Task 16) is threaded unchanged through Task 17's `run_filtering_pipeline` and Task 18's CLI; `NaturalCutParams`/`NaturalCutScratch` (Task 13) are reused as-is by Tasks 14-15 without signature drift.

## Execution handoff

Plan complete and saved to `docs/superpowers/plans/2026-09-17-punch-filtering-implementation.md`. Two execution options:

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

Which approach?
