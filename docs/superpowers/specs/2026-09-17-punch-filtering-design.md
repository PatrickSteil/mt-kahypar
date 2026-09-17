# PUNCH-style graph filtering module — design spec

Date: 2026-09-17
Status: approved for implementation planning

## 1. Goal & scope

Implement the filtering phase of PUNCH (Delling, Goldberg, Werneck et al., *Graph
Partitioning with Natural Cuts*; see `paper/Punch.pdf` §3.1–3.2) as a standalone
module operating on plain undirected weighted graphs, independent of Mt-KaHyPar's
hypergraph/partitioning internals.

Given a graph and a size cap `U`, the module partitions vertices into fragments
(contractible components) such that no fragment exceeds `U`, by finding and
preserving "natural cuts" (tiny cuts + heuristically-detected sparse local cuts)
and contracting everything else.

**Out of scope for this pass:** wiring this into Mt-KaHyPar's coarsening; the
assembly-phase heuristics described in PUNCH §4 (greedy contraction, local search,
multistart/combination) — those partition the *fragment graph* this module
produces, and are a separate future project.

**Decisions locked in during brainstorming:**
- Module lives inside the existing source tree at
  `mt-kahypar/partition/preprocessing/filtering/`, not a new top-level directory.
- Validation starts with small/medium DIMACS road instances (NY, BAY, COL); the
  full-scale Europe/USA Table-1-reproduction and Alps/border visual check are a
  follow-up once the pipeline is proven on small instances.
- The Pritchard–Thurimella 2-edge-cut algorithm is implemented from the actual
  paper (fetched and verified during design, see §4.3), not a simplified
  substitute.
- Delivery is phased with checkpoints: Part 1 (tiny-cut) implemented and tested
  end-to-end first, then Part 2 (natural-cut), then Part 3 (assembly) + validation.

## 2. Module layout & build integration

```
mt-kahypar/partition/preprocessing/filtering/
    filter_graph.h              # plain CSR graph struct, independent of ds::StaticGraph/ds::Graph
    dimacs_io.h/.cpp            # DIMACS .gr reader/writer
    union_find.h                # parallel-friendly union-find, standalone
    graph_contraction.h/.cpp    # shared contract_graph() primitive used by all tiny-cut passes
    parallel_connectivity.h/.cpp# parallel connected components / spanning tree (Shiloach-Vishkin style)
    tiny_cut_detection.h/.cpp   # Part 1: component tree, degree-2 chains, 2-edge-cuts (PT)
    dinic_max_flow.h/.cpp       # local s-t min-cut solver used by Part 2
    natural_cut_detection.h/.cpp# Part 2: BFS core/ring, parallel seed scheduling
    fragment_assembly.h/.cpp    # Part 3: final union-find + invariant helpers
    filtering_pipeline.h/.cpp   # orchestrates Part 1 -> Part 2 -> Part 3, documented extension point
    CMakeLists.txt

tools/
    punch_filter.cc             # standalone CLI

tests/partition/preprocessing/filtering/
    union_find_test.cc
    graph_contraction_test.cc
    tiny_cut_detection_test.cc
    natural_cut_detection_test.cc
    fragment_assembly_test.cc
    filtering_invariant_test.cc  # U-invariant, run on synthetic + small DIMACS graphs
```

Build integration follows the existing pattern used by
`mt-kahypar/partition/preprocessing/community_detection/`:
`target_sources(MtKaHyPar-Sources INTERFACE ${FilteringSources})` and the same
for `MtKaHyPar-ToolsSources`. `tools/punch_filter.cc` is registered in
`tools/CMakeLists.txt` following the existing single-file-tool pattern (e.g.
`verify_partition.cc`), linking against a `MtKaHyPar-BuildTools`-style target.
Tests are registered into the existing monolithic `mtkahypar_tests` gtest binary
via `target_sources(mtkahypar_tests PRIVATE ...)`, matching every other test
directory under `tests/partition/`.

**Independence boundary:** module headers may include TBB, the C++ standard
library, and generic low-level Mt-KaHyPar infrastructure that is not
hypergraph/partitioning-specific — specifically `parallel/atomic_wrapper.h`,
`parallel/stl/thread_locals.h`, `parallel/stl/scalable_vector.h`, and
`datastructures/array.h` (a NUMA-aware array wrapper, not a graph/hypergraph
type). They must never include `ds::StaticGraph`, `ds::StaticHypergraph`,
`ds::Graph`, `HypernodeID`, `Context`, or anything from `partition/` outside
this `filtering/` subdirectory. This is the boundary that keeps the module
standalone and testable without the rest of Mt-KaHyPar.

**Extension point for future integration:** `filtering_pipeline.h` exposes a
single entry point:

```cpp
// Runs the full PUNCH filtering pipeline (tiny-cut + natural-cut detection)
// on `graph`, returning a vertex -> fragment-id mapping plus fragment sizes
// and the set of kept (natural-cut) edges between fragments.
FilteringResult run_filtering_pipeline(const FilterGraph& graph, const FilteringParams& params);
```
`FilterGraph`/`FilteringResult`/`FilteringParams` are the only types a future
integration into Mt-KaHyPar's coarsening would need to construct from/consume
into hypergraph types — documented in a comment block at the top of
`filtering_pipeline.h` so the maintainer has a clear, single seam to hook into
later, without needing to understand the internals of Parts 1–3.

## 3. Data model

```cpp
using FilterNodeID = uint32_t;
using FilterEdgeID = uint32_t;
using NodeWeight    = uint64_t;   // s(v)
using EdgeWeight     = int64_t;   // w(e); DIMACS distances are integers

struct FilterGraph {
  ds::Array<FilterEdgeID> node_begin;   // CSR offsets, size n+1
  ds::Array<FilterNodeID> adj;          // size 2m, neighbor per half-edge
  ds::Array<FilterEdgeID> adj_edge;     // size 2m, half-edge -> canonical edge id
  ds::Array<NodeWeight>   node_weight;  // size n
  ds::Array<EdgeWeight>   edge_weight;  // size m, indexed by canonical edge id
};
```

Canonical edge ids (one id per undirected edge, referenced from both
half-edges) are required throughout: bridge/2-cut detection and the `keep[]`
array both key off edge id, not half-edge/direction.

`FilterNodeID`/`FilterEdgeID` default to 32-bit (Europe's ~18M vertices fits
comfortably); a 64-bit build variant can be added later the same way
`KAHYPAR_USE_64_BIT_IDS` gates the rest of the repo, but is not needed for this
pass.

**DIMACS `.gr` reader:** 9th DIMACS Challenge format — header line `p sp <n> <m>`,
arc lines `a <u> <v> <w>` (1-indexed, and each undirected edge appears as two
directed arc lines with the same weight in both directions). The reader
deduplicates the two arc lines per undirected edge into one canonical edge,
verifying the two directions agree on weight (fail loudly if they don't — a
malformed-input signal, not a case to silently paper over). Vertex weights
default to 1 (`s(v) = 1`), matching PUNCH's usage of vertex count as size on
road networks; `.co` coordinate files are not needed by the graph/filtering
logic itself, only later for the optional Alps/border plotting sanity check.

## 4. Part 1 — tiny-cut detection

### 4.1 Shared primitive: `contract_graph`

All three tiny-cut sub-passes contract vertices. Rather than three bespoke
contraction routines, each pass produces a union-find over the *current*
graph's vertices, and one shared routine rebuilds the next stage's `FilterGraph`
from it: sum vertex weights per equivalence class, merge parallel edges by
summing weights, drop self-loops created by contraction. Parallelized by
bucketing edges by `(rep(u), rep(v))` via TBB parallel sort/hash-aggregation —
the same well-known technique used for Louvain graph coarsening in
`ds::Graph`, reimplemented standalone here (no dependency on that code).

### 4.2 Pass 1 — component tree

Parallel connected components via Shiloach–Vishkin-style pointer-jumping over
an atomic union-find (a classic parallel-connectivity technique, independent of
PT — cite Shiloach & Vishkin 1982, and note Shun/Blelloch-style approaches as an
alternative if benchmarking later shows a need). The same pass yields parent
pointers usable as the component tree. Root = largest component. Top-down:
contract any subtree with total size ≤ `U`; then optionally fold the
contracted vertex into its parent if the subtree's size ≤ τ (default 5) and the
merged size ≤ `U`. Implemented as one union-find pass over the tree, then
`contract_graph`.

### 4.3 Pass 2 — degree-2 chains

One parallel pass marks degree-2 vertices. Chains are maximal runs of marked
vertices; a degree-2 vertex self-elects as a "chain start" if at least one
neighbor is not degree-2 (list-ranking-friendly). Each chain is then walked and
unioned independently and in parallel, stopping the union once total size
would exceed `U`. Then `contract_graph`.

### 4.4 Pass 3 — 2-edge-cuts (Pritchard–Thurimella)

**Source:** Pritchard & Thurimella, *Fast Computation of Small Cuts via Cycle
Space Sampling*, ACM Transactions on Algorithms 7(4), 2011 (arXiv:cs/0702113).

**Algorithm** (verified against brute force on 300 random graphs plus the
simple-cycle edge case during design — see §4.5):

1. Build a spanning tree `T` of the graph produced by pass 2 (reusing the
   connectivity machinery from §4.2, run again on this smaller graph).
2. Assign each **non-tree** edge `x` a uniform random 128-bit label `r_x`.
3. For each tree edge `e`, compute `agg(e)` = XOR of `r_x` over every non-tree
   edge `x` whose fundamental cycle covers `e`. Standard technique: for each
   non-tree edge `x = (u, v)`, XOR `r_x` into `diff[u]` and `diff[v]`; one
   post-order accumulation pass over `T` turns `diff` into `agg` per tree edge.
   `O(n + m)` total, no need for PT's general small-`k` machinery since only
   `k = 2` is required here.
4. **Signature** of every edge: tree edge `e` -> `agg(e)`; **non-tree edge `x`
   -> `r_x` itself.** (This is the critical detail: bucketing only tree-edge
   aggregates misses every 2-cut that includes a non-tree edge — e.g. any two
   edges of a simple cycle, one of the most common structures on road
   networks. Including non-tree edges' own raw labels in the same bucketing
   catches these correctly; see the proof sketch in §4.5.)
5. **Bridges** = tree edges with `agg(e) = 0`. Non-tree edges can never be
   bridges (the spanning tree alone already connects the graph without them).
6. Bucket all non-bridge edges by signature; any bucket with ≥ 2 members is a
   genuine 2-edge-cut equivalence class, matching PUNCH's equivalence relation
   `P` (`(e,f) ∈ P ⟺ e=f, or e,f form a 2-cut and neither is a 1-cut`) exactly.
7. Monte Carlo with one-sided error; at 128 bits the collision probability is
   negligible at road-network scale, but each candidate class is still cheaply
   verified (check the two representative edges' endpoints are actually
   disconnected in `G` minus those two edges via a bounded local traversal)
   before being trusted — an explicit, documented design decision, analogous
   to the max-flow-choice note in §5.
8. For each verified class `S`, compute connected components of `G_S = (V, E∖S)`
   and union any component with size ≤ `U`, using the two-components-at-a-time
   traversal trick from the paper to bound total work to twice the size of the
   smaller components across all classes. Then `contract_graph`.

Output of Part 1: a smaller `FilterGraph` plus a mapping
`orig_vertex -> part1_vertex` (composed across all three passes).

### 4.5 Correctness argument for step 4 (recorded for the implementer)

Edge cuts of size 2 correspond exactly to nonzero elements of size 2 in the
graph's cut space (cographic matroid), which — since the cut space is
orthogonal to the cycle space — can be characterized via the fundamental cycle
basis w.r.t. `T`: `{e, f}` is a cut iff every fundamental cycle `C_x`
intersects `{e, f}` an even number of times (0 or 2), for every non-tree edge
`x`. Since each fundamental cycle `C_x` contains exactly one non-tree edge
(itself) plus tree edges on the path between `x`'s endpoints, this splits into
three cases:
- **Both `e, f` tree edges:** condition becomes "`e` and `f` are covered by
  exactly the same set of non-tree edges' fundamental cycles" — captured by
  `agg(e) = agg(f)` with high probability (random-label XOR argument: if the
  covering sets differ, their symmetric difference is a nonzero linear
  combination of independent random values, hence nonzero whp).
- **One of `e, f` non-tree (say `e = x`):** since `x` only appears in `C_x`,
  the condition reduces to "`f` lies on `x`'s own fundamental cycle, and `f` is
  covered by no *other* non-tree edge's cycle" — i.e. `agg(f) = r_x` exactly.
  This is why non-tree edges' signatures must be their raw label `r_x`: it
  lets them bucket-match against a tree edge whose aggregate happens to equal
  that single term.
- **Both `e, f` non-tree, distinct:** reduces to requiring `f ∈ C_e`, which is
  impossible since `C_e` contains only tree edges besides `e` itself. So two
  distinct non-tree edges never form a 2-cut together — consistent with the
  fact that the spanning tree alone (untouched by removing only non-tree
  edges) keeps the graph connected.

This was validated computationally during design: brute-force enumeration of
all edge pairs on 300 random small connected graphs, plus an explicit 6-cycle
graph (which brute-forces to a single equivalence class containing all 6
edges, since removing any two edges of a cycle disconnects it), matched the
signature-based classification exactly.

## 5. Part 2 — natural-cut detection

Parameters: `α` (default 1.0), `f` (default 10), `C` (coverage, default 2).

**Per-seed local subproblem.** For center `v`: BFS-grow tree `T` from `v` until
`s(T) ≥ αU`. `core` = vertices added before `s(T)` reached `αU/f`. `ring` =
neighbors of `T` in `V ∖ T`. Build a small local flow network in **thread-local
scratch buffers** (never mutating the real graph): contract `core -> s`,
`ring -> t`; edges between two core (or two ring) vertices collapse away;
edges between a `T`-interior vertex and `core`/`ring` retain their weight as
capacity. Run Dinic's max-flow `s -> t`; read the min-cut edge set off the
final residual graph (edges crossing from the side reachable from `s` to the
side that isn't).

**Max-flow choice (explicit design decision):** Dinic's algorithm
(`O(V²E)` worst case) rather than push-relabel with global relabeling. PUNCH
uses push-relabel because it amortizes well when reused across many phases of
a single large flow computation; here, each of the many local subproblems is
small (bounded by BFS growth to size ~`αU`) and independent, so push-relabel's
higher setup/global-relabeling overhead is unlikely to pay off. Documented in
a code comment as swappable pending benchmarking, per the task's request.

**Parallel seed scheduling (new design — PUNCH's is explicitly sequential "for
simplicity", so this needs fresh design, not a port):**
- For `c` in `1..C`: reset `std::atomic<uint8_t> covered[n]` to false, run one
  fully parallel sweep. (Each sweep is an independent random cover, matching
  the paper's semantics of repeating the whole sweep `C` times, while keeping
  each sweep itself maximally parallel.)
- Sweep: TBB parallel work over a **pre-shuffled** vertex order (shuffling up
  front approximates "pick uniformly at random among uncovered" without a
  shared random cursor/contention point). Each task visits its vertex `v`:
  `compare_exchange(covered[v], false, true)`; on failure, skip (already
  claimed or covered by another seed's growth); on success, `v` becomes a seed
  and the task runs BFS growth -> core/ring -> local Dinic's using
  `tls_enumerable_thread_specific` scratch buffers (matching the existing
  pattern in `parallel/stl/thread_locals.h`) to avoid per-task allocation.
- During BFS growth, every visited vertex is marked `covered[v] = true` via a
  plain atomic store (not CAS — idempotent; a benign race here causes at most
  minor redundant work if two threads' growth regions briefly overlap, never
  an incorrectness, since overlap only means a vertex's core/ring membership
  is computed slightly redundantly by two seeds, not that either computation
  is wrong).
- `keep[]` (kept natural-cut edges): `std::atomic<uint8_t>` per edge, set with
  a plain relaxed atomic store — a monotonic boolean flag needs no
  compare-exchange, simplifying over a literal CAS.
- End of sweep: since a handful of vertices can race near the boundary of the
  shuffled scan, a short follow-up pass repeats over any vertices still
  uncovered until none remain (bounded, since coverage is monotonic within a
  sweep).

Output of Part 2: the `keep[]` edge-flag array over Part 1's output graph.

## 6. Part 3 — fragment assembly

Union-find over Part 1's output graph, unioning the endpoints of every edge
**not** flagged in `keep[]`. Resulting components are fragments. Composed with
Part 1's `orig_vertex -> part1_vertex` mapping to produce the final
`orig_vertex -> fragment_id` array, fragment sizes (sum of `s(v)` per
fragment), and the kept-edge set translated back to original vertex ids (for
the later visual sanity check).

**Hard invariant:** no fragment may exceed `U` (guaranteed by `α ≤ 1` per the
paper's own argument). This must be an executable assertion in the test suite,
not just documented — `filtering_invariant_test.cc` runs the full pipeline on
synthetic random graphs and (once available) the small DIMACS instances, and
asserts `max(fragment_sizes) <= U` for every run.

## 7. CLI

`tools/punch_filter.cc`:
```
punch_filter --graph <dimacs.gr> --U <size>
             [--alpha 1.0] [--f 10] [--coverage 2] [--tau 5] [--threads N]
             [--dump-kept-edges <path>]
```
Prints: `|V| -> |V'|` reduction after Part 1, after Part 2, per-stage timings,
fragment size histogram (min/median/max, count > 0.9·U as a near-cap sanity
signal). `--dump-kept-edges` writes `(orig_u, orig_v)` pairs of the final kept
edge set for later plotting against the DIMACS `.co` coordinates.

## 8. Testing

- `union_find_test.cc`, `graph_contraction_test.cc`: correctness of the shared
  primitives on small synthetic graphs.
- `tiny_cut_detection_test.cc`: synthetic cases for each sub-pass — a small
  component hanging off a 1-cut (pass 1), a degree-2 path (pass 2), and for
  pass 3: a bridge (excluded), a simple cycle (all edges one equivalence
  class — the Case-B scenario that the fork's original design missed), and a
  "theta graph" (two 2-cut vertices with three parallel paths between them,
  to confirm classes don't over-merge across genuinely different cut pairs).
- `natural_cut_detection_test.cc`: a small hand-checkable graph with a known
  min-cut value, confirming the local Dinic's + core/ring construction finds
  it.
- `fragment_assembly_test.cc`, `filtering_invariant_test.cc`: as in §6.

## 9. Validation plan

1. Download `USA-road-d.{NY,BAY,COL}.gr.gz` from
   `http://www.diag.uniroma1.it/challenge9/data/USA-road-d/` (confirmed live
   during design).
2. Run `punch_filter` at a few `U` values per instance, comparing `|V'|`
   against the right order of magnitude relative to PUNCH's Table 1 (not an
   exact match, since we don't have their code).
3. Assert the `U`-invariant holds on every run (see §6).
4. Full-scale Europe/USA reproduction (Table 1 numbers, Alps/border visual
   check against Figures 7–8) is an explicit follow-up once the above passes,
   given the multi-GB download and runtime involved.

## 10. Out of scope / explicitly deferred

- Wiring `run_filtering_pipeline` into Mt-KaHyPar's actual coarsening
  (`partition/coarsening/`) — the extension point in §2 exists for this, but
  no caller is added in this pass.
- PUNCH §4's assembly-phase heuristics (greedy contraction, local search,
  multistart/combination) — a separate future project consuming this module's
  output.
- 64-bit vertex/edge IDs, full Europe/USA-scale validation, and Alps/border
  visual plotting — deferred as noted in §9.
