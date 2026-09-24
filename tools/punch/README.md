# PUNCH-style partitioning of road networks

Two steps: `PunchFilter` finds natural cuts and contracts everything else into
fragments (PUNCH filtering phase), then `punch_evolutionary.py` partitions the
resulting fragment graph with Mt-KaHyPar (multistart + recombination, in place
of PUNCH's assembly phase).

```sh
# Lmax = (1 + eps) * ceil(n / k); USA, k = 32, eps = 0.03: Lmax = 770805
build/tools/PunchFilter --graph USA-road-t.USA.gr --U 385407 --cut-side sink \
    --dump-partition usa.frag --dump-fragment-graph usa.frag.graph
tools/punch/punch_evolutionary.py --mtkahypar build/mt-kahypar/application/MtKaHyPar \
    --fragment-graph usa.frag.graph --fragments usa.frag -k 32 -o usa.part32
```

The output is a partition of the original graph (one block ID per vertex); its
cut equals the cut of the fragment-graph partition.

## Filter settings

- **Unit edge weights** (default). The objective is the number of cut edges,
  as in the paper's unweighted setting. `--arc-weights` uses the DIMACS arc
  weights as capacities instead; on `USA-road-t` (travel times) this produced
  far fewer and much worse fragments.
- **`--cut-side sink`** keeps the minimum cut closest to the ring. With unit
  capacities, minimum cuts are rarely unique; the sink-side cut (what a
  first-phase push-relabel yields) is closest to the paper's fragment counts
  and gave better partitions than the source-side cut.
- **U**: the paper uses `U = Lmax / 3` for balanced partitions. `Lmax / 2`
  worked at least as well on USA.

Even so, we get 1.1-1.9x the paper's fragment counts (Table 14 on NY, BAY,
COL, FLA; the gap grows with U). The cause is unknown.

## Results (USA, eps = 0.03, 16 cores)

| k  | method                                   | avg cut | best | time            |
|----|------------------------------------------|---------|------|-----------------|
| 30 | Mt-KaHyPar highest_quality, no fragments | 2279    | 2172 | 39 s            |
| 32 | single run on fragment graph (U=Lmax/3)  | ~2070   | —    | 0.3 s           |
| 32 | best of 160 runs on fragment graph       | —       | 1902 | 36 s            |
| 32 | `punch_evolutionary.py` (U=Lmax/3)       | 1849    | 1837 | ~60 s + filter  |
| 32 | `punch_evolutionary.py` (U=Lmax/2)       | —       | 1829 | ~40 s + filter  |
| 32 | PUNCH (paper, Table 10)                  | 1883    | 1829 | ~115 s total    |

Filtering USA takes about 2.5 minutes. Refining the best fragment-graph
partition on the full graph (Mt-KaHyPar V-cycles) did not improve it.

What did not help much:
- Feeding fragments into Mt-KaHyPar directly (`--initial-partition` in
  fragment mode) gave about -5 to -10% vs. no fragments. Partitioning the
  fragment graph is better, because the fragment restriction is lifted before
  initial partitioning.
- Iterated V-cycles from the best multistart solutions (1909 -> 1891).

## Open ideas

- Evaluate more trials per filter setting (U, coverage, cut side) and at
  k = 2 ... 64 against the paper's Table 10.
- Speed up the filter; it now dominates the running time.
- Keep fragment restrictions active in initial partitioning, so the full
  Mt-KaHyPar pipeline benefits without a separate script.
- Move the recombination into Mt-KaHyPar (C++) instead of driving the binary.
- Explain the remaining fragment-count gap to the paper.
