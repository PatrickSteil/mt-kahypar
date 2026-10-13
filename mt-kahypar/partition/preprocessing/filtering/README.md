# Natural-cut filtering

> **Note:** This summary was written by an AI assistant (Claude) and checked
> only briefly.

The idea: contract the graph into small *fragments* along natural cuts, then
let Mt-KaHyPar partition the graph, with the fragments guiding coarsening.

## Files

| File | What it does |
|---|---|
| `filter_graph.*`, `dimacs_io.*`, `metis_io.*` | Simple CSR graph (`FilterGraph`) and readers for DIMACS `.gr` and METIS files |
| `tiny_cut_detection.*` | PUNCH part 1: contract everything separated by tiny cuts (bridges, degree-2 chains, 2-edge cuts) whose small side weighs at most U |
| `natural_cut_detection.*`, `local_max_flow.*` | PUNCH part 2: BFS trees from random seeds, min cut between core and ring |
| `region_cut_detection.*` | **New idea:** split the graph into regions, one min cut per pair of adjacent regions |
| `fragment_assembly.*`, `filtering_pipeline.*` | PUNCH part 3: the kept cut edges define the fragments |

Tools: `tools/punch_filter.cc` (`PunchFilter`, classic PUNCH) and
`tools/region_filter.cc` (`RegionFilter`, region cuts plus evaluation and dump
options). Tests are in `tests/partition/preprocessing/filtering/`; region cuts
have none yet.

Mt-KaHyPar reads fragments with `--initial-partition` (see `multilevel.cpp`).

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DKAHYPAR_ENABLE_TESTING=ON
cmake --build build -j --target MtKaHyPar PunchFilter RegionFilter mtkahypar_tests
```

## Workflow

```bash
# 1. Fragments (.gr = DIMACS, anything else = METIS)
build/tools/PunchFilter --graph 4elt.metis --U 160 \
    --dump-partition 4elt.frag --dump-fragment-graph 4elt.frag.graph

# 2. Partition, starting from the fragments
build/mt-kahypar/application/MtKaHyPar -h 4elt.metis --input-file-format=metis \
    -k 32 -e 0.03 -o cut -m direct -t 8 --preset-type=default \
    --initial-partition 4elt.frag
```

- The fragment file has one fragment ID per vertex, one per line.
- If all IDs are below k, `--initial-partition` uses the file as a starting
  solution. Otherwise the IDs are fragments: coarsening only contracts within
  a fragment until it stalls. See also `--initial-partition-relax-on-stall`
  and `--initial-partition-vcycle-fragments`.
- PUNCH uses `U = Lmax / 3` with `Lmax = (1 + eps) * ceil(n / k)`.
  `RegionFilter` computes it from `--k` and `--eps`.
- The fragment graph (`--dump-fragment-graph`) is a METIS file. Partitioning
  it directly gives the same cut as on the original graph.

## Region cuts

```bash
build/tools/RegionFilter --graph 4elt.metis --k 32 --eps 0.03 \
    --terminals far --margin 1 --rounds 3 --dump-partition 4elt.region.frag
```

- `--eval-partition <file>` reports the share of a good partition's cut edges
  that lie between fragments (100% = that partition is still reachable).
- The settings to experiment with are in `RegionCutParams`
  (`region_cut_detection.h`).

## Tests

```bash
cd build/tests && ./mtkahypar_tests --gtest_filter='FilterGraph*:DimacsIo*:MetisIo*:AtomicUnionFind*:GraphContraction*:ParallelConnectivity*:CutSignatures*:Contract*Test.*:*TinyCut*:*MaxFlow*:MinCutSides*:*NaturalCut*:AssembleFragments*:FilteringPipeline*'
cd build/tests/interface && ./mtkahypar_interface_test --gtest_filter='*Fragments*:*InitialGraphPartition*'
```

If the tests crash with a TBB error about a "slot that should not exist",
limit the CPUs, e.g. `taskset -c 0-11 ./mtkahypar_tests ...`.
