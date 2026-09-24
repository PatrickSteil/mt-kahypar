#!/usr/bin/env python3
"""Evolutionary k-way partitioning of a PUNCH fragment graph with Mt-KaHyPar.

The fragment graph (written by `PunchFilter --dump-fragment-graph`) has one
vertex per fragment; any partition of it is a partition of the input graph with
the same cut. This script searches it in the style of PUNCH's combination
heuristic and KaFFPaE:

  1. Multistart: `--init-runs` independent Mt-KaHyPar runs; the best distinct
     balanced solutions form the population.
  2. Combine: overlay two parents (tournament selection), contract every cell
     (block in parent 1, block in parent 2) and partition the small cell graph
     with Mt-KaHyPar. Only edges cut by a parent can be cut, and both parents
     are feasible solutions of the cell graph. The best projected solution is
     then refined with a V-cycle on the fragment graph (--initial-partition).
  3. Mutation: a share of the children combine a parent with a fresh random
     solution instead of a second parent, to keep the population diverse.
  4. A child replaces the worst individual if it is balanced, better and not a
     duplicate cut value. Stops after `--generations` or `--patience`
     generations without improvement.

Runs are single-threaded and executed `--jobs` at a time.

Example:
  PunchFilter --graph USA.gr --U 385407 --cut-side sink \
      --dump-partition usa.frag --dump-fragment-graph usa.frag.graph
  punch_evolutionary.py --mtkahypar build/mt-kahypar/application/MtKaHyPar \
      --fragment-graph usa.frag.graph --fragments usa.frag -k 32 -o usa.part32
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor


class Graph:
    """METIS graph with vertex and edge weights (format code 11)."""

    def __init__(self, path):
        with open(path) as f:
            lines = [l for l in f.read().split("\n") if not l.startswith("%")]
        header = lines[0].split()
        self.n = int(header[0])
        fmt = header[2] if len(header) > 2 else "0"
        if fmt not in ("11", "011"):
            raise ValueError(f"{path}: expected METIS format 11 (vertex and edge weights), got {fmt}")
        self.weight, self.adj = [], []
        for line in lines[1:self.n + 1]:
            t = list(map(int, line.split()))
            self.weight.append(t[0])
            self.adj.append([(t[i] - 1, t[i + 1]) for i in range(1, len(t), 2)])
        self.total_weight = sum(self.weight)

    def cut(self, part):
        return sum(w for u in range(self.n) for v, w in self.adj[u] if u < v and part[u] != part[v])

    def is_balanced(self, part, k, eps):
        max_weight = int((1 + eps) * -(-self.total_weight // k))
        block_weight = [0] * k
        for u in range(self.n):
            block_weight[part[u]] += self.weight[u]
        return max(block_weight) <= max_weight


def read_partition(path):
    with open(path) as f:
        return [int(x) for x in f]


def write_partition(path, part):
    with open(path, "w") as f:
        f.write("\n".join(map(str, part)) + "\n")


class MtKaHyPar:
    def __init__(self, binary, k, eps, preset, workdir):
        self.binary, self.k, self.eps, self.preset, self.workdir = binary, k, eps, preset, workdir

    def partition(self, graph_path, seed, extra=()):
        """Runs Mt-KaHyPar with one thread; returns the partition (or None)."""
        out = tempfile.mkdtemp(dir=self.workdir)
        cmd = [self.binary, "-h", graph_path, "--input-file-format=metis", "-k", str(self.k),
               "-e", str(self.eps), "-o", "cut", "-m", "direct", f"--preset-type={self.preset}",
               "-t", "1", f"--seed={seed}", "--write-partition-file=true",
               f"--partition-output-folder={out}", *extra]
        result = subprocess.run(cmd, capture_output=True, text=True)
        files = os.listdir(out)
        part = read_partition(os.path.join(out, files[0])) if result.returncode == 0 and files else None
        shutil.rmtree(out)
        return part


class Evolution:
    def __init__(self, args, graph, mtk, workdir):
        self.args, self.graph, self.mtk, self.workdir = args, graph, mtk, workdir
        self.rng = random.Random(args.seed)
        self.next_seed = args.seed * 1_000_000
        self.population = []  # sorted list of (cut, partition)

    def seed(self):
        self.next_seed += 100
        return self.next_seed

    def evaluate(self, part):
        return self.graph.cut(part), self.graph.is_balanced(part, self.args.k, self.args.epsilon)

    def fresh(self, seed):
        part = self.mtk.partition(self.args.fragment_graph, seed)
        if part is None:
            return None
        cut, balanced = self.evaluate(part)
        return (cut, part) if balanced else None

    def combine(self, job):
        (_, a), (_, b), seed = job
        g = self.graph
        cells = {}
        cell = [cells.setdefault((a[u], b[u]), len(cells)) for u in range(g.n)]
        num_cells = len(cells)
        cell_weight = [0] * num_cells
        cell_adj = [dict() for _ in range(num_cells)]
        for u in range(g.n):
            cell_weight[cell[u]] += g.weight[u]
            for v, w in g.adj[u]:
                if cell[u] != cell[v]:
                    cell_adj[cell[u]][cell[v]] = cell_adj[cell[u]].get(cell[v], 0) + w
        cell_graph = os.path.join(self.workdir, f"cells{seed}.graph")
        with open(cell_graph, "w") as f:
            f.write(f"{num_cells} {sum(len(d) for d in cell_adj) // 2} 11\n")
            for c in range(num_cells):
                f.write(" ".join([str(cell_weight[c])] + [f"{v + 1} {w}" for v, w in cell_adj[c].items()]) + "\n")

        # Partition the cell graph a few times and keep the best projection
        # (balanced first); fall back to parent 1 if none succeeds.
        best = None
        for r in range(self.args.cell_runs):
            cell_part = self.mtk.partition(cell_graph, seed + r)
            if cell_part is None:
                continue
            child = [cell_part[cell[u]] for u in range(g.n)]
            cut, balanced = self.evaluate(child)
            if best is None or (balanced, -cut) > (best[1], -best[0]):
                best = (cut, balanced, child)
        os.remove(cell_graph)
        start = best[2] if best is not None else a

        # Refine the projected solution with a V-cycle on the fragment graph
        start_path = os.path.join(self.workdir, f"start{seed}.part")
        write_partition(start_path, start)
        child = self.mtk.partition(self.args.fragment_graph, seed + 99, [f"--initial-partition={start_path}"])
        os.remove(start_path)
        if child is None:
            return None
        cut, balanced = self.evaluate(child)
        return (cut, child) if balanced else None

    def insert(self, individual):
        if individual is None:
            return False
        cut, _ = individual
        full = len(self.population) >= self.args.population
        if any(cut == c for c, _ in self.population) or (full and cut >= self.population[-1][0]):
            return False
        if full:
            self.population.pop()
        self.population.append(individual)
        self.population.sort(key=lambda x: x[0])
        return True

    def tournament(self):
        x, y = self.rng.sample(self.population, 2)
        return min(x, y, key=lambda i: i[0])

    def run(self):
        args, start = self.args, time.time()
        with ThreadPoolExecutor(args.jobs) as ex:
            for individual in ex.map(self.fresh, [self.seed() for _ in range(args.init_runs)]):
                self.insert(individual)
            if len(self.population) < 2:
                sys.exit("error: multistart produced fewer than two balanced solutions")
            print(f"multistart: best={self.population[0][0]} population={len(self.population)} "
                  f"({time.time() - start:.0f}s)", flush=True)

            stale = 0
            for gen in range(args.generations):
                num_children = args.jobs
                num_mutations = int(round(num_children * args.mutation_share))
                fresh = list(ex.map(self.fresh, [self.seed() for _ in range(num_mutations)]))
                jobs = []
                for i in range(num_children):
                    p1 = self.tournament()
                    if i < num_mutations and fresh[i] is not None:
                        p2 = fresh[i]
                    else:
                        p2 = self.tournament()
                        while p2 is p1:
                            p2 = self.tournament()
                    jobs.append((p1, p2, self.seed()))
                best_before = self.population[0][0]
                for child in ex.map(self.combine, jobs):
                    self.insert(child)
                stale = 0 if self.population[0][0] < best_before else stale + 1
                print(f"generation {gen}: best={self.population[0][0]} worst={self.population[-1][0]} "
                      f"({time.time() - start:.0f}s)", flush=True)
                if stale >= args.patience:
                    break
        return self.population[0]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mtkahypar", required=True, help="path to the MtKaHyPar binary")
    p.add_argument("--fragment-graph", required=True, help="METIS fragment graph (PunchFilter --dump-fragment-graph)")
    p.add_argument("--fragments", help="fragment ID per original vertex (PunchFilter --dump-partition); "
                                       "if given, the output is projected to the original graph")
    p.add_argument("-k", type=int, required=True)
    p.add_argument("-e", "--epsilon", type=float, default=0.03)
    p.add_argument("-o", "--output", required=True, help="output partition file")
    p.add_argument("--preset", default="highest_quality")
    p.add_argument("--jobs", type=int, default=os.cpu_count(), help="parallel single-threaded runs")
    p.add_argument("--init-runs", type=int, default=32)
    p.add_argument("--population", type=int, default=16)
    p.add_argument("--generations", type=int, default=12)
    p.add_argument("--patience", type=int, default=4, help="stop after this many generations without improvement")
    p.add_argument("--mutation-share", type=float, default=0.25)
    p.add_argument("--cell-runs", type=int, default=4, help="partitioning runs per cell graph")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    graph = Graph(args.fragment_graph)
    workdir = tempfile.mkdtemp(prefix="punch_evo_")
    try:
        mtk = MtKaHyPar(args.mtkahypar, args.k, args.epsilon, args.preset, workdir)
        cut, part = Evolution(args, graph, mtk, workdir).run()
    finally:
        shutil.rmtree(workdir)

    if args.fragments:
        part = [part[f] for f in read_partition(args.fragments)]
    write_partition(args.output, part)
    print(f"cut={cut} written to {args.output}")


if __name__ == "__main__":
    main()
