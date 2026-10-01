// Region-based natural cuts (see region_cut_detection.h): tiny cuts, then one
// min cut per pair of adjacent regions, then fragment assembly.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/dimacs_io.h"
#include "mt-kahypar/partition/preprocessing/filtering/fragment_assembly.h"
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"
#include "mt-kahypar/partition/preprocessing/filtering/metis_io.h"
#include "mt-kahypar/partition/preprocessing/filtering/region_cut_detection.h"
#include "mt-kahypar/partition/preprocessing/filtering/tiny_cut_detection.h"

using namespace mt_kahypar::filtering;

namespace {
void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " --graph <dimacs.gr|metis.graph> (--k <blocks> [--eps 0.03] | --U <size>) "
            << "[--region-factor 0.5] [--region-factors f1,f2,...] [--regions <num>] [--lloyd 3] [--terminal-fraction 0.3] [--margin 1] [--pin-ring] [--closed [--no-fallback]] [--terminals far|center|ball] [--ball-fraction 0.3] [--rounds 1] [--min-rounds 1] [--centers random|spread] [--spread-fraction 0.3] "
            << "[--cut-side source|sink|both] [--tau 5] [--no-tiny-cuts] [--seed 0] "
            << "[--dump-partition <path>] [--dump-regions <path>] [--dump-cut-edges <path>] [--dump-pairs <region|auto> <prefix>] [--dump-fragment-graph <path>] [--eval-partition <path>]... [--verbose]\n"
            << "  U = (1 + eps) * ceil(n / k) / 3 (PUNCH's U = Lmax / 3); regions = ceil(n / (region-factor * U)).\n"
            << "  --region-factors: round i uses factor f_(i mod count), with margin * sqrt(f_i / region-factor).\n"
            << "  --min-rounds c keeps only edges cut in at least c rounds; --eval-partition prints, for every c,\n"
            << "  the fragments and the share of the given partition's cut edges that lie between fragments.\n";
}

bool ends_with(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

double seconds_since(std::chrono::steady_clock::time_point t) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}
}  // namespace

int main(int argc, char** argv) {
  std::string graph_path, partition_path, regions_path, fragment_graph_path, cut_edges_path;
  std::string pairs_region, pairs_prefix;
  std::vector<std::string> eval_paths;
  NodeWeight U = 0, tau = 5;
  uint64_t k = 0;
  double eps = 0.03, region_factor = 0.5;
  NodeID num_regions = 0;
  std::vector<double> region_factors;
  int rounds = 1, min_rounds = 1;
  bool tiny_cuts = true, verbose = false, fallback = true;
  RegionCutParams params{0};

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        std::exit(1);
      }
      return argv[++i];
    };
    if (arg == "--graph") graph_path = next();
    else if (arg == "--k") k = std::stoull(next());
    else if (arg == "--eps") eps = std::stod(next());
    else if (arg == "--U") U = std::stoull(next());
    else if (arg == "--region-factor") region_factor = std::stod(next());
    else if (arg == "--regions") num_regions = std::stoul(next());
    else if (arg == "--region-factors") {
      std::stringstream list(next());
      for (std::string f; std::getline(list, f, ',');) region_factors.push_back(std::stod(f));
    }
    else if (arg == "--lloyd") params.lloyd_iterations = std::stoi(next());
    else if (arg == "--terminal-fraction") params.terminal_fraction = std::stod(next());
    else if (arg == "--margin") params.margin = std::stoi(next());
    else if (arg == "--pin-ring") params.pin_ring = true;
    else if (arg == "--closed") params.closed = true;
    else if (arg == "--no-fallback") fallback = false;
    else if (arg == "--terminals") {
      const std::string choice = next();
      if (choice == "far") params.terminals = TerminalChoice::Far;
      else if (choice == "center") params.terminals = TerminalChoice::Center;
      else if (choice == "ball") params.terminals = TerminalChoice::Ball;
      else {
        print_usage(argv[0]);
        return 1;
      }
    }
    else if (arg == "--ball-fraction") params.ball_fraction = std::stod(next());
    else if (arg == "--rounds") rounds = std::stoi(next());
    else if (arg == "--min-rounds") min_rounds = std::stoi(next());
    else if (arg == "--centers") {
      const std::string choice = next();
      if (choice == "random") params.centers = CenterChoice::Random;
      else if (choice == "spread") params.centers = CenterChoice::Spread;
      else {
        print_usage(argv[0]);
        return 1;
      }
    }
    else if (arg == "--spread-fraction") params.spread_fraction = std::stod(next());
    else if (arg == "--eval-partition") eval_paths.push_back(next());
    else if (arg == "--cut-side") {
      const std::string side = next();
      if (side == "source") params.cut_side = CutSide::Source;
      else if (side == "sink") params.cut_side = CutSide::Sink;
      else if (side == "both") params.cut_side = CutSide::Both;
      else {
        print_usage(argv[0]);
        return 1;
      }
    } else if (arg == "--tau") tau = std::stoull(next());
    else if (arg == "--no-tiny-cuts") tiny_cuts = false;
    else if (arg == "--seed") params.seed = std::stoull(next());
    else if (arg == "--dump-partition") partition_path = next();
    else if (arg == "--dump-regions") regions_path = next();
    else if (arg == "--dump-cut-edges") cut_edges_path = next();
    else if (arg == "--dump-pairs") {
      pairs_region = next();
      pairs_prefix = next();
    }
    else if (arg == "--dump-fragment-graph") fragment_graph_path = next();
    else if (arg == "--verbose") verbose = true;
    else {
      print_usage(argv[0]);
      return 1;
    }
  }
  if (graph_path.empty() || (k == 0 && U == 0)) {
    print_usage(argv[0]);
    return 1;
  }

  FilterGraph graph = ends_with(graph_path, ".gr") ? read_dimacs_graph(graph_path) : read_metis_graph(graph_path);
  const size_t n = graph.numNodes();
  NodeWeight total_weight = 0;
  for (size_t v = 0; v < n; ++v) total_weight += graph.node_weight[v];
  if (U == 0) {
    const double lmax = (1.0 + eps) * std::ceil(static_cast<double>(total_weight) / k);
    U = static_cast<NodeWeight>(lmax / 3);
  }
  if (num_regions == 0) {
    num_regions = static_cast<NodeID>(std::ceil(total_weight / (region_factor * U)));
  }
  params.num_regions = num_regions;
  params.max_enclosed_weight = U;
  std::cout << "graph: |V| = " << n << ", |E| = " << graph.numEdges() << "\n"
            << "U = " << U << ", target regions = " << num_regions << "\n";

  const auto start = std::chrono::steady_clock::now();
  ContractionResult contracted;
  if (tiny_cuts) {
    contracted = run_tiny_cut_detection(graph, TinyCutParams{U, tau}, verbose);
  } else {
    contracted.mapping.resize(n);
    std::iota(contracted.mapping.begin(), contracted.mapping.end(), 0);
    contracted.graph = contract_graph_by_mapping(graph, contracted.mapping, n);
  }
  const double tiny_seconds = seconds_since(start);
  std::cout << "tiny cuts: " << contracted.graph.numNodes() << " vertices, " << contracted.graph.numEdges()
            << " edges (" << tiny_seconds << "s)\n";

  const auto region_start = std::chrono::steady_clock::now();
  // Each round uses a new set of regions (different seed; with --centers
  // spread, centered away from the centers of all earlier rounds). The kept
  // edges are those cut in at least min_rounds rounds (default: the union).
  // Region statistics are from round 0.
  std::vector<NodeID> previous_centers;
  params.avoid_centers = &previous_centers;
  // With --region-factors, the margin scales with the region diameter.
  auto params_for_round = [&](int round) {
    RegionCutParams p = params;
    p.seed = params.seed + round;
    if (!region_factors.empty()) {
      const double f = region_factors[round % region_factors.size()];
      p.num_regions = static_cast<NodeID>(std::ceil(total_weight / (f * U)));
      p.margin = static_cast<int>(std::lround(params.margin * std::sqrt(f / region_factor)));
    }
    if (verbose) std::cerr << "[round " << round << "] regions " << p.num_regions << ", margin " << p.margin << "\n";
    return p;
  };
  RegionCutStats stats;
  std::vector<char> keep = run_region_cut_detection(contracted.graph, params_for_round(0), &stats, verbose);
  std::vector<uint16_t> times_cut(keep.begin(), keep.end());
  previous_centers = stats.regions.centers;
  // For --dump-regions / --dump-cut-edges: regions of every round, and per
  // contracted edge a bitmask of the rounds whose cuts contain it.
  std::vector<std::vector<NodeID>> round_regions;
  std::vector<uint32_t> round_mask;
  const bool record_rounds = !regions_path.empty() || !cut_edges_path.empty();
  auto record_round = [&](int round, const std::vector<char>& round_keep, const Regions& regions) {
    if (!record_rounds) return;
    round_regions.push_back(regions.region_of);
    round_mask.resize(round_keep.size(), 0);
    for (size_t e = 0; e < round_keep.size(); ++e)
      if (round_keep[e]) round_mask[e] |= 1u << std::min(round, 30);
  };
  record_round(0, keep, stats.regions);
  int64_t flow_sum = 0;
  size_t num_cuts = stats.region_edges.size();
  for (int64_t f : stats.flow_value) flow_sum += f;
  for (int round = 1; round < rounds; ++round) {
    const RegionCutParams round_params = params_for_round(round);
    RegionCutStats round_stats;
    std::vector<char> round_keep = run_region_cut_detection(contracted.graph, round_params, &round_stats, verbose);
    record_round(round, round_keep, round_stats.regions);
    for (size_t e = 0; e < keep.size(); ++e) times_cut[e] += round_keep[e];
    previous_centers.insert(previous_centers.end(), round_stats.regions.centers.begin(),
                            round_stats.regions.centers.end());
    for (int64_t f : round_stats.flow_value) flow_sum += f;
    num_cuts += round_stats.flow_value.size();
    for (size_t v = 0; v < stats.covered.size(); ++v) stats.covered[v] |= round_stats.covered[v];
  }
  for (size_t e = 0; e < keep.size(); ++e) keep[e] = times_cut[e] >= min_rounds;

  for (const std::string& eval_path : eval_paths) {
    // Proxy for fragment quality: which share of a reference partition's cut
    // edges lie between fragments (100% = the partition is still reachable
    // on the fragment graph), for every min_rounds threshold.
    std::vector<NodeID> reference(n);
    std::ifstream in(eval_path);
    for (size_t v = 0; v < n; ++v) {
      if (!(in >> reference[v])) {
        std::cerr << "error: " << eval_path << " has fewer than " << n << " entries\n";
        return 1;
      }
    }
    for (int c = 1; c <= rounds; ++c) {
      std::vector<char> keep_c(keep.size());
      for (size_t e = 0; e < keep.size(); ++e) keep_c[e] = times_cut[e] >= c;
      const FilteringResult r = assemble_fragments(contracted.graph, keep_c, contracted.mapping);
      size_t reference_cut = 0, covered = 0;
      for (NodeID a = 0; a < n; ++a) {
        for (EdgeID pos = graph.node_begin[a]; pos < graph.node_begin[a + 1]; ++pos) {
          const NodeID b = graph.adj[pos];
          if (b <= a || reference[a] == reference[b]) continue;
          reference_cut += graph.edge_weight[graph.adj_edge[pos]];
          if (r.fragment_id[a] != r.fragment_id[b]) covered += graph.edge_weight[graph.adj_edge[pos]];
        }
      }
      std::cout << "eval " << eval_path << " min-rounds " << c << ": fragments " << r.fragment_size.size() << ", reference cut covered "
                << covered << "/" << reference_cut << " (" << 100.0 * covered / std::max<size_t>(1, reference_cut)
                << "%)\n";
    }
  }

  size_t num_covered = 0, fallback_cuts = 0;
  if (params.closed) {
    // Coverage rule: every vertex not enclosed by a bounded closed region cut
    // becomes a PUNCH natural-cut seed, so every fragment weighs at most U.
    for (char c : stats.covered) num_covered += c;
    if (fallback) {
      std::vector<NodeID> order(contracted.graph.numNodes());
      std::iota(order.begin(), order.end(), 0);
      std::mt19937_64 rng(params.seed);
      std::shuffle(order.begin(), order.end(), rng);
      NaturalCutParams punch_params{U, 1.0, 10.0, 1, params.cut_side};
      NaturalCutScratch punch_scratch;
      for (NodeID v : order) {
        if (stats.covered[v]) continue;
        for (EdgeID e : compute_natural_cut(contracted.graph, v, punch_params, punch_scratch, stats.covered))
          keep[e] = 1;
        ++fallback_cuts;
      }
    }
  }
  const double region_seconds = seconds_since(region_start);
  FilteringResult result = assemble_fragments(contracted.graph, keep, contracted.mapping);
  const double total_seconds = seconds_since(start);

  // Region statistics (region weight = number of original vertices).
  const NodeID k_regions = stats.regions.num_regions;
  std::vector<NodeWeight> region_weight(k_regions, 0);
  for (NodeID v = 0; v < contracted.graph.numNodes(); ++v)
    region_weight[stats.regions.region_of[v]] += contracted.graph.node_weight[v];
  const auto [rmin, rmax] = std::minmax_element(region_weight.begin(), region_weight.end());
  // Dangling cut edges: marked, but both endpoints in the same fragment.
  std::vector<NodeID> fragment_of(contracted.graph.numNodes());
  for (size_t v = 0; v < n; ++v) fragment_of[contracted.mapping[v]] = result.fragment_id[v];
  size_t marked = 0, dangling = 0;
  for (NodeID u = 0; u < contracted.graph.numNodes(); ++u) {
    for (EdgeID pos = contracted.graph.node_begin[u]; pos < contracted.graph.node_begin[u + 1]; ++pos) {
      const NodeID v = contracted.graph.adj[pos];
      if (v <= u || !keep[contracted.graph.adj_edge[pos]]) continue;
      ++marked;
      dangling += fragment_of[u] == fragment_of[v];
    }
  }
  std::cout << "regions: " << k_regions << " (size min=" << *rmin << " avg=" << total_weight / k_regions
            << " max=" << *rmax << "), region edges: " << stats.region_edges.size()
            << ", avg region degree: " << 2.0 * stats.region_edges.size() / k_regions << "\n"
            << "min cuts: " << num_cuts << " in " << rounds << " round(s), avg value "
            << static_cast<double>(flow_sum) / std::max<size_t>(1, num_cuts)
            << ", marked edges (union) " << marked << ", dangling " << dangling << " (" << region_seconds << "s)\n";
  if (params.closed) {
    std::cout << "coverage: " << num_covered << "/" << contracted.graph.numNodes()
              << " vertices enclosed by bounded region cuts, " << fallback_cuts << " PUNCH fallback cuts\n";
  }

  // Fragment statistics and total cut (original edges between fragments).
  const size_t num_fragments = result.fragment_size.size();
  const auto [fmin, fmax] = std::minmax_element(result.fragment_size.begin(), result.fragment_size.end());
  size_t over_U = 0;
  for (NodeWeight s : result.fragment_size) over_U += s > U;
  std::vector<NodeID> mapping(result.fragment_id.begin(), result.fragment_id.end());
  FilterGraph fragment_graph = contract_graph_by_mapping(graph, mapping, num_fragments);
  int64_t total_cut = 0;
  for (size_t e = 0; e < fragment_graph.numEdges(); ++e) total_cut += fragment_graph.edge_weight[e];
  std::cout << "fragments: " << num_fragments << " (size min=" << *fmin << " avg=" << total_weight / num_fragments
            << " max=" << *fmax << ", >U: " << over_U << ")\n"
            << "fragment graph: " << fragment_graph.numEdges() << " edges, total cut = " << total_cut << "\n"
            << "time: " << total_seconds << "s\n";

  if (!partition_path.empty()) {
    std::ofstream out(partition_path);
    for (NodeID fragment : result.fragment_id) out << fragment << "\n";
    std::cout << "wrote fragment id per vertex to " << partition_path << "\n";
  }
  if (!regions_path.empty()) {
    // One line per original vertex: its region in round 0, 1, ...
    std::ofstream out(regions_path);
    for (NodeID v : contracted.mapping) {
      for (size_t r = 0; r < round_regions.size(); ++r) out << (r ? " " : "") << round_regions[r][v];
      out << "\n";
    }
    std::cout << "wrote regions of " << round_regions.size() << " round(s) to " << regions_path << "\n";
  }
  if (!cut_edges_path.empty()) {
    // Every original edge whose contracted edge is marked: "u v rounds
    // boundary" (0-based ids; rounds = bitmask of the rounds that cut it, bit
    // 31 = PUNCH fallback; boundary = 1 if it lies between two fragments,
    // 0 if it dangles inside a fragment).
    std::ofstream out(cut_edges_path);
    size_t written = 0;
    for (NodeID a = 0; a < n; ++a) {
      for (EdgeID pos = graph.node_begin[a]; pos < graph.node_begin[a + 1]; ++pos) {
        const NodeID b = graph.adj[pos];
        const NodeID u = contracted.mapping[a], v = contracted.mapping[b];
        if (b <= a || u == v) continue;
        for (EdgeID cpos = contracted.graph.node_begin[u]; cpos < contracted.graph.node_begin[u + 1]; ++cpos) {
          if (contracted.graph.adj[cpos] != v) continue;
          const EdgeID e = contracted.graph.adj_edge[cpos];
          if (keep[e]) {
            const uint32_t mask = round_mask[e] ? round_mask[e] : (1u << 31);
            out << a << " " << b << " " << mask << " " << (result.fragment_id[a] != result.fragment_id[b]) << "\n";
            ++written;
          }
          break;
        }
      }
    }
    std::cout << "wrote " << written << " cut edges to " << cut_edges_path << "\n";
  }
  if (!pairs_prefix.empty()) {
    // Every subproblem of one round-0 region a (with each neighbor b):
    // <prefix>.<a>_<b>.txt with "pair a b flow", then per original vertex of
    // the flow subgraph "v id role source_side" (role: 0 free in a u b,
    // 1 source S1, 2 sink S2, 3 margin, 4 ring->s, 5 ring->t), then the
    // original cut edges "e u v".
    const Regions& regions = stats.regions;
    std::vector<std::vector<NodeID>> neighbors(regions.num_regions);
    for (auto [a, b] : stats.region_edges) {
      neighbors[a].push_back(b);
      neighbors[b].push_back(a);
    }
    NodeID a = 0;
    if (pairs_region == "auto") {
      for (NodeID r = 0; r < regions.num_regions; ++r)
        if (neighbors[r].size() > neighbors[a].size()) a = r;
    } else {
      a = std::stoul(pairs_region);
    }
    const size_t cn = contracted.graph.numNodes();
    std::vector<int8_t> role(cn, -1);
    std::vector<char> side(cn, 0), is_cut(contracted.graph.numEdges(), 0);
    for (NodeID b : neighbors[a]) {
      PairSubproblem sub = solve_region_pair(contracted.graph, regions, a, b, params);
      for (size_t i = 0; i < sub.vertices.size(); ++i) {
        role[sub.vertices[i]] = static_cast<int8_t>(sub.role[i]);
        side[sub.vertices[i]] = sub.source_side[i];
      }
      for (EdgeID e : sub.cut_edges) is_cut[e] = 1;
      const std::string path = pairs_prefix + "." + std::to_string(a) + "_" + std::to_string(b) + ".txt";
      std::ofstream out(path);
      out << "pair " << a << " " << b << " " << sub.flow << "\n";
      for (NodeID v = 0; v < n; ++v) {
        const NodeID u = contracted.mapping[v];
        if (role[u] >= 0) out << "v " << v << " " << int(role[u]) << " " << int(side[u]) << "\n";
      }
      for (NodeID x = 0; x < n; ++x) {
        for (EdgeID pos = graph.node_begin[x]; pos < graph.node_begin[x + 1]; ++pos) {
          const NodeID y = graph.adj[pos];
          const NodeID u = contracted.mapping[x], w = contracted.mapping[y];
          if (y <= x || u == w || role[u] < 0 || role[w] < 0) continue;
          for (EdgeID cpos = contracted.graph.node_begin[u]; cpos < contracted.graph.node_begin[u + 1]; ++cpos) {
            if (contracted.graph.adj[cpos] == w) {
              if (is_cut[contracted.graph.adj_edge[cpos]]) out << "e " << x << " " << y << "\n";
              break;
            }
          }
        }
      }
      for (NodeID v : sub.vertices) role[v] = -1;
      for (EdgeID e : sub.cut_edges) is_cut[e] = 0;
      std::cout << "wrote subproblem (" << a << ", " << b << "): " << sub.vertices.size() << " vertices, flow "
                << sub.flow << " to " << path << "\n";
    }
  }
  if (!fragment_graph_path.empty()) {
    write_metis_graph(fragment_graph, fragment_graph_path);
    std::cout << "wrote fragment graph to " << fragment_graph_path << "\n";
  }
  return 0;
}
