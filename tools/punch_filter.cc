// tools/punch_filter.cc
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "mt-kahypar/partition/preprocessing/filtering/dimacs_io.h"
#include "mt-kahypar/partition/preprocessing/filtering/filtering_pipeline.h"
#include "mt-kahypar/partition/preprocessing/filtering/graph_contraction.h"
#include "mt-kahypar/partition/preprocessing/filtering/metis_io.h"

using namespace mt_kahypar::filtering;

namespace {
void print_usage(const char* prog) {
  std::cerr
      << "Usage: " << prog << " --graph <dimacs.gr> --U <size> "
      << "[--alpha 1.0] [--f 10] [--coverage 2] [--cut-side source|sink|both] [--arc-weights] [--tau 5] "
      << "[--dump-kept-edges <path>] [--dump-partition <path>] [--dump-fragment-graph <path>] [--verbose]\n";
}
}  // namespace

int main(int argc, char** argv) {
  std::string graph_path;
  std::string dump_path;
  std::string partition_path;
  std::string fragment_graph_path;
  NodeWeight U = 0;
  NodeWeight tau = 5;
  double alpha = 1.0;
  double f = 10.0;
  int coverage = 2;
  CutSide cut_side = CutSide::Source;
  bool has_U = false;
  bool use_arc_weights = false;
  bool verbose = false;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        std::exit(1);
      }
      return argv[++i];
    };
    if (arg == "--graph")
      graph_path = next();
    else if (arg == "--U") {
      U = std::stoull(next());
      has_U = true;
    } else if (arg == "--alpha")
      alpha = std::stod(next());
    else if (arg == "--f")
      f = std::stod(next());
    else if (arg == "--coverage")
      coverage = std::stoi(next());
    else if (arg == "--cut-side") {
      const std::string side = next();
      if (side == "source") cut_side = CutSide::Source;
      else if (side == "sink") cut_side = CutSide::Sink;
      else if (side == "both") cut_side = CutSide::Both;
      else {
        print_usage(argv[0]);
        return 1;
      }
    } else if (arg == "--tau")
      tau = std::stoull(next());
    else if (arg == "--dump-kept-edges")
      dump_path = next();
    else if (arg == "--dump-partition")
      partition_path = next();
    else if (arg == "--dump-fragment-graph")
      fragment_graph_path = next();
    else if (arg == "--arc-weights")
      use_arc_weights = true;
    else if (arg == "--verbose")
      verbose = true;
    else {
      print_usage(argv[0]);
      return 1;
    }
  }

  if (graph_path.empty() || !has_U) {
    print_usage(argv[0]);
    return 1;
  }

  FilterGraph graph = read_dimacs_graph(graph_path, use_arc_weights);
  std::cout << "Loaded graph: |V| = " << graph.numNodes()
            << ", |E| = " << graph.numEdges() << "\n";

  FilteringParams params{U, tau, alpha, f, coverage};
  params.cut_side = cut_side;
  params.verbose = verbose;

  const auto start = std::chrono::steady_clock::now();
  FilteringResult result = run_filtering_pipeline(graph, params);
  const auto end = std::chrono::steady_clock::now();
  const double seconds = std::chrono::duration<double>(end - start).count();

  std::cout << "Filtering complete in " << seconds << "s\n";
  std::cout << "|V'| (fragments) = " << result.fragment_size.size() << "\n";

  NodeWeight min_size =
      result.fragment_size.empty() ? 0 : result.fragment_size[0];
  NodeWeight max_size = 0;
  size_t near_cap = 0;
  for (NodeWeight size : result.fragment_size) {
    min_size = std::min(min_size, size);
    max_size = std::max(max_size, size);
    if (size > static_cast<NodeWeight>(0.9 * static_cast<double>(U)))
      ++near_cap;
  }
  std::cout << "fragment size: min=" << min_size << " max=" << max_size
            << " count>0.9*U=" << near_cap << "\n";

  if (!dump_path.empty()) {
    std::ofstream out(dump_path);
    for (const auto& [u, v] : result.kept_edges) out << u << " " << v << "\n";
    std::cout << "Wrote " << result.kept_edges.size() << " kept edges to "
              << dump_path << "\n";
  }

  if (!partition_path.empty()) {
    std::ofstream out(partition_path);
    for (NodeID fragment : result.fragment_id) out << fragment << "\n";
    std::cout << "Wrote partition (" << result.fragment_id.size()
              << " vertices, " << result.fragment_size.size()
              << " fragments) to " << partition_path << "\n";
  }

  if (!fragment_graph_path.empty()) {
    // Fragment graph: one vertex per fragment (weight = fragment size), edge
    // weights = number (or total weight) of original edges between fragments.
    // Any partition of it is a partition of the input graph with the same cut.
    std::vector<NodeID> mapping(result.fragment_id.begin(), result.fragment_id.end());
    FilterGraph fragment_graph =
        contract_graph_by_mapping(graph, mapping, result.fragment_size.size());
    write_metis_graph(fragment_graph, fragment_graph_path);
    std::cout << "Wrote fragment graph (" << fragment_graph.numNodes() << " vertices, "
              << fragment_graph.numEdges() << " edges) to " << fragment_graph_path << "\n";
  }

  return 0;
}
