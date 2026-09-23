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
      iss >> u1 >> v1 >> w;  // w is parsed but ignored -- filtering minimizes border nodes, not cut weight
      if (u1 == 0 || v1 == 0 || u1 > n || v1 > n) {
        throw std::runtime_error("DIMACS arc references an out-of-range node id: " + path);
      }
      const NodeID u = static_cast<NodeID>(u1 - 1);
      const NodeID v = static_cast<NodeID>(v1 - 1);
      if (u == v) continue;
      const uint64_t key = undirected_key(u, v);
      edges_by_key.emplace(key, EdgeListEntry{u, v, static_cast<EdgeWeight>(1)});
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
