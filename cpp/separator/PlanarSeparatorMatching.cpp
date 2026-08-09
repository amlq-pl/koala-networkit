#include "separator/PlanarSeparatorMatching.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/structures/UnionFind.hpp"

namespace Koala {
PlanarSeparatorMatching::PlanarSeparatorMatching(NetworKit::Graph &G)
    : graph(G) {}

void PlanarSeparatorMatching::run() { hasRun = true; }

std::vector<NetworKit::Edge>
PlanarSeparatorMatching::reduce_procedure(NetworKit::Graph &subgraph) {
  NetworKit::count N = subgraph.numberOfNodes();
}
} // namespace Koala
