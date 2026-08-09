#include "separator/EpsilonPlanarSeparator.hpp"
#include "networkit/Globals.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "separator/PlanarSeparator.hpp"
#include <queue>
#include <unordered_set>
#include <vector>

namespace {
std::pair<std::vector<NetworKit::Graph>, std::vector<NetworKit::Graph>>
findCostlyComponents(const NetworKit::Graph &G, std::vector<double> &costs,
                     double epsilon) {
  std::vector<NetworKit::Graph> costlyCCs;
  std::vector<NetworKit::Graph> rest;
  NetworKit::ConnectedComponents components(G);
  components.run();
  auto ccs = components.getComponents();

  for (auto &cc : ccs) {
    double cost = 0.0;
    for (auto v : cc)
      cost += costs[v];

    if (cost > epsilon)
      costlyCCs.push_back(NetworKit::GraphTools::subgraphFromNodes(
          G, std::unordered_set<NetworKit::node>(cc.begin(), cc.end())));
    else
      rest.push_back(NetworKit::GraphTools::subgraphFromNodes(
          G, std::unordered_set<NetworKit::node>(cc.begin(), cc.end())));
  }

  return {costlyCCs, rest};
}

} // namespace

namespace Koala {
EpsilonPlanarSeparator::EpsilonPlanarSeparator(
    const NetworKit::Graph &G, double epsilon,
    std::optional<std::map<NetworKit::node, double>> costs)
    : graph(G), epsilon(epsilon) {
  vertexCost.assign(graph.upperNodeIdBound(), 0.0);
  if (costs.has_value()) {
    for (const auto &[v, c] : *costs) {
      vertexCost[v] = c;
    }
  } else {
    NetworKit::count n = graph.numberOfNodes();
    double uniform = n > 0 ? 1.0 / static_cast<double>(n) : 0.0;
    graph.forNodes([&](NetworKit::node v) { vertexCost[v] = uniform; });
  }
}

void EpsilonPlanarSeparator::run() {
  separator.clear();
  connectedComponents.clear();

  std::queue<NetworKit::Graph> Q;
  processComponents(graph, vertexCost, epsilon, Q);

  while (!Q.empty()) {
    NetworKit::Graph K = std::move(Q.front());
    Q.pop();

    PlanarSeparator sep(K, vertexCost);
    sep.run();
    const auto &partition = sep.getPartition();

    for (auto v : partition.separator)
      separator.push_back(v);

    auto graphA = NetworKit::GraphTools::subgraphFromNodes(
        K, std::unordered_set<NetworKit::node>(partition.A.begin(),
                                               partition.A.end()));
    auto graphB = NetworKit::GraphTools::subgraphFromNodes(
        K, std::unordered_set<NetworKit::node>(partition.B.begin(),
                                               partition.B.end()));

    processComponents(graphA, vertexCost, epsilon, Q);
    processComponents(graphB, vertexCost, epsilon, Q);
  }
  hasRun = true;
}

void EpsilonPlanarSeparator::processComponents(
    const NetworKit::Graph &graph, std::vector<double> &vertexCost,
    double epsilon, std::queue<NetworKit::Graph> &Q) {
  auto [costlyCCs, cheapCCs] = findCostlyComponents(graph, vertexCost, epsilon);

  for (auto &cc : cheapCCs)
    connectedComponents.push_back(cc);

  for (auto &cc : costlyCCs) {
    Q.emplace(cc);
  }
}
} // namespace Koala
