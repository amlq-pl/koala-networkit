#include "separator/MISP.hpp"
#include "separator/EpsilonPlanarSeparator.hpp"

namespace Koala {
template <typename T>
MISP<T>::MISP(const NetworKit::Graph &G, double epsilon,
              ComponentSolver componentSolver,
              std::optional<std::map<NetworKit::node, double>> costs)
    : graph(G), epsilon(epsilon), componentSolver(componentSolver),
      costs(costs) {}

template <typename T> void MISP<T>::run() {
  EpsilonPlanarSeparator sep(graph, epsilon, costs);
  sep.run();

  for (const auto &cc : sep.connectedComponents) {
    std::vector<T> componentSolution = componentSolver(cc);
    maximum_independent_set.insert(maximum_independent_set.end(),
                                   componentSolution.begin(),
                                   componentSolution.end());
  }
  hasRun = true;
}
} // namespace Koala
