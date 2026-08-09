#include "networkit/Globals.hpp"
#include "networkit/base/Algorithm.hpp"
#include "networkit/graph/Graph.hpp"
#include <map>
#include <vector>
namespace Koala {
class EpsilonPlanarSeparator : public NetworKit::Algorithm {
public:
  std::vector<NetworKit::node> separator;
  std::vector<NetworKit::Graph> connectedComponents;

  void run() override;

  explicit EpsilonPlanarSeparator(
      const NetworKit::Graph &G, double epsilon,
      std::optional<std::map<NetworKit::node, double>> costs = std::nullopt);

private:
  const NetworKit::Graph &graph;
  double epsilon;
  std::vector<double> vertexCost;

  void processComponents(const NetworKit::Graph &graph,
                         std::vector<double> &vertexCost, double epsilon,
                         std::queue<NetworKit::Graph> &Q);
};
} // namespace Koala
