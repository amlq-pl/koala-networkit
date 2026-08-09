#pragma once

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#include <graph/PlanarGraphTools.hpp>

namespace Koala {

class PlanarSeparator : public NetworKit::Algorithm {
public:
  struct Partition {
    std::vector<NetworKit::node> separator;
    std::vector<NetworKit::node> A;
    std::vector<NetworKit::node> B;
  };

  explicit PlanarSeparator(const NetworKit::Graph &graph,
                           std::vector<double> &vertexCost);

  void run() override;

  const Partition &getPartition() const;

  const std::vector<NetworKit::node> &getSeparator() const;
  const std::vector<NetworKit::node> &getPartitionA() const;
  const std::vector<NetworKit::node> &getPartitionB() const;

private:
  const NetworKit::Graph &graph;

  std::vector<double> &vertexCost;

  Partition partition;
  enum class Side : uint8_t { OUTSIDE = 0, ON_CYCLE = 1, INSIDE = 2 };

  void cleanPartitions();

  bool areConnectedComponentsEligibleForPartition(
      NetworKit::ConnectedComponents &components);

  NetworKit::Graph
  findLargestCostComponent(NetworKit::ConnectedComponents &components);

  std::vector<Side>
  markInsideOutside(const NetworKit::Graph &G,
                    const std::vector<NetworKit::node> &cycle);

  void findSeparatorFromComponents(NetworKit::ConnectedComponents &components);

  void extractSeparatorAndPartitions(const NetworKit::Graph &G,
                                     const std::vector<NetworKit::node> &lvl,
                                     NetworKit::node l0, NetworKit::node l2,
                                     const std::vector<NetworKit::node> &cycle,
                                     NetworKit::node x);
  void fallbackLevelSeparator(const std::vector<NetworKit::node> &lvl,
                              NetworKit::node l1);
};

} // namespace Koala
