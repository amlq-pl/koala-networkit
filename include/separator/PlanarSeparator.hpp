#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <graph/PlanarGraphTools.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <vector>

namespace Koala {

class PlanarSeparator : public NetworKit::Algorithm {
 public:
    explicit PlanarSeparator(const NetworKit::Graph& graph);

    void run() override;

    const std::vector<NetworKit::node>& getSeparator() const;
    const std::vector<NetworKit::node>& getPartitionA() const;
    const std::vector<NetworKit::node>& getPartitionB() const;

 private:
    const NetworKit::Graph& graph;

    std::vector<NetworKit::node> separator;
    std::vector<NetworKit::node> partitionA;
    std::vector<NetworKit::node> partitionB;

    void cleanPartitions();

    bool areConnectedComponentsEligibleForPartition(
        NetworKit::ConnectedComponents& components);

    void findSeparatorFromComponents(
        NetworKit::ConnectedComponents& components);

    std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
        performBFSAndFindSpanningTree(const NetworKit::Graph& G, NetworKit::node startNode = 0);

    std::vector<NetworKit::node> findNumberOfVerticesAtLevel(
        const std::vector<NetworKit::node>& lvl);
};

}  // namespace Koala
