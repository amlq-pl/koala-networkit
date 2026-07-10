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

    // Side of a vertex with respect to the current cycle in step 9.
    enum class Side : uint8_t { OUTSIDE = 0, ON_CYCLE = 1, INSIDE = 2 };

    // Classify every vertex of H as OUTSIDE / ON_CYCLE / INSIDE with respect to
    // the given cycle, using a single BFS seeded from a vertex on the inside
    // arc of the embedding. Runs in O(|V(H)| + |E(H)|).
    std::vector<Side> markInsideOutside(
        const NetworKit::Graph& H,
        const PlanarGraphTools::planar_embedding_t& embeddingH,
        const std::vector<NetworKit::node>& cycle,
        bool insideIsClockwiseArc);

    // Step 9, Case B handler. Called when neither (v1, y) nor (y, w1) is a
    // tree edge. Picks whichever of the two triangle edges has more cost
    // inside its sub-cycle, updates `cycle`, `side`, `insideCost`, `v1`,
    // `w1` in place. This is the slow (O(n) per call) but simple version;
    // it can be replaced by the alternating-scan O(1)-amortized variant
    // without changing the rest of the algorithm.
    void applyStep9CaseB(
        const NetworKit::Graph& H,
        const std::vector<NetworKit::node>& parentH,
        NetworKit::node y,
        NetworKit::node& v1,
        NetworKit::node& w1,
        std::vector<NetworKit::node>& cycle,
        std::vector<Side>& side,
        int& insideCost);

    // Step 10: build the final separator / partitionA / partitionB sets in G
    // from the cycle found in H. Handles synthetic vertex x, levels l0/l2,
    // levels < l0 (contracted), levels > l2 (removed), and other components.
    void extractSeparatorAndPartitions(
        const NetworKit::Graph& G,
        const std::vector<NetworKit::node>& lvl,
        NetworKit::node l0,
        NetworKit::node l2,
        const std::vector<NetworKit::node>& cycle,
        NetworKit::node x,
        const std::vector<Side>& side);

    // Fallback used when the cycle-finding part of the algorithm degenerates
    // (graph too small / H has no non-tree edge): use BFS level l1 as the
    // separator. Levels < l1 go to A, levels > l1 go to B. Vertices outside
    // the largest CC go to B. This always satisfies |A|, |B| <= n/2 <= 2n/3.
    void fallbackLevelSeparator(
        const std::vector<NetworKit::node>& lvl,
        NetworKit::node l1);
};

}  // namespace Koala
