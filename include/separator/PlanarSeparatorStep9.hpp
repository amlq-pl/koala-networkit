#pragma once

#include <networkit/graph/Graph.hpp>
#include <unordered_map>
#include <vector>

namespace Koala {

// Step 9 of the Lipton-Tarjan planar separator, exposed on its own so it can be
// benchmarked with a *precomputed* embedding -- independent of the (slow)
// planarity routine. Given the triangulated graph H with its rotation embedding
// (embeddingH / idxOf), a spanning tree (parentH) and subtree costs (costsH),
// it shrinks the fundamental cycle of the non-tree edge (v1, w1) into a
// balanced separating cycle and returns its vertices (empty if degenerate).
// `x` is the contracted cost-0 super-vertex, or NetworKit::none if none.
std::vector<NetworKit::node> shrinkFundamentalCycle(
    const NetworKit::Graph &H,
    std::unordered_map<NetworKit::node, std::vector<NetworKit::node>>
        &embeddingH,
    std::vector<std::unordered_map<NetworKit::node, int>> &idxOf,
    std::vector<NetworKit::node> &parentH, std::vector<double> &costsH,
    std::vector<double> &vertexCost, NetworKit::node x, NetworKit::node v1,
    NetworKit::node w1);

} // namespace Koala
