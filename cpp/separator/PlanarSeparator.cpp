#include "separator/PlanarSeparator.hpp"

#include <algorithm>
#include <cmath>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <queue>
#include <unordered_set>

namespace Koala {

PlanarSeparator::PlanarSeparator(const NetworKit::Graph &graph)
    : graph(graph) {}

void PlanarSeparator::run() {
  cleanPartitions();
  // Step 1: Find a planar embedding of the graph G
  auto embedding = PlanarGraphTools::findPlanarEmbedding(graph);

  // Step 2: Find connected components of the graph G
  auto components = NetworKit::ConnectedComponents(graph);
  components.run();

  if (areConnectedComponentsEligibleForPartition(components)) {
    findSeparatorFromComponents(components);
  } else {
    // Step 3: Perform BFS
    auto G = components.extractLargestConnectedComponent(graph, false);
    NetworKit::count n = G.numberOfNodes();

    // extractLargestConnectedComponent(..., false) preserves the original
    // node ids (with gaps), so node 0 may not belong to G. Start the BFS
    // from a node that actually exists in the component.
    NetworKit::node root = NetworKit::none;
    G.forNodes([&](NetworKit::node v) {
      if (root == NetworKit::none)
        root = v;
    });

    auto [lvl, parent] = performBFSAndFindSpanningTree(G, root);

    auto verticesAtLevel = findNumberOfVerticesAtLevel(lvl);

    // prefix sum for easier counting of vertices at levels
    std::vector<NetworKit::count> prefixSum(verticesAtLevel.size());
    std::partial_sum(verticesAtLevel.begin(), verticesAtLevel.end(),
                     prefixSum.begin());

    // Step 4: Find level l1, and k
    // l1 is the smallest level such that the number of vertices in levels 0..l1
    // is at least n/2
    NetworKit::node l1 =
        std::lower_bound(prefixSum.begin(), prefixSum.end(), n / 2) -
        prefixSum.begin();
    // k is the number of vertices in levels 0..l1
    NetworKit::count k = prefixSum[l1];

    // Step 5: Find levels l0 and l2
    // l0 is the largest level such that the number of vertices in level l0 +
    // 2*(l1-l0) is at most 2*sqrt(k)
    NetworKit::node l0 = 0;
    for (int i = 0; i <= l1; i++) {
      if (verticesAtLevel[i] + 2 * (l1 - i) <= 2 * sqrt(k)) {
        l0 = i;
      } else
        break;
    }
    // l2 is the smallest level such that the number of vertices in level l2 +
    // 2*(l2-l1-1) is at most 2*sqrt(n-k)
    NetworKit::node l2 = l1 + 1;
    for (int i = l1 + 1; i < verticesAtLevel.size(); i++) {
      if (verticesAtLevel[i] + 2 * (i - l1 - 1) <= 2 * sqrt(n - k)) {
        l2 = i;
        break;
      }
    }

    // Step 6: Remove vertices at level >= l2, contract levels <= l0 into vertex
    // x
    std::unordered_set<NetworKit::node> verticesToKeep;
    G.forNodes([&](NetworKit::node v) {
      if (lvl[v] < l2)
        verticesToKeep.insert(v);
    });

    auto H = NetworKit::GraphTools::subgraphFromNodes(G, verticesToKeep);
    NetworKit::node x = H.addNode();

    std::vector<bool> isConnectedToX(H.upperNodeIdBound(), false);

    // Walk around the subtree (levels 0..l0) using the embedding
    // and rewire boundary edges to x
    std::unordered_set<NetworKit::node> subtreeNodes;
    G.forNodes([&](NetworKit::node v) {
      if (lvl[v] <= l0)
        subtreeNodes.insert(v);
    });

    for (auto v : subtreeNodes) {
      for (auto w : embedding[v]) {
        if (subtreeNodes.count(w) > 0)
          continue; // both in subtree, skip
        if (!H.hasNode(w))
          continue; // deleted (level >= l2), skip
        if (!isConnectedToX[w]) {
          isConnectedToX[w] = true;
          H.addEdge(x, w);
        }
      }
    }

    // Remove the contracted subtree nodes from H
    for (auto v : subtreeNodes) {
      H.removeNode(v);
    }

    // Step 7 compute spanngin tree of H and costs of subtrees hanging off it
    auto [lvlH, parentH] = performBFSAndFindSpanningTree(H, x);
    std::vector<std::vector<NetworKit::node>> childrenH(
        H.upperNodeIdBound(), std::vector<NetworKit::node>());
    H.forNodes([&](NetworKit::node v) {
      if (parentH[v] != NetworKit::none) {
        childrenH[parentH[v]].push_back(v);
      }
    });
    std::vector<NetworKit::node> costsH(H.upperNodeIdBound(), 0);
    // Iterative post-order accumulation of subtree sizes, rooted at x.
    // Recursion would overflow the stack on large graphs (tens of
    // thousands of vertices).
    {
      std::vector<std::pair<NetworKit::node, size_t>> stk;
      costsH[x] = 1;
      stk.emplace_back(x, 0);
      while (!stk.empty()) {
        NetworKit::node v = stk.back().first;
        size_t &idx = stk.back().second;
        if (idx < childrenH[v].size()) {
          NetworKit::node c = childrenH[v][idx++];
          costsH[c] = 1;
          stk.emplace_back(c, 0);
        } else {
          NetworKit::node finished = v;
          stk.pop_back();
          if (!stk.empty())
            costsH[stk.back().first] += costsH[finished];
        }
      }
    }

    H = PlanarGraphTools::makeMaximalPlanar(H);
    auto embeddingH = PlanarGraphTools::findPlanarEmbedding(H);

    // Step 8
    NetworKit::node v1 = NetworKit::none;
    NetworKit::node w1 = NetworKit::none;

    // find non tree edge (v1, w1)
    H.forEdges([&](NetworKit::node v, NetworKit::node w) {
      if (v1 != NetworKit::none)
        return; // already found
      if (parentH[v] != w && parentH[w] != v) {
        v1 = v;
        w1 = w;
      }
    });

    // Degenerate H (no non-tree edge): fall back to a level-based separator.
    if (v1 == NetworKit::none) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    // Collect ancestors of v1 and w1
    std::vector<NetworKit::node> pathV, pathW;

    // Walk v1 to root
    for (auto v = v1; v != NetworKit::none; v = parentH[v])
      pathV.push_back(v);

    // Walk w1 to root
    for (auto w = w1; w != NetworKit::none; w = parentH[w])
      pathW.push_back(w);

    // Find LCA — first common vertex
    std::unordered_set<NetworKit::node> ancestorsV(pathV.begin(), pathV.end());
    NetworKit::node lca = NetworKit::none;
    for (auto w : pathW) {
      if (ancestorsV.count(w)) {
        lca = w;
        break;
      }
    }

    // build cycle
    std::vector<NetworKit::node> cycle;
    for (auto v : pathV) {
      cycle.push_back(v);
      if (v == lca)
        break;
    }

    std::vector<NetworKit::node> pathWToLCA;
    for (auto w : pathW) {
      if (w == lca)
        break;
      pathWToLCA.push_back(w);
    }
    std::reverse(pathWToLCA.begin(), pathWToLCA.end());
    cycle.insert(cycle.end(), pathWToLCA.begin(), pathWToLCA.end());

    // Degenerate cycle: fall back to a level-based separator.
    if (cycle.size() < 3) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    // walk the cycle and find the costs
    int inside = 0;
    int outside = 0;

    // check which direction the cycle goes and ensure it's clockwise
    auto &neighborsOfCycle1 = embeddingH[cycle[1]];
    int posA = std::find(neighborsOfCycle1.begin(), neighborsOfCycle1.end(),
                         cycle[0]) -
               neighborsOfCycle1.begin();
    int posC = std::find(neighborsOfCycle1.begin(), neighborsOfCycle1.end(),
                         cycle[2]) -
               neighborsOfCycle1.begin();

    // Walk clockwise from posA; if we hit posC, the cycle is clockwise
    bool clockwise = false;
    for (int i = 1; i < (int)neighborsOfCycle1.size(); i++) {
      int pos = (posA + i) % neighborsOfCycle1.size();
      if (pos == posC) {
        clockwise = true;
        break;
      }
    }

    if (!clockwise) {
      std::reverse(cycle.begin(), cycle.end());
    }

    int cycleSize = cycle.size();

    for (int i = 0; i < cycleSize; i++) {
      auto prev = cycle[(cycleSize + i - 1) % cycleSize];
      auto v = cycle[i];
      auto u = cycle[(i + 1) % cycleSize];

      auto neighborsV = embeddingH[v];
      int posU = std::find(neighborsV.begin(), neighborsV.end(), u) -
                 neighborsV.begin();
      int posPrev = std::find(neighborsV.begin(), neighborsV.end(), prev) -
                    neighborsV.begin();

      auto isInside = [&](int pos) {
        // inside = clockwise arc from posU to posPrev (exclusive)
        if (posU < posPrev) {
          return pos > posU && pos < posPrev;
        } else {
          return pos > posU || pos < posPrev;
        }
      };

      for (int j = 0; j < (int)neighborsV.size(); j++) {
        auto w = neighborsV[j];
        if (w == prev || w == u)
          continue; // skip cycle neighbors

        int cost = 0;
        if (parentH[w] == v) {
          cost = costsH[w];
        } else if (parentH[v] == w) {
          cost = costsH[x] - costsH[v];
        } else {
          continue; // non-tree edge, no cost
        }

        if (isInside(j)) {
          inside += cost;
        } else {
          outside += cost;
        }
      }
    }

    bool insideIsClockwiseArc = true;

    if (outside > inside) {
      std::swap(inside, outside);
      insideIsClockwiseArc = false;
    }

    // if the cost of inside > 2/3 we move to step 9
    //

    // Step 9: iteratively shrink the cycle until insideCost <= 2/3 * total.

    int i = 1;

    // ---- Step 10: extract separator and partitions back to G ----
    extractSeparatorAndPartitions(G, lvl, l0, l2, cycle, x, side);
  }
  hasRun = true;
}

bool PlanarSeparator::areConnectedComponentsEligibleForPartition(
    NetworKit::ConnectedComponents &components) {
  auto largestComponent =
      components.extractLargestConnectedComponent(graph, false);
  NetworKit::count n = this->graph.numberOfNodes();

  return largestComponent.numberOfNodes() < 2 * n / 3;
}

void PlanarSeparator::findSeparatorFromComponents(
    NetworKit::ConnectedComponents &components) {
  auto largestComponent =
      components.extractLargestConnectedComponent(graph, false);
  NetworKit::count n = this->graph.numberOfNodes();

  // in both cases the separator will be empty, so we do not assign it
  // explicitly

  if (largestComponent.numberOfNodes() < n / 3) {
    int totalSize = 0;
    for (auto component : components.getComponents()) {
      if (totalSize < n / 3) {
        partitionA.insert(partitionA.end(), component.begin(), component.end());
        totalSize += component.size();
      } else {
        partitionB.insert(partitionB.end(), component.begin(), component.end());
      }
    }
  } else {
    // largestComponent has cost > 2/3, put its nodes in A, rest in B
    std::unordered_set<NetworKit::node> largestNodes;
    largestComponent.forNodes(
        [&](NetworKit::node v) { largestNodes.insert(v); });

    partitionA.assign(largestNodes.begin(), largestNodes.end());

    graph.forNodes([&](NetworKit::node v) {
      if (largestNodes.find(v) == largestNodes.end()) {
        partitionB.push_back(v);
      }
    });
  }
}

void PlanarSeparator::cleanPartitions() {
  separator.clear();
  partitionA.clear();
  partitionB.clear();
}

std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
PlanarSeparator::performBFSAndFindSpanningTree(const NetworKit::Graph &G,
                                               NetworKit::node startNode) {
  NetworKit::count n = G.upperNodeIdBound();
  std::vector<bool> vis(n, false);
  std::queue<NetworKit::node> Q;

  std::vector<NetworKit::node> lvl(n, NetworKit::none);
  std::vector<NetworKit::node> parent(n, NetworKit::none);

  parent[startNode] = NetworKit::none;
  lvl[startNode] = 0;
  vis[startNode] = true;
  Q.push(startNode);

  while (!Q.empty()) {
    auto u = Q.front();
    Q.pop();

    G.forNeighborsOf(u, [&](NetworKit::node v) {
      if (!vis[v]) {
        vis[v] = true;
        parent[v] = u;
        lvl[v] = lvl[u] + 1;
        Q.push(v);
      }
    });
  }

  return {lvl, parent};
}

std::vector<NetworKit::node> PlanarSeparator::findNumberOfVerticesAtLevel(
    const std::vector<NetworKit::node> &lvl) {
  NetworKit::count maxLevel = 0;
  for (auto l : lvl) {
    if (l != NetworKit::none)
      maxLevel = std::max(maxLevel, (NetworKit::count)l);
  }
  std::vector<NetworKit::count> verticesAtLevel(maxLevel + 1, 0);
  for (auto l : lvl) {
    if (l != NetworKit::none)
      verticesAtLevel[l]++;
  }
  return verticesAtLevel; // BUT: return type is wrong — see below
}

const std::vector<NetworKit::node> &PlanarSeparator::getSeparator() const {
  assureFinished();
  return separator;
}

const std::vector<NetworKit::node> &PlanarSeparator::getPartitionA() const {
  assureFinished();
  return partitionA;
}

const std::vector<NetworKit::node> &PlanarSeparator::getPartitionB() const {
  assureFinished();
  return partitionB;
}

} // namespace Koala
