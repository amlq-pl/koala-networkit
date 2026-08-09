#include "separator/PlanarSeparator.hpp"
#include "shortest_path/planar/SuitableRDivision.hpp"

#include <algorithm>
#include <cmath>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <queue>
#include <unordered_set>

using cycle_t = std::vector<NetworKit::node>;

struct CostComputation {
  double insideCost;
  double outsideCost;
  bool insideIsClockwiseArc;
};

namespace {

// BFS from `startNode`, returning (level, parent) vectors indexed by node id.
std::pair<std::vector<NetworKit::node>, std::vector<NetworKit::node>>
performBFSAndFindSpanningTree(const NetworKit::Graph &G,
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

// Number of vertices at each BFS level, indexed by level.
std::vector<NetworKit::node>
findNumberOfVerticesAtLevel(const std::vector<NetworKit::node> &lvl) {
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
  return verticesAtLevel;
}

bool isTreeEdge(NetworKit::node u, NetworKit::node v,
                std::vector<NetworKit::node> &parent) {
  return u == parent[v] || v == parent[u];
}

std::pair<NetworKit::node, NetworKit::node>
findNonTreeEdge(NetworKit::Graph &H, std::vector<NetworKit::node> &parent) {
  NetworKit::node v1 = NetworKit::none;
  NetworKit::node w1 = NetworKit::none;

  // find non tree edge (v1, w1)
  H.forEdges([&](NetworKit::node v, NetworKit::node w) {
    if (v1 != NetworKit::none)
      return; // already found
    if (parent[v] != w && parent[w] != v) {
      v1 = v;
      w1 = w;
    }
  });

  return {v1, w1};
}

std::optional<cycle_t>
buildFundamentalCycle(NetworKit::node v1, NetworKit::node w1,
                      std::vector<NetworKit::node> &parent) {
  // Degenerate H (no non-tree edge): fall back to a level-based separator.
  if (v1 == NetworKit::none) {
    return std::nullopt;
  }

  // Collect ancestors of v1 and w1
  std::vector<NetworKit::node> pathV, pathW;

  // Walk v1 to root
  for (auto v = v1; v != NetworKit::none; v = parent[v])
    pathV.push_back(v);

  // Walk w1 to root
  for (auto w = w1; w != NetworKit::none; w = parent[w])
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
  cycle_t cycle;
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

  return cycle;
}

bool isOnInsideArc(int pos, int posNext, int posPrev) {
  if (posNext < posPrev)
    return pos > posNext && pos < posPrev;
  else
    return pos > posNext || pos < posPrev;
}

int rotationPosition(planar_embedding_t &embedding, NetworKit::node v,
                     NetworKit::node nodeToFindPosition) {
  auto &neighborsV = embedding[v];
  return std::find(neighborsV.begin(), neighborsV.end(), nodeToFindPosition) -
         neighborsV.begin();
}

int wrapIndex(int i, int size) { return ((i % size) + size) % size; }

NetworKit::node findNodeInCycleAtIthPosition(cycle_t &cycle, int i) {
  return cycle[wrapIndex(i, static_cast<int>(cycle.size()))];
}

CostComputation computeSidesCost(cycle_t &cycle, planar_embedding_t &embedding,
                                 std::vector<NetworKit::node> &parent,
                                 std::vector<double> &costs,
                                 NetworKit::node root) {
  int arcTrueCost = 0;
  int arcFalseCost = 0;

  int cycleSize = cycle.size();

  for (int i = 0; i < cycleSize; i++) {
    auto prev = findNodeInCycleAtIthPosition(cycle, i - 1);
    auto v = findNodeInCycleAtIthPosition(cycle, i);
    auto u = findNodeInCycleAtIthPosition(cycle, i + 1);

    auto &neighborsV = embedding[v];
    int posU = rotationPosition(embedding, v, u);
    int posPrev = rotationPosition(embedding, v, prev);

    for (int j = 0; j < static_cast<int>(neighborsV.size()); j++) {
      auto w = neighborsV[j];
      if (w == prev || w == u)
        continue; // skip cycle neighbors

      int cost = 0;
      if (parent[w] == v) {
        cost = costs[w];
      } else if (parent[v] == w) {
        cost = costs[root] - costs[v];
      } else {
        continue; // non-tree edge, no cost
      }

      if (isOnInsideArc(j, posU, posPrev)) {
        arcTrueCost += cost;
      } else {
        arcFalseCost += cost;
      }
    }
  }

  double insideCost = fmax(arcTrueCost, arcFalseCost);
  double outsideCost = fmin(arcTrueCost, arcFalseCost);
  bool isTrueArcInside = arcTrueCost >= arcFalseCost;

  return {insideCost, outsideCost, isTrueArcInside};
}

NetworKit::node findTriangleApex(planar_embedding_t &embedding,
                                 NetworKit::node v, NetworKit::node w,
                                 NetworKit::node otherCycleNeighbor,
                                 bool isTrueArcInside) {
  auto &neighborsV = embedding[v];
  int deg = static_cast<int>(neighborsV.size());

  int posW = rotationPosition(embedding, v, w);
  int posOther = rotationPosition(embedding, v, otherCycleNeighbor);

  int posAfter = wrapIndex(posW + 1, deg);
  int posBefore = wrapIndex(posW - 1, deg);

  bool afterIsInside =
      (isOnInsideArc(posAfter, posOther, posW) == isTrueArcInside);

  return neighborsV[afterIsInside ? posAfter : posBefore];
}
} // namespace

namespace Koala {

PlanarSeparator::PlanarSeparator(const NetworKit::Graph &graph,
                                 std::vector<double> &costs)
    : graph(graph), vertexCost(costs) {}

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
    std::vector<double> costsH(H.upperNodeIdBound(), 0.0);
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
    auto [v1, w1] = findNonTreeEdge(H, parentH);

    // Degenerate H (no non-tree edge): fall back to a level-based separator.
    if (v1 == NetworKit::none) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    auto cycleOpt = buildFundamentalCycle(v1, w1, parentH);

    // Degenerate cycle: fall back to a level-based separator.
    if (!cycleOpt) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    auto cycle = cycleOpt.value();

    if (cycle.size() < 3) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    auto computeParams =
        computeSidesCost(cycle, embeddingH, parentH, costsH, x);

    // if the cost of inside > 2/3 we move to step 9

    cycle_t ci = cycle;
    NetworKit::node vi = v1;
    NetworKit::node wi = w1;
    bool isTrueArcInsideI = computeParams.insideIsClockwiseArc;
    double cost = computeParams.insideCost;
    if (computeParams.insideCost > 2.0 / 3) {
      while (cost > 2.0 / 3) {
        NetworKit::node otherNeighbor = findNodeInCycleAtIthPosition(ci, 1);
        NetworKit::node y = findTriangleApex(embeddingH, vi, wi, otherNeighbor,
                                             isTrueArcInsideI);
        if (isTreeEdge(vi, y, parentH) || isTreeEdge(wi, y, parentH)) {
          NetworKit::node v_next;
          NetworKit::node w_next;

          if (!isTreeEdge(vi, y, parentH)) {
            v_next = vi;
            w_next = y;
          } else {
            v_next = y;
            w_next = wi;
          }

          auto cycleNext =
              buildFundamentalCycle(v_next, w_next, parentH).value();

          auto costParamsNext =
              computeSidesCost(cycleNext, embeddingH, parentH, costsH, x);

          cost = costParamsNext.insideCost;
          isTrueArcInsideI = costParamsNext.insideIsClockwiseArc;
          ci = cycleNext;
          vi = v_next;
          wi = w_next;
        } else {
          auto R1 = buildFundamentalCycle(vi, y, parentH).value();
          auto R2 = buildFundamentalCycle(y, wi, parentH).value();
          auto costParamsR1 =
              computeSidesCost(R1, embeddingH, parentH, costsH, x);
          auto costParamsR2 =
              computeSidesCost(R2, embeddingH, parentH, costsH, x);

          cost = fmax(costParamsR1.insideCost, costParamsR2.insideCost);
          isTrueArcInsideI = costParamsR1.insideCost >= costParamsR2.insideCost
                                 ? costParamsR1.insideIsClockwiseArc
                                 : costParamsR2.insideIsClockwiseArc;
          ci = costParamsR1.insideCost >= costParamsR2.insideCost ? R1 : R2;
          vi = costParamsR1.insideCost >= costParamsR2.insideCost ? vi : y;
          wi = costParamsR1.insideCost >= costParamsR2.insideCost ? y : wi;
        }
      }
      cycle = ci;
    }

    // Step 9: iteratively shrink the cycle until insideCost <= 2/3 * total.

    // ---- Step 10: extract separator and partitions back to G ----
    extractSeparatorAndPartitions(G, lvl, l0, l2, cycle, x);
  }
  hasRun = true;
}

std::vector<PlanarSeparator::Side>
PlanarSeparator::markInsideOutside(const NetworKit::Graph &G,
                                   const std::vector<NetworKit::node> &cycle) {
  std::vector<PlanarSeparator::Side> side(G.upperNodeIdBound());
  std::vector<bool> vis(G.upperNodeIdBound(), false);
  std::queue<NetworKit::node> Q;
  for (auto v : cycle) {
    side[v] = PlanarSeparator::Side::ON_CYCLE;
    vis[v] = true;
  }

  G.forNodes([&](NetworKit::node v) {

  });
}

void PlanarSeparator::extractSeparatorAndPartitions(
    const NetworKit::Graph &G, const std::vector<NetworKit::node> &lvl,
    NetworKit::node l0, NetworKit::node l2,
    const std::vector<NetworKit::node> &cycle, NetworKit::node x) {}

bool PlanarSeparator::areConnectedComponentsEligibleForPartition(
    NetworKit::ConnectedComponents &components) {

  return highestCost < 2.0 / 3;
}

NetworKit::Graph PlanarSeparator::findLargestCostComponent(
    NetworKit::ConnectedComponents &components) {
  NetworKit::Graph mostCostly;
  double highestCost = -1.0;

  for (auto cc : components.getComponents()) {
    double ccCost = 0.0;
    for (auto v : cc) {
      ccCost += vertexCost[v];
    }

    if (ccCost > highestCost) {
      mostCostly = cc;
      highestCost = ccCost;
    }
  }
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
        partition.A.insert(partition.A.end(), component.begin(),
                           component.end());
        totalSize += component.size();
      } else {
        partition.B.insert(partition.B.end(), component.begin(),
                           component.end());
      }
    }
  } else {
    // largestComponent has cost > 2/3, put its nodes in A, rest in B
    std::unordered_set<NetworKit::node> largestNodes;
    largestComponent.forNodes(
        [&](NetworKit::node v) { largestNodes.insert(v); });

    partition.A.assign(largestNodes.begin(), largestNodes.end());

    graph.forNodes([&](NetworKit::node v) {
      if (largestNodes.find(v) == largestNodes.end()) {
        partition.B.push_back(v);
      }
    });
  }
}

void PlanarSeparator::cleanPartitions() {
  partition.separator.clear();
  partition.A.clear();
  partition.B.clear();
}

const PlanarSeparator::Partition &PlanarSeparator::getPartition() const {
  assureFinished();
  return partition;
}

const std::vector<NetworKit::node> &PlanarSeparator::getSeparator() const {
  assureFinished();
  return partition.separator;
}

const std::vector<NetworKit::node> &PlanarSeparator::getPartitionA() const {
  assureFinished();
  return partition.A;
}

const std::vector<NetworKit::node> &PlanarSeparator::getPartitionB() const {
  assureFinished();
  return partition.B;
}

} // namespace Koala
