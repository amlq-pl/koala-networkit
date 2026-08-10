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
  double arcTrueCost = 0.0;
  double arcFalseCost = 0.0;

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

      double cost = 0.0;
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
} // namespace

namespace Koala {

PlanarSeparator::PlanarSeparator(const NetworKit::Graph &graph,
                                 const std::vector<double> &costs)
    : graph(graph), vertexCost(costs) {}

PlanarSeparator::PlanarSeparator(const NetworKit::Graph &graph)
    : graph(graph), vertexCost(graph.upperNodeIdBound(), 1.0) {}

void PlanarSeparator::run() {
  cleanPartitions();

  // Nothing to separate in an empty graph.
  if (graph.numberOfNodes() == 0) {
    hasRun = true;
    return;
  }

  // Normalize the vertex costs so they sum to 1, matching the paper's model
  // ("nonnegative vertex costs summing to no more than one"). Every balance
  // bound is then the literal constant 2/3, never a multiple of an arbitrary
  // total.
  double totalCost = 0.0;
  graph.forNodes([&](NetworKit::node v) { totalCost += vertexCost[v]; });
  if (totalCost > 0.0) {
    graph.forNodes([&](NetworKit::node v) { vertexCost[v] /= totalCost; });
  }

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

    // Cost accumulated per BFS level (the weighted analogue of the counts).
    std::vector<double> costAtLevel(verticesAtLevel.size(), 0.0);
    double costG = 0.0;
    G.forNodes([&](NetworKit::node v) {
      costAtLevel[lvl[v]] += vertexCost[v];
      costG += vertexCost[v];
    });
    std::vector<double> prefixCost(costAtLevel.size());
    std::partial_sum(costAtLevel.begin(), costAtLevel.end(),
                     prefixCost.begin());

    // Step 4: Find level l1, and k
    // l1 is the smallest level such that the total cost of levels 0..l1 is at
    // least half of the total cost of G.
    NetworKit::node l1 =
        std::lower_bound(prefixCost.begin(), prefixCost.end(), 0.5 * costG) -
        prefixCost.begin();
    if (l1 >= verticesAtLevel.size())
      l1 = verticesAtLevel.size() - 1;
    // k is the number of vertices in levels 0..l1 (drives the sqrt bounds).
    NetworKit::count k = prefixSum[l1];

    // Step 5: Find levels l0 and l2
    // l0 is the largest level <= l1 such that |L(l0)| + 2*(l1-l0) <= 2*sqrt(k).
    // The condition is not monotone, so scan downward from l1 and take the
    // first (largest) level that satisfies it.
    NetworKit::node l0 = 0;
    for (int i = static_cast<int>(l1); i >= 0; i--) {
      if (verticesAtLevel[i] + 2.0 * (static_cast<int>(l1) - i) <=
          2.0 * std::sqrt(static_cast<double>(k))) {
        l0 = static_cast<NetworKit::node>(i);
        break;
      }
    }
    // l2 is the smallest level > l1 such that
    // |L(l2)| + 2*(l2-l1-1) <= 2*sqrt(n-k). If none exists (e.g. l1 is the last
    // level) it defaults to l1+1, a virtual empty ring below the middle.
    NetworKit::node l2 = l1 + 1;
    for (int i = static_cast<int>(l1) + 1;
         i < static_cast<int>(verticesAtLevel.size()); i++) {
      if (verticesAtLevel[i] + 2.0 * (i - static_cast<int>(l1) - 1) <=
          2.0 * std::sqrt(static_cast<double>(n - k))) {
        l2 = static_cast<NetworKit::node>(i);
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
      // x is the contracted super-vertex (levels <= l0); it has no cost of its
      // own. Every other node contributes its own vertex cost.
      costsH[x] = 0.0;
      stk.emplace_back(x, 0);
      while (!stk.empty()) {
        NetworKit::node v = stk.back().first;
        size_t &idx = stk.back().second;
        if (idx < childrenH[v].size()) {
          NetworKit::node c = childrenH[v][idx++];
          costsH[c] = vertexCost[c];
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

    // Steps 8-9 (O(n^2) variant). By Lemma 2, among all nontree edges of the
    // triangulated graph, the fundamental cycle that minimizes the larger of
    // its two side costs separates the graph so that neither side exceeds 2/3.
    // The linear-time variant reaches that cycle by iteratively shrinking one
    // candidate; here we instead price every nontree edge's fundamental cycle
    // and keep the best one. The triangulation edges (including those closing
    // the outer face) are exactly the extra nontree edges that make a balanced
    // cycle available even when the spanning tree is shallow (e.g. a wheel
    // rooted at its hub, whose balanced separator is a hub-to-rim "diameter").
    // Cost: O(n) nontree edges times O(n) to price each cycle, i.e. O(n^2).
    double bestInside =
        computeSidesCost(cycle, embeddingH, parentH, costsH, x).insideCost;
    H.forEdges([&](NetworKit::node a, NetworKit::node b) {
      if (isTreeEdge(a, b, parentH))
        return;
      auto candidateOpt = buildFundamentalCycle(a, b, parentH);
      if (!candidateOpt || candidateOpt->size() < 3)
        return;
      double inside =
          computeSidesCost(candidateOpt.value(), embeddingH, parentH, costsH, x)
              .insideCost;
      if (inside < bestInside) {
        bestInside = inside;
        cycle = candidateOpt.value();
      }
    });

    // ---- Step 10: extract separator and partitions back to G ----
    extractSeparatorAndPartitions(G, lvl, l0, l2, cycle, x);
  }
  hasRun = true;
}

void PlanarSeparator::assignComponentsToSides(
    const std::unordered_set<NetworKit::node> &separatorNodes) {
  // Connected components of graph \ separatorNodes, each with its total cost.
  // Keeping every component whole guarantees that no edge ever crosses between
  // A and B.
  std::vector<std::pair<double, std::vector<NetworKit::node>>> components;
  std::vector<bool> visited(graph.upperNodeIdBound(), false);
  graph.forNodes([&](NetworKit::node s) {
    if (visited[s] || separatorNodes.count(s) > 0)
      return;
    std::vector<NetworKit::node> comp;
    double compCost = 0.0;
    std::queue<NetworKit::node> Q;
    Q.push(s);
    visited[s] = true;
    while (!Q.empty()) {
      NetworKit::node u = Q.front();
      Q.pop();
      comp.push_back(u);
      compCost += vertexCost[u];
      graph.forNeighborsOf(u, [&](NetworKit::node w) {
        if (!visited[w] && separatorNodes.count(w) == 0) {
          visited[w] = true;
          Q.push(w);
        }
      });
    }
    components.emplace_back(compCost, std::move(comp));
  });

  partition.separator.assign(separatorNodes.begin(), separatorNodes.end());

  // Longest-processing-time greedy: assign the costliest component to the
  // lighter side. The separator (Lemma 2 / Lemma 3) guarantees every component
  // has cost at most 2/3 and the total is at most 1, so both sides stay within
  // the 2/3 bound.
  std::sort(components.begin(), components.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });

  double costA = 0.0;
  double costB = 0.0;
  for (auto &[compCost, comp] : components) {
    if (costA <= costB) {
      partition.A.insert(partition.A.end(), comp.begin(), comp.end());
      costA += compCost;
    } else {
      partition.B.insert(partition.B.end(), comp.begin(), comp.end());
      costB += compCost;
    }
  }
}

void PlanarSeparator::extractSeparatorAndPartitions(
    const NetworKit::Graph &G, const std::vector<NetworKit::node> &lvl,
    NetworKit::node l0, NetworKit::node l2,
    const std::vector<NetworKit::node> &cycle, NetworKit::node x) {
  // The separator is the two level rings L(l0) and L(l2) (which keep the
  // above / middle / below regions apart) together with the cycle vertices
  // (which split the middle into its inside and outside). The contracted
  // super-vertex x is virtual and is never part of the graph.
  std::unordered_set<NetworKit::node> separatorNodes;
  G.forNodes([&](NetworKit::node v) {
    if (lvl[v] == l0 || lvl[v] == l2)
      separatorNodes.insert(v);
  });
  for (auto c : cycle) {
    if (c != x)
      separatorNodes.insert(c);
  }

  assignComponentsToSides(separatorNodes);
}

void PlanarSeparator::fallbackLevelSeparator(
    const std::vector<NetworKit::node> &lvl, NetworKit::node l1) {
  // Degenerate case (no usable fundamental cycle): fall back to the BFS level
  // l1 as the separator. It splits the component into levels < l1 and > l1,
  // each holding less than half of the cost.
  std::unordered_set<NetworKit::node> separatorNodes;
  for (NetworKit::node v = 0; v < lvl.size(); v++) {
    if (lvl[v] == l1 && graph.hasNode(v))
      separatorNodes.insert(v);
  }

  assignComponentsToSides(separatorNodes);
}

bool PlanarSeparator::areConnectedComponentsEligibleForPartition(
    NetworKit::ConnectedComponents &components) {
  // Costs are normalized to sum to 1, so "eligible" means no single component
  // has cost exceeding 2/3: then the whole components can be split into two
  // sides each within 2/3, with an empty separator. Otherwise the heavy
  // component must be separated further.
  double highestCost = 0.0;
  for (const auto &cc : components.getComponents()) {
    double ccCost = 0.0;
    for (auto v : cc)
      ccCost += vertexCost[v];
    highestCost = std::max(highestCost, ccCost);
  }

  return highestCost <= 2.0 / 3.0;
}

void PlanarSeparator::findSeparatorFromComponents(
    NetworKit::ConnectedComponents &components) {
  (void)components;
  // No single component exceeds 2/3 of the total cost, so an empty separator
  // suffices: distribute whole components across the two sides.
  assignComponentsToSides(std::unordered_set<NetworKit::node>{});
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
