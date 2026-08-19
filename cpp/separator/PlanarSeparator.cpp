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

int wrapIndex(int i, int size) { return ((i % size) + size) % size; }

NetworKit::node findNodeInCycleAtIthPosition(cycle_t &cycle, int i) {
  return cycle[wrapIndex(i, static_cast<int>(cycle.size()))];
}

CostComputation
computeSidesCost(cycle_t &cycle,
                 std::vector<std::unordered_map<NetworKit::node, int>> &idxOf,
                 planar_embedding_t &embedding,
                 std::vector<NetworKit::node> &parent,
                 std::vector<double> &costs, NetworKit::node root) {
  double arcTrueCost = 0.0;
  double arcFalseCost = 0.0;

  int cycleSize = cycle.size();

  for (int i = 0; i < cycleSize; i++) {
    auto prev = findNodeInCycleAtIthPosition(cycle, i - 1);
    auto v = findNodeInCycleAtIthPosition(cycle, i);
    auto u = findNodeInCycleAtIthPosition(cycle, i + 1);

    auto neighborsV = embedding[v];
    int posU = idxOf[v][u];
    int posPrev = idxOf[v][prev];

    for (int j = 0; j < static_cast<int>(neighborsV.size()); j++) {
      auto w = neighborsV[j];
      if (w == prev || w == u)
        continue;

      double cost = 0.0;
      if (parent[w] == v) {
        cost = costs[w];
      } else if (parent[v] == w) {
        cost = costs[root] - costs[v];
      } else {
        continue;
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

std::optional<cycle_t> shrinkFundamentalCycle(
    NetworKit::Graph &H, planar_embedding_t &embeddingH,
    std::vector<std::unordered_map<NetworKit::node, int>> &idxOf,
    std::vector<NetworKit::node> &parentH, std::vector<double> &costsH,
    std::vector<double> &vertexCost, NetworKit::node x, NetworKit::node v1,
    NetworKit::node w1) {
  auto cycleOpt = buildFundamentalCycle(v1, w1, parentH);

  if (!cycleOpt) {
    return std::nullopt;
  }

  auto cycle = cycleOpt.value();

  if (cycle.size() < 3) {
    return std::nullopt;
  }

  auto cc = computeSidesCost(cycle, idxOf, embeddingH, parentH, costsH, x);

  // here we enter step 9
  cycle_t cycle_i = std::move(cycle);

  int dir = cc.insideIsClockwiseArc ? -1 : 1;
  auto getApex = [&](NetworKit::node v, NetworKit::node u) {
    const auto &neighborsV = embeddingH[v];
    int deg = static_cast<int>(neighborsV.size());
    int j = idxOf[v][u];
    return neighborsV[wrapIndex(j + dir, deg)];
  };

  NetworKit::count nodeBound = H.upperNodeIdBound();
  auto vcost = [&](NetworKit::node v) { return v == x ? 0.0 : vertexCost[v]; };
  double threshold = 2.0 / 3.0 * costsH[x];

  std::vector<int> is_on_cycle(nodeBound, 0);
  for (auto c : cycle_i) {
    is_on_cycle[c] = 1;
  }

  std::vector<NetworKit::node> cycleNext(nodeBound, NetworKit::none);
  std::vector<NetworKit::node> cyclePrev(nodeBound, NetworKit::none);
  for (size_t i = 0; i < cycle_i.size(); i++) {
    cycleNext[cycle_i[i]] = cycle_i[wrapIndex(i + 1, cycle_i.size())];
    cyclePrev[cycle_i[i]] = cycle_i[wrapIndex(i - 1, cycle_i.size())];
  }

  NetworKit::node curVi = v1;
  NetworKit::node curWi = w1;
  double insideCost = cc.insideCost;

  bool forward = (cycleNext[curVi] == curWi);
  auto cNext = [&](NetworKit::node v) {
    return forward ? cycleNext[v] : cyclePrev[v];
  };

  long iterGuard = 0;
  const long iterCap = 8 * static_cast<long>(nodeBound);

  while (insideCost > threshold) {
    if (++iterGuard > iterCap)
      return std::nullopt;

    NetworKit::node y = getApex(curVi, curWi);

    std::vector<NetworKit::node> pPath;
    NetworKit::node z = y;
    while (z != NetworKit::none && !is_on_cycle[z]) {
      pPath.push_back(z);
      z = parentH[z];
    }
    if (z == NetworKit::none)
      break;
    pPath.push_back(z);

    std::vector<NetworKit::node> seqV;
    NetworKit::node curr = z;
    while (curr != curVi) {
      seqV.push_back(curr);
      curr = cNext(curr);
    }
    seqV.push_back(curVi);
    for (size_t i = 0; i + 1 < pPath.size(); i++) {
      seqV.push_back(pPath[i]);
    }

    std::vector<NetworKit::node> seqW;
    curr = curWi;
    while (curr != z) {
      seqW.push_back(curr);
      curr = cNext(curr);
    }
    seqW.push_back(z);
    for (int i = static_cast<int>(pPath.size()) - 2; i >= 0; i--) {
      seqW.push_back(pPath[i]);
    }

    double costV = 0.0, costW = 0.0;
    int idxV = 0, idxW = 0;
    bool doneV = false, doneW = false;

    auto stepSeq = [&](const std::vector<NetworKit::node> &seq, int &idx,
                       double &cost, bool &done) {
      if (idx >= static_cast<int>(seq.size())) {
        done = true;
        return;
      }
      int sz = seq.size();
      if (sz <= 2) {
        cost = 0;
        idx = sz;
        done = true;
        return;
      }

      NetworKit::node prv = seq[(idx - 1 + sz) % sz];
      NetworKit::node v = seq[idx];
      NetworKit::node nxt = seq[(idx + 1) % sz];

      int start = idxOf[v][nxt];
      int end = idxOf[v][prv];
      int deg = embeddingH[v].size();

      for (int step = 1; step < deg; step++) {
        int j = wrapIndex(start + dir * step, deg);
        if (j == end)
          break;
        NetworKit::node w = embeddingH[v][j];
        if (parentH[w] == v) {
          cost += costsH[w];
        }
      }
      idx++;
      if (idx >= sz)
        done = true;
    };

    while (!doneV && !doneW) {
      stepSeq(seqV, idxV, costV, doneV);
      stepSeq(seqW, idxW, costW, doneW);
    }

    double Pcost = 0.0;
    for (size_t i = 0; i + 1 < pPath.size(); i++) {
      Pcost += vcost(pPath[i]);
    }

    if (doneV) {
      costW = insideCost - Pcost - costV;
    } else {
      costV = insideCost - Pcost - costW;
    }

    bool keepV = (costV >= costW);
    insideCost = keepV ? costV : costW;

    const auto &discardedSeq = keepV ? seqW : seqV;
    const auto &keptSeq = keepV ? seqV : seqW;

    std::unordered_set<NetworKit::node> keptSet(keptSeq.begin(), keptSeq.end());
    for (NetworKit::node p : discardedSeq) {
      if (!keptSet.count(p)) {
        is_on_cycle[p] = 0;
      }
    }

    for (size_t i = 0; i < keptSeq.size(); i++) {
      NetworKit::node v = keptSeq[i];
      NetworKit::node prv = keptSeq[(i - 1 + keptSeq.size()) % keptSeq.size()];
      NetworKit::node nxt = keptSeq[(i + 1) % keptSeq.size()];
      if (forward) {
        cycleNext[v] = nxt;
        cyclePrev[v] = prv;
      } else {
        cyclePrev[v] = nxt;
        cycleNext[v] = prv;
      }
      is_on_cycle[v] = 1;
    }

    if (keepV) {
      curWi = y;
    } else {
      curVi = y;
    }
  }

  cycle_t finalCycle;
  NetworKit::node curr = curVi;
  do {
    finalCycle.push_back(curr);
    curr = cNext(curr);
  } while (curr != curVi);

  return finalCycle;
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

  if (graph.numberOfNodes() == 0) {
    hasRun = true;
    return;
  }

  // normalize total cost to 1
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

    NetworKit::node root = NetworKit::none;
    G.forNodes([&](NetworKit::node v) {
      if (root == NetworKit::none)
        root = v;
    });

    auto [lvl, parent] = performBFSAndFindSpanningTree(G, root);

    auto verticesAtLevel = findNumberOfVerticesAtLevel(lvl);

    std::vector<NetworKit::count> prefixSum(verticesAtLevel.size());
    std::partial_sum(verticesAtLevel.begin(), verticesAtLevel.end(),
                     prefixSum.begin());

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
    NetworKit::node l0 = 0;
    for (int i = static_cast<int>(l1); i >= 0; i--) {
      if (verticesAtLevel[i] + 2.0 * (static_cast<int>(l1) - i) <=
          2.0 * std::sqrt(static_cast<double>(k))) {
        l0 = static_cast<NetworKit::node>(i);
        break;
      }
    }
    // l2 is the smallest level > l1 such that
    // |L(l2)| + 2*(l2-l1-1) <= 2*sqrt(n-k)
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

    // Step 7 compute spanning tree of H and costs of subtrees hanging off it
    auto [lvlH, parentH] = performBFSAndFindSpanningTree(H, x);
    std::vector<std::vector<NetworKit::node>> childrenH(
        H.upperNodeIdBound(), std::vector<NetworKit::node>());
    H.forNodes([&](NetworKit::node v) {
      if (parentH[v] != NetworKit::none) {
        childrenH[parentH[v]].push_back(v);
      }
    });
    std::vector<double> costsH(H.upperNodeIdBound(), 0.0);
    {
      std::vector<std::pair<NetworKit::node, size_t>> stk;
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

    std::vector<std::unordered_map<NetworKit::node, int>> idxOf(
        H.upperNodeIdBound());
    H.forNodes([&](NetworKit::node u) {
      const auto &rotation = embeddingH[u];
      for (int i = 0; i < static_cast<int>(rotation.size()); ++i) {
        idxOf[u][rotation[i]] = i;
      }
    });

    // Step 8
    auto [v1, w1] = findNonTreeEdge(H, parentH);

    // Degenerate H (no non-tree edge): fall back to a level-based separator.
    if (v1 == NetworKit::none) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }

    // ---- Step 9: shrink the fundamental cycle to a balanced one (O(n)). ----
    auto finalCycleOpt = shrinkFundamentalCycle(H, embeddingH, idxOf, parentH,
                                                costsH, vertexCost, x, v1, w1);
    if (!finalCycleOpt || finalCycleOpt->size() < 3) {
      fallbackLevelSeparator(lvl, l1);
      hasRun = true;
      return;
    }
    cycle_t finalCycle = std::move(*finalCycleOpt);

    // ---- Step 10: extract separator and partitions back to G ----
    extractSeparatorAndPartitions(G, lvl, l0, l2, finalCycle, x);

    double cA = 0.0, cB = 0.0;
    for (auto v : partition.A)
      cA += vertexCost[v];
    for (auto v : partition.B)
      cB += vertexCost[v];
    if (cA > 2.0 / 3.0 || cB > 2.0 / 3.0) {
      cleanPartitions();
      fallbackLevelSeparator(lvl, l1);
    }
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

// if the cycle is degenerate, fallback to lvl1 as the separator
void PlanarSeparator::fallbackLevelSeparator(
    const std::vector<NetworKit::node> &lvl, NetworKit::node l1) {
  std::unordered_set<NetworKit::node> separatorNodes;
  for (NetworKit::node v = 0; v < lvl.size(); v++) {
    if (lvl[v] == l1 && graph.hasNode(v))
      separatorNodes.insert(v);
  }

  assignComponentsToSides(separatorNodes);
}

bool PlanarSeparator::areConnectedComponentsEligibleForPartition(
    NetworKit::ConnectedComponents &components) {
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
