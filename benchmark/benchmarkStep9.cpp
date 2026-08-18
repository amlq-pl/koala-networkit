// Isolated time-complexity benchmark for Step 9 of Koala::PlanarSeparator
// (Koala::shrinkFundamentalCycle).
//
// The full PlanarSeparator::run() is dominated (>95%) by the Boyer-Myrvold
// planar-embedding routine, which is empirically super-linear. That noise makes
// it impossible to see the complexity of the fundamental-cycle shrink (Step 9)
// from the end-to-end benchmark. Here we compute the embedding, spanning tree
// and subtree costs ONCE, OUTSIDE the timed region, and time only the shrink:
//
//     p = log(t / t_prev) / log(n / n_prev)
//
// If Step 9 is truly O(n), p should hover around 1.
//
// Usage:
//     benchmark_step9 [startN] [maxN] [reps]
//
//   startN   first target vertex count            (default 1000)
//   maxN     hard cap on target vertex count       (default 4000000)
//   reps     timed repetitions per size (min)      (default 5)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <graph/PlanarGraphTools.hpp>
#include <separator/PlanarSeparatorStep9.hpp>

using NetworKit::node;
using NetworKit::none;

namespace {

// Deterministic triangulated square grid on side*side vertices: every unit
// square is split by its (r,c)-(r+1,c+1) diagonal.
NetworKit::Graph buildTriangulatedGrid(int side) {
  auto idx = [side](int r, int c) { return r * side + c; };
  NetworKit::Graph G(side * side, false, false);
  for (int r = 0; r < side; r++)
    for (int c = 0; c < side; c++) {
      if (c + 1 < side)
        G.addEdge(idx(r, c), idx(r, c + 1));
      if (r + 1 < side)
        G.addEdge(idx(r, c), idx(r + 1, c));
      if (r + 1 < side && c + 1 < side)
        G.addEdge(idx(r, c), idx(r + 1, c + 1));
    }
  return G;
}

// BFS spanning tree of H from `root`, filling parent[] and depth[].
void bfsTree(const NetworKit::Graph &H, node root, std::vector<node> &parent,
             std::vector<node> &depth) {
  node nb = H.upperNodeIdBound();
  parent.assign(nb, none);
  depth.assign(nb, 0);
  std::vector<char> seen(nb, 0);
  std::queue<node> q;
  seen[root] = 1;
  q.push(root);
  while (!q.empty()) {
    node v = q.front();
    q.pop();
    H.forNeighborsOf(v, [&](node w) {
      if (!seen[w]) {
        seen[w] = 1;
        parent[w] = v;
        depth[w] = depth[v] + 1;
        q.push(w);
      }
    });
  }
}

// Subtree costs, mirroring PlanarSeparator's Step 7: the root plays the role of
// the contracted super-vertex x (cost 0), every other vertex has cost 1, so
// costs[root] equals the total interior weight (n - 1).
std::vector<double> subtreeCosts(const NetworKit::Graph &H, node root,
                                 const std::vector<node> &parent,
                                 const std::vector<double> &vcost) {
  node nb = H.upperNodeIdBound();
  std::vector<std::vector<node>> children(nb);
  H.forNodes([&](node v) {
    if (parent[v] != none)
      children[parent[v]].push_back(v);
  });
  std::vector<double> costs(nb, 0.0);
  std::vector<std::pair<node, size_t>> stk;
  stk.emplace_back(root, 0);
  while (!stk.empty()) {
    node v = stk.back().first;
    size_t &i = stk.back().second;
    if (i < children[v].size()) {
      node c = children[v][i++];
      costs[c] = vcost[c];
      stk.emplace_back(c, 0);
    } else {
      node f = v;
      stk.pop_back();
      if (!stk.empty())
        costs[stk.back().first] += costs[f];
    }
  }
  return costs;
}

// Vertices of the fundamental cycle of non-tree edge (u, v).
std::vector<node> fundamentalCycle(node u, node v,
                                   const std::vector<node> &parent,
                                   const std::vector<node> &depth) {
  std::vector<node> pu, pv;
  node a = u, b = v;
  while (depth[a] > depth[b]) {
    pu.push_back(a);
    a = parent[a];
  }
  while (depth[b] > depth[a]) {
    pv.push_back(b);
    b = parent[b];
  }
  while (a != b) {
    pu.push_back(a);
    a = parent[a];
    pv.push_back(b);
    b = parent[b];
  }
  pu.push_back(a); // LCA
  pu.insert(pu.end(), pv.rbegin(), pv.rend());
  return pu;
}

// Number of vertices in the heavier of the two regions the fundamental cycle of
// (u, v) cuts H into -- exactly the count the shrink starts with "inside", so it
// tells us how much work the shrink will do.
size_t heavierSide(const NetworKit::Graph &H, node u, node v,
                   const std::vector<node> &parent,
                   const std::vector<node> &depth) {
  node nb = H.upperNodeIdBound();
  auto cyc = fundamentalCycle(u, v, parent, depth);
  std::vector<char> onCycle(nb, 0);
  for (node c : cyc)
    onCycle[c] = 1;
  node seed = none;
  H.forNodes([&](node w) {
    if (seed == none && !onCycle[w])
      seed = w;
  });
  if (seed == none)
    return 0;
  std::vector<char> seen(nb, 0);
  std::queue<node> q;
  seen[seed] = 1;
  q.push(seed);
  size_t reached = 0;
  while (!q.empty()) {
    node w = q.front();
    q.pop();
    reached++;
    H.forNeighborsOf(w, [&](node z) {
      if (!seen[z] && !onCycle[z]) {
        seen[z] = 1;
        q.push(z);
      }
    });
  }
  size_t other = H.numberOfNodes() - cyc.size() - reached;
  return std::max(reached, other);
}

long long parseArg(int argc, char **argv, int index, long long fallback) {
  return index < argc ? std::strtoll(argv[index], nullptr, 10) : fallback;
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 1 &&
      (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
    std::cerr << "Usage: " << argv[0] << " [startN] [maxN] [reps]\n"
              << "  Times ONLY Step 9 (shrinkFundamentalCycle) with a\n"
              << "  precomputed embedding and prints n, m, work and exp_p.\n";
    return 0;
  }

  const long long startN = parseArg(argc, argv, 1, 1000);
  const long long maxN = parseArg(argc, argv, 2, 4000000);
  const int reps = static_cast<int>(parseArg(argc, argv, 3, 5));
  const double growth = 2.0;

  std::cout << "# Step 9 (shrinkFundamentalCycle) time-complexity benchmark\n";
  std::cout << "# embedding/tree/costs are precomputed OUTSIDE the timed loop\n";
  std::cout << "#" << std::setw(11) << "n" << std::setw(13) << "m"
            << std::setw(12) << "inside" << std::setw(12) << "thresh"
            << std::setw(12) << "work" << std::setw(14) << "time_s"
            << std::setw(9) << "exp_p" << "\n";

  double prevTime = 0.0;
  long long prevN = 0;

  for (long long targetN = startN; targetN <= maxN;
       targetN = static_cast<long long>(std::ceil(targetN * growth))) {
    int side = static_cast<int>(std::llround(std::sqrt(
        static_cast<double>(targetN))));
    if (side < 3)
      side = 3;
    if (side > 40000)
      side = 40000;

    // ---- All setup below is UNTIMED. ----
    NetworKit::Graph G = buildTriangulatedGrid(side);
    NetworKit::Graph H = Koala::PlanarGraphTools::makeMaximalPlanar(G);
    auto embeddingH = Koala::PlanarGraphTools::findPlanarEmbedding(H);

    node nb = H.upperNodeIdBound();
    std::vector<std::unordered_map<node, int>> idxOf(nb);
    H.forNodes([&](node u) {
      const auto &rot = embeddingH[u];
      for (int i = 0; i < static_cast<int>(rot.size()); ++i)
        idxOf[u][rot[i]] = i;
    });

    const node root = 0;
    std::vector<node> parent, depth;
    bfsTree(H, root, parent, depth);
    std::vector<double> vcost(nb, 1.0);
    vcost[root] = 0.0; // root == super-vertex x, cost 0
    std::vector<double> costs = subtreeCosts(H, root, parent, vcost);

    // Pick the non-tree edge whose fundamental cycle encloses the most
    // vertices, so the shrink actually peels Theta(n) of them. Rank all
    // non-tree edges cheaply by tree-depth sum, then flood-score the top few.
    std::vector<std::pair<node, node>> cands;
    H.forEdges([&](node a, node b) {
      if (parent[a] != b && parent[b] != a)
        cands.emplace_back(a, b);
    });
    std::sort(cands.begin(), cands.end(), [&](const auto &e1, const auto &e2) {
      return depth[e1.first] + depth[e1.second] >
             depth[e2.first] + depth[e2.second];
    });
    const size_t K = std::min<size_t>(cands.size(), 40);
    node v1 = none, w1 = none;
    size_t inside = 0;
    for (size_t i = 0; i < K; i++) {
      size_t h = heavierSide(H, cands[i].first, cands[i].second, parent, depth);
      if (h > inside) {
        inside = h;
        v1 = cands[i].first;
        w1 = cands[i].second;
      }
    }
    if (v1 == none)
      continue;

    // ---- Timed region: only Step 9. ----
    double bestT = std::numeric_limits<double>::infinity();
    for (int r = 0; r < reps; r++) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto cycle = Koala::shrinkFundamentalCycle(H, embeddingH, idxOf, parent,
                                                 costs, vcost, root, v1, w1);
      auto t1 = std::chrono::high_resolution_clock::now();
      bestT = std::min(bestT, std::chrono::duration<double>(t1 - t0).count());
      (void)cycle;
    }

    const long long n = static_cast<long long>(H.numberOfNodes());
    const long long m = static_cast<long long>(H.numberOfEdges());
    const double thresh = 2.0 / 3.0 * costs[root];
    const double work = std::max(0.0, static_cast<double>(inside) - thresh);

    double expP = 0.0;
    if (prevN > 0 && prevTime > 0.0 && bestT > 0.0)
      expP = std::log(bestT / prevTime) /
             std::log(static_cast<double>(n) / static_cast<double>(prevN));

    std::cout << std::setw(12) << n << std::setw(13) << m << std::setw(12)
              << inside << std::setw(12) << std::fixed << std::setprecision(0)
              << thresh << std::setw(12) << work << std::setw(14)
              << std::setprecision(6) << bestT << std::setw(9)
              << std::setprecision(3) << expP << "\n"
              << std::flush;

    prevTime = bestT;
    prevN = n;
  }

  return 0;
}
