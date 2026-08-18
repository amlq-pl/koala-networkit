#include "separator/PlanarSeparatorMatching.hpp"
#include "matching/MaximumMatching.hpp"
#include "networkit/Globals.hpp"
#include "networkit/graph/Graph.hpp"
#include <queue>

namespace {
std::vector<NetworKit::Edge> getExactMatching(NetworKit::Graph &graph) {
  EdmondsMaximumMatching m(graph, false);
  m.run();
}

void removeVertex();
void removeEdge();
void contract();
} // namespace

namespace Koala {
PlanarSeparatorMatching::PlanarSeparatorMatching(NetworKit::Graph &G)
    : graph(G) {}

void PlanarSeparatorMatching::run() { hasRun = true; }

struct Action {
  int deg;
  NetworKit::node v = 0, u = 0, w = 0;
};

std::vector<NetworKit::Edge>
PlanarSeparatorMatching::reduce_procedure(NetworKit::Graph &graph) {
  NetworKit::count n = graph.numberOfNodes();
  double loglog =
      std::max(1.0, std::log2(std::log2((double)std::max<size_t>(4, n))));

  int currentGraphSize = n;
  std::vector<Action> stk;
  std::queue<NetworKit::node> Q;
  std::vector<std::unordered_set<NetworKit::node>> adj;

  graph.forNodes([&](NetworKit::node v) {
    graph.forNeighborsOf(v, [&](NetworKit::node u) { adj[v].insert(u); });
  });

  graph.forNodes([&](NetworKit::node v) {
    if (adj[v].size() <= 2)
      Q.push(v);
  });

  while (!Q.empty()) {
    if (currentGraphSize <= loglog) {
      getExactMatching(adj);
      break;
    }
    auto v = Q.front();
    Q.pop();

    if (adj[v].size() > 2)
      continue;

    if (adj[v].size() == 0) {
      stk.push_back({0, v});
      currentGraphSize--;
    } else if (adj[v].size() == 1) {
      NetworKit::node u = *adj[v].begin();
      stk.push_back({1, v, u});
      for (auto x : adj[u]) {
        adj[x].erase(u);
      }
      adj[u].clear();
      currentGraphSize -= 2;
    } else {
      auto it = adj[v].begin();
      NetworKit::node u = *it;
      NetworKit::node w = *std::next(it);
      NetworKit::node x = adj.size();
      std::unordered_set<NetworKit::node> seen;
      for (auto t : adj[u]) {
        if (t != v) {
          adj[t].erase(u);
          adj[t].insert(x);
          if (seen.find(t) == seen.end()) {
            adj[x].insert(t);
          }
        }
      }

      for (auto t : adj[w]) {
        if (t != v) {
          adj[t].erase(w);
          adj[t].insert(x);
          if (seen.find(t) == seen.end()) {
            adj[x].insert(t);
          }
        }
      }
      adj[v].clear();
      adj[u].clear();
      adj[w].clear();

      stk.push_back({2, v, u, w});
      currentGraphSize -= 2;
    }
  }

  while (!stk.empty()) {
    auto action = stk.back();
    stk.pop_back();

    if (action.deg == 0) {

    } else if (action.deg == 1) {

    } else if (action.deg == 2) {
    }
  }
}

} // namespace Koala
