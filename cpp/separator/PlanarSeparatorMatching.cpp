#include "separator/PlanarSeparatorMatching.hpp"
#include "matching/MaximumMatching.hpp"
#include "networkit/Globals.hpp"
#include "networkit/graph/EdgeUtils.hpp"
#include "networkit/graph/Graph.hpp"
#include "separator/GraphUtils.hpp"
#include "separator/MISP.hpp"
#include <queue>

namespace {
std::vector<NetworKit::Edge> getExactMatching(NetworKit::Graph &graph) {
  Koala::EdmondsMaximumMatching m(graph, false);
  m.run();
  return m.getMatching();
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
  std::vector<std::unordered_set<NetworKit::node>> adj(n);

  graph.forNodes([&](NetworKit::node v) {
    graph.forNeighborsOf(v, [&](NetworKit::node u) { adj[v].insert(u); });
  });

  graph.forNodes([&](NetworKit::node v) {
    if (adj[v].size() <= 2)
      Q.push(v);
  });

  std::unordered_set<NetworKit::Edge> S;
  std::unordered_set<NetworKit::node> matched_nodes;

  bool stop = false;
  while (!stop) {
    if (currentGraphSize <= loglog) {
      NetworKit::Graph induced(adj.size());
      for (size_t i = 0; i < adj.size(); ++i) {
        if (adj[i].size() > 0) {
          for (auto u : adj[i]) {
            if (i < u)
              induced.addEdge(i, u);
          }
        }
      }
      auto matching = getExactMatching(induced);
      for (NetworKit::Edge e : matching) {
        S.insert(e);
        matched_nodes.insert(e.u);
        matched_nodes.insert(e.v);
      }
      stop = true;
    } else if (!Q.empty()) {
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
          if (adj[x].size() <= 2)
            Q.push(x);
        }
        adj[v].clear();
        adj[u].clear();
        currentGraphSize -= 2;
      } else if (adj[v].size() == 2) {
        auto it = adj[v].begin();
        NetworKit::node u = *it;
        NetworKit::node w = *std::next(it);
        NetworKit::node x = adj.size();
        adj.push_back({});
        std::unordered_set<NetworKit::node> seen;
        for (auto t : adj[u]) {
          if (t != v) {
            adj[t].erase(u);
            adj[t].insert(x);
            if (seen.find(t) == seen.end()) {
              adj[x].insert(t);
              seen.insert(t);
            }
          }
        }

        for (auto t : adj[w]) {
          if (t != v) {
            adj[t].erase(w);
            adj[t].insert(x);
            if (seen.find(t) == seen.end()) {
              adj[x].insert(t);
              seen.insert(t);
            }
          }
        }
        adj[v].clear();
        adj[u].clear();
        adj[w].clear();
        if (adj[x].size() <= 2)
          Q.push(x);

        stk.push_back({2, v, u, w});
        currentGraphSize -= 2;
      }
    } else {
      NetworKit::Graph induced(adj.size());
      for (size_t i = 0; i < adj.size(); ++i) {
        if (adj[i].size() > 0) {
          for (auto u : adj[i]) {
            if (i < u)
              induced.addEdge(i, u);
          }
        }
      }
      double epsilon = loglog;
      MISP<NetworKit::Edge> mispAlgo(
          induced, epsilon,
          [&](const NetworKit::Graph &c) -> std::vector<NetworKit::Edge> {
            auto copyC = c;
            auto matching = getExactMatching(copyC);
            return matching;
          });
      mispAlgo.run();
      auto matching = mispAlgo.maximum_independent_set;
      for (NetworKit::Edge e : matching) {
        S.insert(e);
        matched_nodes.insert(e.v);
        matched_nodes.insert(e.u);
      }
      stop = true;
    }
  }

  while (!stk.empty()) {
    auto action = stk.back();
    stk.pop_back();

    if (action.deg == 0) {
      // just skip
    } else if (action.deg == 1) {
      S.insert(NetworKit::Edge(action.v, action.u));
      matched_nodes.insert(action.v);
      matched_nodes.insert(action.u);
    } else if (action.deg == 2) {
      if (matched_nodes.count(action.u) == 0) {
        S.insert(NetworKit::Edge(action.v, action.u));
        matched_nodes.insert(action.v);
        matched_nodes.insert(action.u);
      } else {
        S.insert(NetworKit::Edge(action.v, action.w));
        matched_nodes.insert(action.v);
        matched_nodes.insert(action.w);
      }
    }
  }
  matching_set = std::vector<NetworKit::Edge>(S.begin(), S.end());
  hasRun = true;
}

} // namespace Koala
