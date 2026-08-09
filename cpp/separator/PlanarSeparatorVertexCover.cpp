#include "separator/PlanarSeparatorVertexCover.hpp"
#include "networkit/Globals.hpp"
#include "networkit/graph/EdgeUtils.hpp"
#include "networkit/graph/Graph.hpp"
#include <vector>

namespace Koala {
PlanarSeparatorVertexCover::PlanarSeparatorVertexCover(NetworKit::Graph &G)
    : graph(G) {}

void PlanarSeparatorVertexCover::run() {
  std::vector<bool> VC(graph.upperNodeIdBound());
  std::vector<bool> U(graph.upperNodeIdBound());
  bool stop = false;

  std::vector<std::vector<int>> adj;

  while (!stop) {
    std::vector<int> deg_residual(graph.upperNodeIdBound());
    graph.forNodes(
        [&](NetworKit::node v) { deg_residual[v] = graph.degree(v); });
    std::queue<NetworKit::node> Q;
    graph.forNodes([&](NetworKit::node v) {
      if (deg_residual[v] < 2)
        Q.push(v);
    });

    while (!Q.empty()) {
      auto v = Q.back();
      Q.pop();

      if (deg_residual[v] == 0) {
        U[v] = true;
      }

      if (deg_residual[v] == 1) {
        auto u = graph.getIthNeighbor(v, 0);
        U[v] = true;
        U[u] = true;
        VC[u] = true;

        graph.forNeighborsOf(u, [&](NetworKit::node t) {
          if (--deg_residual[t] < 2)
            Q.push(t);
        });
      }
    }

    // here deg(v) >= 2 for each v
    std::vector<bool> vis(graph.upperNodeIdBound(), true);
    graph.forNodes([&](NetworKit::node t) {
      if (!U[t])
        vis[t] = false;
    });

    for (int i = 0; i < static_cast<int>(vis.size()); i++) {
      if (!vis[i]) {
        std::vector<std::pair<NetworKit::node, int>> seen_deg_pair;
        Q.push(i);

        while (!Q.empty()) {
          auto v = Q.front();
          Q.pop();

          vis[v] = true;

          seen_deg_pair.push_back({v, deg_residual[v]});

          graph.forNeighborsOf(v, [&](NetworKit::node t) {
            if (!(vis[t] || U[t])) {
              Q.push(t);
            }
          });
        }

        bool all_degs_two = true;
        for (auto [v, deg] : seen_deg_pair) {
          if (deg != 2) {
            all_degs_two = false;
            break;
          }
        }

        if (all_degs_two) {
          auto first = seen_deg_pair.front().first;
          auto cur = first;
          auto prev = NetworKit::none;
          bool go = true;
          int i = 1;

          while (go) {
            auto neighborRange = graph.neighborRange(cur);
            NetworKit::node next = NetworKit::none;
            for (auto v : neighborRange) {
              if (v != prev) {
                next = v;
                break;
              }
            }

            if (i++ % 2 == 1)
              VC[cur] = true;
            U[cur] = true;

            prev = cur;
            cur = next;
            if (cur == first)
              go = false;
          }
        }
      }
      // bipartite finding
      std::vector<int> label(graph.upperNodeIdBound(), 0);
      graph.forNodes([&](NetworKit::node t) {
        if (!U[t] && deg_residual[t] > 2) {
          label[t] = '+';
          Q.push(t);
        }
      });
      while (!Q.empty()) {
        auto v = Q.front();
        Q.pop();

        graph.forNeighborsOf(v, [&](NetworKit::node t) {
          if (!U[t] && label[t] == 0) {
            if (label[v] == '+')
              label[t] = '-';
            else
              label[t] = '+';

            Q.push(t);
          }
        });
      }

      std::vector<bool> isInX(graph.upperNodeIdBound(), false);
      std::vector<bool> isInY(graph.upperNodeIdBound(), false);
      std::vector<int> degB(graph.upperNodeIdBound(), 0);

      for (int i = 0; i < static_cast<int>(label.size()); i++) {
        if (label[i] == '+') {
          isInX[i] = true;
          graph.forNeighborsOf(i, [&](NetworKit::node t) {
            if (label[t] == '-') {
              degB[i]++;
              degB[t]++;
            }
          });
        }
        if (label[i] == '-')
          isInY[i] = true;
      }

      graph.forNodes([&](NetworKit::node v) {
        if (!U[v] && degB[v] < 2)
          Q.push(v);
      });

      while (!Q.empty()) {
        auto v = Q.front();
        Q.pop();

        isInX[v] = false;
        isInY[v] = false;
        auto u = graph.getIthNeighbor(v, 0);

        graph.forNeighborsOf(u, [&](NetworKit::node t) {
          if ((isInX[t] || isInY[t]) && --degB[t] < 2)
            Q.push(t);
        });
      }
    }
  }

  hasRun = true;
}

void PlanarSeparatorVertexCover::prep() {}
} // namespace Koala
