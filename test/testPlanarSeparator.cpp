#include <gtest/gtest.h>
#include <separator/PlanarSeparator.hpp>
#include <unordered_set>
#include <cmath>
#include <list>
#include <random>
#include <tuple>
#include <string>
#include "helpers.hpp"

// ---------------------------------------------------------------------------
// Helper: verify the four separator properties.
//   1. S, A, B are pairwise disjoint and cover V(G).
//   2. No edge has one endpoint in A and the other in B.
//   3. |A| <= 2n/3 and |B| <= 2n/3.
//   4. |S| is "reasonable" — we use a generous constant times sqrt(n).
// ---------------------------------------------------------------------------
static void verifySeparator(
        const NetworKit::Graph& G,
        const std::vector<NetworKit::node>& sep,
        const std::vector<NetworKit::node>& partA,
        const std::vector<NetworKit::node>& partB,
        double sizeConstant = 8.0) {
    const NetworKit::count n = G.numberOfNodes();

    std::unordered_set<NetworKit::node> S(sep.begin(), sep.end());
    std::unordered_set<NetworKit::node> A(partA.begin(), partA.end());
    std::unordered_set<NetworKit::node> B(partB.begin(), partB.end());

    // Disjointness.
    for (auto v : S) {
        EXPECT_EQ(A.count(v), 0u) << "Vertex " << v << " is in both S and A";
        EXPECT_EQ(B.count(v), 0u) << "Vertex " << v << " is in both S and B";
    }
    for (auto v : A) {
        EXPECT_EQ(B.count(v), 0u) << "Vertex " << v << " is in both A and B";
    }

    // Coverage.
    EXPECT_EQ(S.size() + A.size() + B.size(), n)
        << "Sizes do not sum to n: |S|=" << S.size()
        << ", |A|=" << A.size() << ", |B|=" << B.size() << ", n=" << n;

    G.forNodes([&](NetworKit::node v) {
        EXPECT_TRUE(S.count(v) || A.count(v) || B.count(v))
            << "Vertex " << v << " is missing from S, A, and B";
    });

    // No A-B edges.
    G.forEdges([&](NetworKit::node u, NetworKit::node v) {
        const bool inA_u = A.count(u), inB_u = B.count(u);
        const bool inA_v = A.count(v), inB_v = B.count(v);
        const bool crosses = (inA_u && inB_v) || (inB_u && inA_v);
        EXPECT_FALSE(crosses) << "Edge (" << u << "," << v
                              << ") crosses partitions A and B";
    });

    // Balance bound.
    const NetworKit::count maxSide = (2 * n + 2) / 3;  // ceil(2n/3)
    EXPECT_LE(A.size(), maxSide) << "|A| exceeds 2n/3";
    EXPECT_LE(B.size(), maxSide) << "|B| exceeds 2n/3";

    // Size bound on S — generous, mostly a sanity check.
    const double bound =
        sizeConstant * std::sqrt(static_cast<double>(std::max<NetworKit::count>(n, 1)));
    EXPECT_LE(static_cast<double>(S.size()), bound + 1.0)
        << "Separator size " << S.size() << " exceeds " << bound;
}

// ---------------------------------------------------------------------------
struct SeparatorParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
};

class PlanarSeparatorTest : public testing::TestWithParam<SeparatorParameters> {};

TEST_P(PlanarSeparatorTest, TestSeparatorProperties) {
    const SeparatorParameters& parameters = GetParam();
    std::cout << "\nTest name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);

    Koala::PlanarSeparator algorithm(G);
    algorithm.run();

    verifySeparator(G, algorithm.getSeparator(),
                    algorithm.getPartitionA(), algorithm.getPartitionB());
}

// ---------------------------------------------------------------------------
// Helpers to construct standard families of planar graphs.
// ---------------------------------------------------------------------------
static std::list<std::tuple<int, int, int>> makeGridEdges(int rows, int cols) {
    std::list<std::tuple<int, int, int>> E;
    auto idx = [cols](int r, int c) { return r * cols + c; };
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (c + 1 < cols) E.emplace_back(idx(r, c), idx(r, c + 1), 1);
            if (r + 1 < rows) E.emplace_back(idx(r, c), idx(r + 1, c), 1);
        }
    }
    return E;
}

static std::list<std::tuple<int, int, int>> makePathEdges(int n) {
    std::list<std::tuple<int, int, int>> E;
    for (int i = 0; i + 1 < n; i++) E.emplace_back(i, i + 1, 1);
    return E;
}

static std::list<std::tuple<int, int, int>> makeCycleEdges(int n) {
    auto E = makePathEdges(n);
    if (n > 1) E.emplace_back(n - 1, 0, 1);
    return E;
}

static std::list<std::tuple<int, int, int>> makeWheelEdges(int rim) {
    // vertex 0 is the hub, 1..rim form the rim cycle
    std::list<std::tuple<int, int, int>> E;
    for (int i = 1; i <= rim; i++) E.emplace_back(0, i, 1);
    for (int i = 1; i < rim; i++) E.emplace_back(i, i + 1, 1);
    E.emplace_back(rim, 1, 1);
    return E;
}

static std::list<std::tuple<int, int, int>> makeTwoTriangles() {
    return {
        {0, 1, 1}, {1, 2, 1}, {2, 0, 1},
        {3, 4, 1}, {4, 5, 1}, {5, 3, 1}
    };
}

INSTANTIATE_TEST_SUITE_P(
    PlanarSeparatorTests,
    PlanarSeparatorTest,
    testing::Values(
        SeparatorParameters{"path_5",         5, makePathEdges(5)},
        SeparatorParameters{"path_10",       10, makePathEdges(10)},
        SeparatorParameters{"cycle_6",        6, makeCycleEdges(6)},
        SeparatorParameters{"cycle_12",      12, makeCycleEdges(12)},
        SeparatorParameters{"triangle",       3, makeCycleEdges(3)},
        SeparatorParameters{"grid_3x3",       9, makeGridEdges(3, 3)},
        SeparatorParameters{"grid_4x4",      16, makeGridEdges(4, 4)},
        SeparatorParameters{"grid_5x5",      25, makeGridEdges(5, 5)},
        SeparatorParameters{"grid_6x6",      36, makeGridEdges(6, 6)},
        SeparatorParameters{"grid_3x7",      21, makeGridEdges(3, 7)},
        SeparatorParameters{"wheel_6",        7, makeWheelEdges(6)},
        SeparatorParameters{"wheel_12",      13, makeWheelEdges(12)},
        SeparatorParameters{"two_triangles",  6, makeTwoTriangles()},
        SeparatorParameters{"triangle_prism", 6, {
            {0, 1, 1}, {1, 2, 1}, {2, 0, 1},
            {3, 4, 1}, {4, 5, 1}, {5, 3, 1},
            {0, 3, 1}, {1, 4, 1}, {2, 5, 1}
        }}
    )
);

// ---------------------------------------------------------------------------
// Large randomized test: triangulated grid (planar) with random vertex
// deletions. The resulting graph is still planar; nodes are relabelled to
// be contiguous starting at 0.
// ---------------------------------------------------------------------------
static NetworKit::Graph buildTriangulatedGridMinusRandomNodes(
        int rows, int cols, double deletionRate, uint64_t seed) {
    // Build the full triangulated grid: each cell gets one diagonal.
    // Vertex (r, c) has id r * cols + c.
    const int total = rows * cols;
    auto idx = [cols](int r, int c) { return r * cols + c; };

    std::vector<std::pair<int, int>> edges;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (c + 1 < cols) edges.emplace_back(idx(r, c), idx(r, c + 1));
            if (r + 1 < rows) edges.emplace_back(idx(r, c), idx(r + 1, c));
            if (r + 1 < rows && c + 1 < cols) {
                // Diagonal: alternate orientation for variety.
                if ((r + c) % 2 == 0) {
                    edges.emplace_back(idx(r, c), idx(r + 1, c + 1));
                } else {
                    edges.emplace_back(idx(r, c + 1), idx(r + 1, c));
                }
            }
        }
    }

    // Randomly mark vertices for deletion.
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<bool> deleted(total, false);
    for (int v = 0; v < total; v++) {
        if (dist(rng) < deletionRate) deleted[v] = true;
    }

    // Relabel surviving vertices to a contiguous 0..N-1 range.
    std::vector<int> newId(total, -1);
    int next = 0;
    for (int v = 0; v < total; v++) {
        if (!deleted[v]) newId[v] = next++;
    }
    const int N = next;

    NetworKit::Graph G(N, false, false);
    for (auto [u, v] : edges) {
        if (!deleted[u] && !deleted[v]) {
            G.addEdge(newId[u], newId[v]);
        }
    }
    return G;
}

TEST(PlanarSeparatorLargeTest, TriangulatedGridMinusRandom) {
    NetworKit::Graph G = buildTriangulatedGridMinusRandomNodes(
        15, 15, /*deletionRate=*/0.10, /*seed=*/12345);
    std::cout << "\nLarge test: n=" << G.numberOfNodes()
              << ", m=" << G.numberOfEdges() << std::endl;

    Koala::PlanarSeparator algorithm(G);
    algorithm.run();

    verifySeparator(G, algorithm.getSeparator(),
                    algorithm.getPartitionA(), algorithm.getPartitionB());

    std::cout << "  |S|=" << algorithm.getSeparator().size()
              << "  |A|=" << algorithm.getPartitionA().size()
              << "  |B|=" << algorithm.getPartitionB().size() << std::endl;
}

TEST(PlanarSeparatorLargeTest, TriangulatedGrid200ish) {
    // ~200-node triangulated grid (14x14 = 196), no deletions.
    NetworKit::Graph G = buildTriangulatedGridMinusRandomNodes(
        300, 100, /*deletionRate=*/0.0, /*seed=*/777);
    std::cout << "\nLarge test (no deletions): n=" << G.numberOfNodes()
              << ", m=" << G.numberOfEdges() << std::endl;

    Koala::PlanarSeparator algorithm(G);
    algorithm.run();

    verifySeparator(G, algorithm.getSeparator(),
                    algorithm.getPartitionA(), algorithm.getPartitionB());

    std::cout << "  |S|=" << algorithm.getSeparator().size()
              << "  |A|=" << algorithm.getPartitionA().size()
              << "  |B|=" << algorithm.getPartitionB().size() << std::endl;
}
