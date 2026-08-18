#include <gtest/gtest.h>
#include <separator/PlanarSeparator.hpp>
#include <array>
#include <unordered_set>
#include <cmath>
#include <list>
#include <random>
#include <tuple>
#include <vector>
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

// ---------------------------------------------------------------------------
// Additional planar families, sizes and randomized graphs.
// ---------------------------------------------------------------------------

// Random maximal-planar graph (Apollonian network): start from a triangle and
// repeatedly insert a vertex into a random triangular face, joining it to the
// three corners. The result is always planar and fully triangulated.
static NetworKit::Graph makeApollonian(int extra, uint64_t seed) {
    NetworKit::Graph G(3 + extra, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    std::vector<std::array<NetworKit::node, 3>> faces = {{0, 1, 2}};
    std::mt19937_64 rng(seed);
    NetworKit::node next = 3;
    for (int i = 0; i < extra; i++) {
        std::uniform_int_distribution<size_t> pick(0, faces.size() - 1);
        size_t fi = pick(rng);
        auto f = faces[fi];
        NetworKit::node w = next++;
        G.addEdge(w, f[0]);
        G.addEdge(w, f[1]);
        G.addEdge(w, f[2]);
        faces[fi] = {f[0], f[1], w};
        faces.push_back({f[1], f[2], w});
        faces.push_back({f[2], f[0], w});
    }
    return G;
}

// Stacked triangular prisms ("tube"): k triangle rings, adjacent rings joined.
static NetworKit::Graph makePrismStack(int rings) {
    NetworKit::Graph G(3 * rings, false, false);
    auto id = [](int r, int i) { return 3 * r + i; };
    for (int r = 0; r < rings; r++) {
        G.addEdge(id(r, 0), id(r, 1));
        G.addEdge(id(r, 1), id(r, 2));
        G.addEdge(id(r, 2), id(r, 0));
        if (r + 1 < rings)
            for (int i = 0; i < 3; i++) G.addEdge(id(r, i), id(r + 1, i));
    }
    return G;
}

static void runAndVerify(const NetworKit::Graph& G) {
    Koala::PlanarSeparator algorithm(G);
    algorithm.run();
    verifySeparator(G, algorithm.getSeparator(), algorithm.getPartitionA(),
                    algorithm.getPartitionB());
}

TEST(PlanarSeparatorMore, TriangulatedGridsManySizes) {
    for (int s = 3; s <= 22; s++) {
        NetworKit::Graph G = buildTriangulatedGridMinusRandomNodes(s, s, 0.0, 42);
        runAndVerify(G);
    }
}

TEST(PlanarSeparatorMore, NonSquareTriangulatedGrids) {
    const std::vector<std::pair<int, int>> shapes = {
        {3, 10}, {5, 17}, {4, 25}, {2, 40}, {7, 13}, {10, 30}};
    for (auto [r, c] : shapes) {
        NetworKit::Graph G = buildTriangulatedGridMinusRandomNodes(r, c, 0.0, 7);
        runAndVerify(G);
    }
}

TEST(PlanarSeparatorMore, WheelsAndPrisms) {
    for (int rim : {6, 12, 20, 30, 50})
        runAndVerify(build_graph(rim + 1, makeWheelEdges(rim), false));
    for (int k : {3, 5, 10, 20, 40})
        runAndVerify(makePrismStack(k));
}

TEST(PlanarSeparatorMore, ApollonianRandomManySeeds) {
    for (uint64_t seed = 1; seed <= 40; seed++) {
        int extra = 10 + static_cast<int>(seed % 25) * 7; // n ~ 13..185
        NetworKit::Graph G = makeApollonian(extra, seed);
        runAndVerify(G);
    }
}

TEST(PlanarSeparatorMore, RandomPlanarWithDeletions) {
    for (uint64_t seed = 1; seed <= 40; seed++) {
        int rows = 8 + static_cast<int>(seed % 12);
        int cols = 8 + static_cast<int>((seed * 7) % 12);
        double del = 0.05 * static_cast<double>(seed % 6);
        NetworKit::Graph G =
            buildTriangulatedGridMinusRandomNodes(rows, cols, del, seed);
        if (G.numberOfNodes() < 3) continue;
        runAndVerify(G);
    }
}
