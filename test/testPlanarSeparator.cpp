#include <gtest/gtest.h>
#include <separator/PlanarSeparator.hpp>
#include <unordered_set>
#include <cmath>
#include "helpers.hpp"

// ---- Helper: verify separator properties ----
// Call after algorithm.run() to check:
//   1. S ∪ A ∪ B = V  (all vertices accounted for)
//   2. A ∩ B = ∅, A ∩ S = ∅, B ∩ S = ∅
//   3. No edge between A and B (S is a valid separator)
//   4. |A| <= 2n/3 and |B| <= 2n/3  (balance)
//   5. |S| <= c * sqrt(n)  (size bound)
//
// void verifySeparator(const NetworKit::Graph& G,
//                      const std::vector<NetworKit::node>& separator,
//                      const std::vector<NetworKit::node>& partA,
//                      const std::vector<NetworKit::node>& partB) {
//     NetworKit::count n = G.numberOfNodes();
//     // Check union covers all vertices
//     // Check pairwise disjoint
//     // Check no A-B edge
//     // Check |A|, |B| <= 2n/3
//     // Check |S| <= 2*sqrt(2)*sqrt(n)  (or some reasonable constant)
// }

// ---- Test cases ----
// Suggested test graphs:
//
// 1. Trivial: single vertex, single edge, triangle
// 2. Path graph (P_n): separator should be a single middle vertex
// 3. Cycle graph (C_n): separator = 2 vertices
// 4. Grid graph (k x k): separator ~ sqrt(n)
// 5. Wheel graph: center vertex is natural separator
// 6. Larger planar graphs from input/*.g6 files

struct SeparatorParameters {
    std::string name;
    int N;
    std::list<std::tuple<int, int, int>> E;
};

class PlanarSeparatorTest : public testing::TestWithParam<SeparatorParameters> {};

TEST_P(PlanarSeparatorTest, TestSeparatorProperties) {
    SeparatorParameters const& parameters = GetParam();
    std::cout << "\nTest name: " << parameters.name << std::endl;

    NetworKit::Graph G = build_graph(parameters.N, parameters.E, false);

    Koala::PlanarSeparator algorithm(G);
    algorithm.run();

    // auto& sep = algorithm.getSeparator();
    // auto& partA = algorithm.getPartitionA();
    // auto& partB = algorithm.getPartitionB();
    // verifySeparator(G, sep, partA, partB);
}

// ---- Test data ----
// Example: path graph P5:  0-1-2-3-4
INSTANTIATE_TEST_SUITE_P(
    PlanarSeparatorTests,
    PlanarSeparatorTest,
    testing::Values(
        SeparatorParameters{"path_5", 5, {{0,1,1},{1,2,1},{2,3,1},{3,4,1}}},
        SeparatorParameters{"cycle_6", 6, {{0,1,1},{1,2,1},{2,3,1},{3,4,1},{4,5,1},{5,0,1}}},
        SeparatorParameters{"grid_3x3", 9, {
            {0,1,1},{1,2,1},{3,4,1},{4,5,1},{6,7,1},{7,8,1},
            {0,3,1},{1,4,1},{2,5,1},{3,6,1},{4,7,1},{5,8,1}
        }},
        SeparatorParameters{"triangle", 3, {{0,1,1},{1,2,1},{2,0,1}}}
    )
);
