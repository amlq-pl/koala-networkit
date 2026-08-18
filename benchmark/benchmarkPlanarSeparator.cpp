// Time-complexity benchmark for Koala::PlanarSeparator.
//
// It builds planar triangulated grids of geometrically increasing size,
// times PlanarSeparator::run() on each, and prints a table whose last column
// is the empirical exponent p in an assumed Theta(n^p) running time:
//
//     p = log(t / t_prev) / log(n / n_prev)
//
// For the current O(n^2) implementation p converges towards 2; once Step 9 is
// made linear the same benchmark (unchanged) should show p trending to 1.
//
// Output rows are whitespace-separated so they can be plotted directly
// (gnuplot/awk/python skip the leading '#' comment lines).
//
// Usage:
//     benchmark_planar_separator [startN] [maxN] [reps] [timeBudget]
//
//   startN      first target vertex count           (default 1000)
//   maxN        hard cap on target vertex count      (default 5000000)
//   reps        timed repetitions per size (median)  (default 3)
//   timeBudget  stop once a single run exceeds this  (default 10 seconds)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <separator/PlanarSeparator.hpp>

namespace {

// Build a triangulated (planar) grid with optional random vertex deletions.
// Vertex (r, c) has id r * cols + c before relabelling; survivors are
// relabelled to a contiguous 0..N-1 range. Replicated from
// test/testPlanarSeparator.cpp so this benchmark stays self-contained.
NetworKit::Graph buildTriangulatedGrid(int rows, int cols, double deletionRate,
                                       uint64_t seed) {
  const int total = rows * cols;
  auto idx = [cols](int r, int c) { return r * cols + c; };

  std::vector<std::pair<int, int>> edges;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      if (c + 1 < cols)
        edges.emplace_back(idx(r, c), idx(r, c + 1));
      if (r + 1 < rows)
        edges.emplace_back(idx(r, c), idx(r + 1, c));
      if (r + 1 < rows && c + 1 < cols) {
        // Diagonal: alternate orientation for variety.
        if ((r + c) % 2 == 0)
          edges.emplace_back(idx(r, c), idx(r + 1, c + 1));
        else
          edges.emplace_back(idx(r, c + 1), idx(r + 1, c));
      }
    }
  }

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<bool> deleted(total, false);
  for (int v = 0; v < total; v++) {
    if (dist(rng) < deletionRate)
      deleted[v] = true;
  }

  std::vector<int> newId(total, -1);
  int next = 0;
  for (int v = 0; v < total; v++) {
    if (!deleted[v])
      newId[v] = next++;
  }
  const int N = next;

  NetworKit::Graph G(N, false, false);
  for (auto [u, v] : edges) {
    if (!deleted[u] && !deleted[v])
      G.addEdge(newId[u], newId[v]);
  }
  return G;
}

struct RunResult {
  double timeSeconds;
  std::size_t separatorSize;
  std::size_t sizeA;
  std::size_t sizeB;
};

// Runs PlanarSeparator up to `reps` times (stopping early once a single run
// exceeds `budget`) and returns the median wall-clock time plus the partition
// sizes of the last run.
RunResult timePlanarSeparator(const NetworKit::Graph &G, int reps,
                              double budget) {
  std::vector<double> times;
  std::size_t sepSize = 0, sizeA = 0, sizeB = 0;
  for (int i = 0; i < reps; i++) {
    Koala::PlanarSeparator algorithm(G);
    auto start = std::chrono::high_resolution_clock::now();
    algorithm.run();
    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration<double>(end - start).count();
    times.push_back(t);
    sepSize = algorithm.getSeparator().size();
    sizeA = algorithm.getPartitionA().size();
    sizeB = algorithm.getPartitionB().size();
    if (t > budget)
      break; // no point repeating an already-slow run
  }
  std::sort(times.begin(), times.end());
  return {times[times.size() / 2], sepSize, sizeA, sizeB};
}

long long parseArg(int argc, char **argv, int index, long long fallback) {
  if (index >= argc)
    return fallback;
  return std::strtoll(argv[index], nullptr, 10);
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 1 &&
      (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
    std::cerr
        << "Usage: " << argv[0]
        << " [startN] [maxN] [reps] [timeBudget]\n"
           "  Sweeps triangulated planar grids of growing size and prints\n"
           "  n, m, time_s, |S|, |A|, |B| and the empirical exponent p.\n";
    return 0;
  }

  const long long startN = parseArg(argc, argv, 1, 1000);
  const long long maxN = parseArg(argc, argv, 2, 5000000);
  const int reps = static_cast<int>(parseArg(argc, argv, 3, 3));
  const double timeBudget =
      argc > 4 ? std::strtod(argv[4], nullptr) : 10.0;
  const double growth = 2.0; // target vertex count roughly doubles each step
  const uint64_t seed = 12345;

  std::cout << "# PlanarSeparator time-complexity benchmark\n";
  std::cout << "# startN=" << startN << " maxN=" << maxN << " reps=" << reps
            << " timeBudget=" << timeBudget << "s\n";
  std::cout << "#" << std::setw(11) << "n" << std::setw(13) << "m"
            << std::setw(14) << "time_s" << std::setw(9) << "|S|"
            << std::setw(10) << "|A|" << std::setw(10) << "|B|"
            << std::setw(9) << "exp_p" << "\n";

  double prevTime = 0.0;
  long long prevN = 0;

  for (long long targetN = startN; targetN <= maxN;
       targetN = static_cast<long long>(std::ceil(targetN * growth))) {
    // Square-ish triangulated grid with about `targetN` vertices. The side is
    // clamped so side*side stays within int range.
    int side = static_cast<int>(std::llround(std::sqrt(
        static_cast<double>(targetN))));
    if (side < 2)
      side = 2;
    if (side > 40000)
      side = 40000;

    NetworKit::Graph G = buildTriangulatedGrid(side, side, 0.0, seed);
    const long long n = static_cast<long long>(G.numberOfNodes());
    const long long m = static_cast<long long>(G.numberOfEdges());

    RunResult r = timePlanarSeparator(G, reps, timeBudget);

    double expP = 0.0;
    if (prevN > 0 && prevTime > 0.0 && r.timeSeconds > 0.0)
      expP = std::log(r.timeSeconds / prevTime) /
             std::log(static_cast<double>(n) / static_cast<double>(prevN));

    std::cout << std::setw(12) << n << std::setw(13) << m << std::setw(14)
              << std::fixed << std::setprecision(6) << r.timeSeconds
              << std::setw(9) << r.separatorSize << std::setw(10) << r.sizeA
              << std::setw(10) << r.sizeB << std::setw(9)
              << std::setprecision(3) << expP << "\n"
              << std::flush;

    prevTime = r.timeSeconds;
    prevN = n;

    if (r.timeSeconds > timeBudget) {
      std::cout << "# stopping: run time " << r.timeSeconds
                << "s exceeded budget " << timeBudget << "s\n";
      break;
    }
  }

  return 0;
}
