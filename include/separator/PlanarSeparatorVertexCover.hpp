#include "networkit/base/Algorithm.hpp"
#include "networkit/graph/Graph.hpp"

namespace Koala {
class PlanarSeparatorVertexCover : NetworKit::Algorithm {
public:
  explicit PlanarSeparatorVertexCover(NetworKit::Graph &G);
  void run() override;

private:
  NetworKit::Graph graph;
  NetworKit::Graph residual;
  void prep();
};
} // namespace Koala
