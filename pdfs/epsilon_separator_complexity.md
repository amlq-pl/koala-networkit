# EpsilonPlanarSeparator — why it is O(n²) now and how to get O(n log n)

## The algorithm being implemented

`EpsilonPlanarSeparator` implements **Algorithm 2 (Finding ϵ-planar separator)**:

```
1: C ← ∅
2: while G − C has a component with cost greater than ϵ do
3:     Let K be such a connected component
4:     Apply Theorem 2.1.1 to K producing a partition A₁, B₁, C₁
5:     C ← C ∪ C₁
6: end while
7: return C
```

Here Theorem 2.1.1 is the balanced planar separator (`PlanarSeparator`), which
splits `K` into `A`, `B`, `C` with `cost(A), cost(B) ≤ 2/3·cost(K)`.

Expected complexity: **O(n log n)** (or O(n log(1/ϵ))) — the separator work per
recursion level is linear in the pieces (which are disjoint, so ≤ n per level),
and there are O(log(1/ϵ)) levels because each split reduces cost by a 2/3 factor.

## Why the current implementation is O(n²)

The regression comes from **node-id gaps not being squashed**, not from any
single quadratic operation.

`NetworKit::GraphTools::subgraphFromNodes` **preserves the original node ids**.
So a sub-component with only `m` real vertices still reports a large
`upperNodeIdBound()` (close to `n_original`).

`PlanarSeparator` sizes almost all of its scratch storage by
`upperNodeIdBound()`, not by `numberOfNodes()`:

- `cpp/separator/PlanarSeparator.cpp:18` — BFS arrays `vis`, `lvl`, `parent`
  (`n = G.upperNodeIdBound()`)
- `findNumberOfVerticesAtLevel` iterates the full `lvl` vector
- `cpp/separator/PlanarSeparator.cpp:145` — `isConnectedToX`
- `cpp/separator/PlanarSeparator.cpp:176` — `childrenH`
- `cpp/separator/PlanarSeparator.cpp:182` — `costsH`

Consequence: each `PlanarSeparator::run()` on a component of true size `m` costs
**Θ(upperNodeIdBound) = Θ(n_original)** instead of Θ(m). Algorithm 2 processes
Θ(n) components down the recursion, so the total becomes **Θ(n²)**.

If the graphs handed to `PlanarSeparator` are compacted so that
`upperNodeIdBound() == numberOfNodes()`, each call costs Θ(m). Summed over the
O(log(1/ϵ)) balanced levels (pieces are disjoint, ≤ n per level) this is
**O(n log(1/ϵ))** — the intended O(n log n).

> Note: this is orthogonal to the *balance* concern (the Step-9 cycle-shrinking
> rebalancing in `PlanarSeparator` is not implemented). Both matter for the final
> bound, but **the id-gap issue is the one that turns the algorithm quadratic.**

## The fix — reuse an existing repo helper

The repo already uses exactly the right compaction pattern in
`include/matching/gaussian_matching/utils.hpp` (`reindexGraph`):

```cpp
auto indexes = NetworKit::GraphTools::getContinuousNodeIds(G); // oldId -> newId (0..m-1)
auto G1      = NetworKit::GraphTools::getCompactedGraph(G, indexes);
```

`getCompactedGraph` yields a graph with `upperNodeIdBound() == numberOfNodes()`.
Carry a `newId -> originalId` label vector so results stay in original ids.

### Proposed design (contained entirely in `EpsilonPlanarSeparator`)

Store `{compactedGraph, labels}` in the work queue instead of a bare `Graph`,
compact each piece as soon as it is created, and translate results back through
`labels`. No change to `PlanarSeparator` is required — its
`upperNodeIdBound`-based sizing becomes tight automatically once inputs are
compact.

```cpp
struct LabeledGraph {
  NetworKit::Graph graph;                  // compacted: upperNodeIdBound == numberOfNodes
  std::vector<NetworKit::node> toOriginal; // compacted id -> original id
};

static LabeledGraph compact(const NetworKit::Graph &sub,
                            const std::vector<NetworKit::node> &parentToOriginal) {
  auto idMap = NetworKit::GraphTools::getContinuousNodeIds(sub);   // sub-id -> new-id
  auto compacted = NetworKit::GraphTools::getCompactedGraph(sub, idMap);
  std::vector<NetworKit::node> toOriginal(compacted.upperNodeIdBound());
  for (const auto &[subId, newId] : idMap)
    toOriginal[newId] = parentToOriginal[subId];                  // compose the maps
  return {std::move(compacted), std::move(toOriginal)};
}
```

Then in `run()`:

- Build a local cost vector `localCost[i] = vertexCost[lg.toOriginal[i]]`
  (compact-indexed) to hand to `PlanarSeparator`.
- Push separator vertices back as `separator.push_back(lg.toOriginal[v])`.
- In `processComponents`, compact each connected component with
  `compact(sub, lg.toOriginal)` before enqueuing.

Every graph that reaches `PlanarSeparator` is then compact ⇒ per-call cost Θ(m)
⇒ overall O(n log(1/ϵ)).

## Open decision (affects output semantics)

The public `connectedComponents` output currently returns `Graph`s in original
ids. After compaction they would be in local ids. Cheapest correct options:

- Expose the parallel `toOriginal` labels alongside each component, or
- Return original-id node sets (`std::vector<std::vector<node>>`),

rather than rebuilding subgraphs on the original graph (which would reintroduce
an `upperNodeIdBound` cost per component and risk quadratic behavior again).

## References

- `cpp/separator/EpsilonPlanarSeparator.cpp` — Algorithm 2 driver
- `cpp/separator/PlanarSeparator.cpp` — Theorem 2.1.1 (balanced separator)
- `include/matching/gaussian_matching/utils.hpp` — existing `reindexGraph` pattern
- `NetworKit::GraphTools::getContinuousNodeIds` / `getCompactedGraph`
