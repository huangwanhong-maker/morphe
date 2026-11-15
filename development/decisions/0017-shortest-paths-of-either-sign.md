# 0017. Shortest Paths with Weights of Either Sign

- **Status:** Accepted, provisionally (2026-11-15).
- **Date:** 2026-11-15
- **Decided by:** Claude, for the evidence foundation tekmerion's timelines (the infrastructure's ADR 0115 and tekmerion's decision 0005), continuing the user's request for evidence and its reconstruction.

## 1. Context

Decision 0011 gave shortest paths with weights that are not negative (Dijkstra). A set of bounds on differences of times is consistent exactly when its distance graph has no negative cycle, and its tightest bounds are its shortest paths. The arcs of a lower bound are negative. The general method is a graph's, and so morphe's.

## 2. Decision

`combinatorics/shortest_paths.hpp` gains:

- `bellman_ford`: from given sources with offsets, or from every vertex at 0 when none is given, which finds any negative cycle and gives potentials of the whole graph;
- `NegativeCycle`: the cycle as a list of vertices whose arcs are the graph's, with its weight;
- `floyd_warshall`: every pair, with the next vertex on each shortest path, the negative cycle named by Bellman–Ford when one exists.

Weights are doubles, as Dijkstra's are: they are data of the graph, not algebraic identities (invariant 3).

## 3. Consequences

`check_morphe_shortest_paths` gains four checks. Bellman–Ford equals Dijkstra on non-negative weights. On 400 random graphs with negative arcs, a negative cycle is found exactly when a search of every simple cycle finds one, the cycle returned is made of the graph's arcs with the weight stated, and otherwise Floyd–Warshall equals Bellman–Ford from every vertex, each path as long as its distance.
