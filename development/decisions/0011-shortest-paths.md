# 0011. Shortest Paths

- **Status:** Accepted, provisionally (2026-11-07).
- **Decided by:** Claude, for the urban-science domain's `access` (the infrastructure's ADR 0078), since a shortest path names no object in nature.

## 1. Decision

`combinatorics/shortest_paths.hpp` holds Dijkstra's algorithm with a binary heap on a weighted adjacency: from one source or several with offsets (the distance to the nearest), with a limit, reporting each vertex's distance, parent and origin. Weights are real numbers, so that this is among the few floating-point structures of morphe, and is checked by `check()` rather than exactly.

## 2. Verification

`check_morphe_shortest_paths`: Manhattan distances on a grid, nearest sources on a line, a limit, equality with Bellman–Ford on random graphs, and refusals of negative weights.
