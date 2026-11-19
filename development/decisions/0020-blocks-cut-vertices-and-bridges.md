# 0020. Blocks, Cut Vertices and Bridges

- **Status:** Accepted, provisionally (2026-11-19).
- **Date:** 2026-11-19
- **Decided by:** Claude, for the infrastructure's diplomacy domain (its ADR 0132). An undertaking among many parties must say whether it can continue without a particular coordinator: whether one participant, or one tie, holds the relations of authority together.

## 1. Context

The infrastructure's diplomacy domain reads an undertaking as a typed multigraph of participants. One of its diagnostics asks whether the undertaking can continue without a particular participant: the "succession" failure of a coordinating function, in the user's monograph. That is the question whether the participant is a cut vertex of the graph of the relations that carry authority, and whether a tie is a bridge. It names no participant or relation, so it belongs here. The combinatorics module had spanning forests and components, and no decomposition into biconnected parts.

## 2. Decision

`combinatorics/blocks.hpp`: the blocks, cut vertices and bridges of a multigraph read as undirected. They are found in one iterative depth-first search by Hopcroft and Tarjan's lowpoints, in O(n + m):

- a vertex is a cut vertex when removing it leaves more components among the others;
- an edge is a bridge when removing it leaves more components;
- every edge that is not a loop lies in exactly one block, and two blocks share at most one vertex, a cut vertex.

Loops lie in no block. Parallel edges are told apart by their identities, so that two edges joining the same pair are never a bridge.

The checks hold it to small graphs whose answers are known, to removal on 400 random multigraphs, and to a path of 200,000 vertices (no recursion to exhaust).

## 3. Consequences

- The block–cut tree, two-edge-connected components and the directed analogues (dominators, strong articulation points) are left until a module needs them.
- 8 checks in `check_morphe_blocks`.
