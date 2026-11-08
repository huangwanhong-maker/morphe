# 0012. Centralities

- **Status:** Accepted, provisionally (2026-11-08).
- **Decided by:** Claude, for urban science's configuration of streets (the infrastructure's ADR 0084), since a centrality names no object in nature.

## 1. Decision

`combinatorics/centrality.hpp` holds closeness and betweenness from every source within a radius, by Brandes's algorithm (2001), on a graph whose arcs carry a cost and a length:
- **Least paths.** They are least in cost, then in length, to a relative tolerance. Costs may be zero (a street running straight on, under an angular cost) and lengths must be positive. Each arc then strictly raises (cost, length), so the predecessors of least paths are acyclic and their counts finite, even where zero-cost arcs form cycles. When cost is length, the least paths are the shortest paths.
- **The radius.** It is measured in length along the least paths, and the search does not pass beyond it.
- **Closeness.** For each source: the vertices reached, the sum of their least costs (total depth) and harmonic closeness.
- **Betweenness.** The share of least paths through each vertex, over ordered pairs or, halved, over unordered ones.
- **Splitting the sources.** Disjoint sets of sources add, so a caller may spread them over threads.

Like shortest paths (decision 0011), this is among the few floating-point structures of morphe, checked by `check()`.

## 2. Verification

`check_morphe_centrality`:
- **Closed forms.** A path's betweenness i(n − 1 − i), and a star's centre on all C(k, 2) pairs.
- **networkx 3.4.2.** Betweenness, total distance, harmonic closeness and vertices reached, on a random graph with real weights, a grid with unit weights and a geometric graph with integer weights (`tests/fixtures/centrality.txt`, written by `make_centrality_fixtures.py`). They agree to 2 × 10⁻¹⁶ relative.
- **Brute force.** Every simple path is enumerated on twenty random graphs whose costs include zeros: least paths by (cost, length), their betweenness and total cost, exactly. The same is done within a radius where cost is length.
- **Additivity and refusals.** Sources split in two give the whole when added, and negative costs and lengths that are not positive are refused.

## 3. Limits

Under a cost that is not length, the radius bounds the search's region, not a constrained least path (which would be NP-hard in general). A vertex whose least path leaves the radius is reached, if at all, by the least path within the region searched. This is the convention of Space Syntax's analyses at a metric radius.
