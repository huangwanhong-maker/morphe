# 0019. Flows in Networks

- **Status:** Accepted, provisionally (2026-11-18).
- **Date:** 2026-11-18
- **Decided by:** Claude, at the user's request for the risks of cities (the infrastructure's ADR 0129): which roads matter most to reaching a hospital, how many separate ways lead out of a district, and how quickly a place can be cleared towards its shelters.

## 1. Context

The infrastructure's new domains of urban risk and disaster risk need three things of a graph that name no road, shelter or hazard: how much can pass from a set of places to another, and which arcs bound it (a cut); how to send people at least cost to destinations of limited capacity (an assignment to shelters); and how much can arrive by each instant when moving takes time and arcs pass so much a step (a flow over time, whose value bounds any evacuation from above). Edge-disjoint paths, the redundancy of a place, are a flow with unit capacities. The README planned "paths, cuts, flows, matchings" in `combinatorics`; shortest paths are there (decisions 0011, 0017).

## 2. Decision

`combinatorics/flows.hpp`:

- `max_flow`: Dinic's algorithm, with the least cut (the vertices still reached from the source in the residual network) as its certificate;
- `min_cost_flow`: successive shortest paths under potentials, with a price below which a unit is worth shipping, which gives the cheapest flow of any value and Ford and Fulkerson's temporally repeated flows;
- `decompose`: a flow as simple paths with their amounts;
- `flow_over_time`: the time-expanded network of whole steps, with supplies released at their steps, waiting with a storage, sinks of limited capacity over the horizon, and the most that can reach a sink by every step, found on one network by opening the sinks step by step and augmenting what was found.

Capacities, costs and amounts are doubles, data of the network as weights are (decision 0017, invariant 3), compared with a tolerance relative to the largest finite capacity. Each value of `safe` is the greatest for its own step; one flow need not attain them all at once, and the header says so. How a continuous problem is bounded from above by the discrete one (transit times rounded down, the horizon up) is stated in the header and is the caller's to do.

## 3. Consequences

`check_morphe_flows` (13 checks) holds the greatest flow to every cut listed on random networks, the cheapest flow to every permutation of random assignments with its potentials as a certificate, and the time-expanded value at every horizon to the temporally repeated flow computed from the cheapest flow, the theorem of Ford and Fulkerson, on random networks. A flow over time of n vertices and H steps has (H + 1) n vertices; at a city's scale a step of ten seconds or more keeps it within memory. Matchings and flows with gains are not written.
