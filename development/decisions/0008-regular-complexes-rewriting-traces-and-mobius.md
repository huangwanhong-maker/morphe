# 0008. Regular Complexes and Their Rewriting, Traces, and Möbius Inversion

- **Status:** Accepted, provisionally. Regularity is decided exactly up to dimension 3 and by homology mod 2 above it. Rewriting reads a complex as its incidence graph, with injective matches. The count of linearisations is exact up to 40 events.
- **Date:** 2026-10-31
- **Decided by:**
  - **The request.** The user asked for the formalism of their manuscript on urban landscapes to be implemented, saying that "some module may resuable, like CW-complex".
  - **The rule.** General methods are written in the foundations.
  - **What Claude did.** It wrote the headers, their checks and this record.
  - **The composition.** The infrastructure's [ADR 0062](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0062-landscapes-as-cellular-configurations.md).

## 1. Context

A place whose structure changes is a cell complex rewritten by events: a street subdivided where a passage meets it, a block split, two rooms merged. morphe had CW complexes for computing invariants (decision 0004's topology) and graph rewriting by double pushouts (decision 0007). It had no complex that could be rewritten and stay a complex. Nor did it have the order in which a record of such events must have happened, or Möbius inversion over a poset, which says what several interventions do together beyond each alone.

## 2. Decision

### 2.1 Regular Complexes (`topology/regular_complex.hpp`)

A regular complex is presented by its face poset: cells with dimensions, facets with incidence signs, types and persistent identifiers (Björner 1984, *Eur. J. Combin.* 5, 7: a regular CW complex is determined by its face poset). `admissible` decides regularity exactly up to dimension 3:
- a 1-cell has two distinct ends;
- the facets of a 2-cell form one cycle;
- the facets of a 3-cell form a closed connected surface, every vertex link one cycle, of Euler characteristic 2.

Above dimension 3 the test is the homology mod 2 of each boundary, and is marked inexact. `orient` propagates orientations across the poset's diamonds, `chain_complex` gives the cellular chains, and `cubical_grid` builds boxes of cubes, a building's rooms among them.

### 2.2 Rewriting a Complex (`topology/complex_rewriting.hpp`)

A complex is rewritten as its incidence graph by decision 0007's double pushout, then read back and tested for regularity. A rewrite that would leave a non-regular complex is refused with the reason.
- **Rules built at a site:** `attach_cell`, `attach_edge`, `remove_cell`, `subdivide_edge`, `merge_at_vertex`, `partition_cell`, `split_face`, `merge_cells`, `relabel`. Each records the cells created, deleted and relabelled, and its footprint (what it read, what it wrote).
- **`isomorphic`**, of complexes, oriented or not, by a search for a match. The canonical forms of `combinatorics/canonical.hpp` are not fine enough for unsigned incidence graphs.

### 2.3 Traces (`formal/traces.hpp`)

Footprints are sets of what an event reads and writes; two events are dependent when one writes what the other reads or writes (Mazurkiewicz 1977; Diekert and Rozenberg 1995). `causal_order` gives a record's order (dependence, closed transitively), its covers and its Foata steps (Cartier and Foata 1969). `linearisations` counts the orders of the record that are the same history, exactly with morphe's integers.

### 2.4 Möbius Inversion (`combinatorics/mobius.hpp`)

The Möbius function of a finite poset, the inversion of a function on it, and the subset lattice's zeta and Möbius transforms (Rota 1964, *Z. Wahrsch.* 2, 340). `anchored_interactions` gives, for an outcome of any subset of n interventions, what each subset adds beyond its proper subsets.

## 3. Alternatives Considered

- **Combinatorial maps or generalised maps** as the representation, which make 2- and 3-manifolds' rewriting local. They do not hold the non-manifold complexes a district's network gives (a junction of five streets, a dead end), and the face poset does.
- **Testing regularity by homology alone.** It cannot tell a sphere from a homology sphere in dimension 3, and the exact test there is cheap.
- **Isomorphism by canonical forms.** They are too coarse on unsigned incidence graphs; a search is exact, and complexes of a few hundred cells are small enough.

## 4. Consequences

- The topology module now depends on `formal` (traces and footprints), and CMake builds `formal` before `topology`.
- Checks: `check_regular_complex` (38), `check_complex_rewriting` (31), `check_traces` (11), `check_mobius` (11).
- The urban-science domain builds its districts and buildings on these headers (the infrastructure's ADR 0062).
