# combinatorics: Graphs and Combinatorial Structures

*Status: draft (2026-10-24): multigraphs, spanning forests, paths and distances; canonical forms and automorphism groups of labelled multigraphs.* Enumeration and generating functions, designs, graphs and hypergraphs as C-sets, flows, matchings, colourings, isomorphism with pruning by automorphisms, matroids, exact combinatorial optimisation and the enumeration of Feynman graphs with their automorphisms are planned ([decision 0003](../development/decisions/0003-logic-and-combinatorics.md), [decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md)).

## 1. Purpose

Combinatorial structures and their algorithms. The first is the graph, as the 1-skeleton of a complex: its components and spanning forests, through which the fundamental group of a complex is presented.

## 2. Concepts

- **A multigraph** has vertices and edges (tail, head), loops and parallel edges allowed: the plain C-set over the schema E ⇉ V ([idea of 2026-10-21](../development/ideas/2026-10-21-c-sets-as-common-data-model.md)), held as arrays until a second schema calls for the generic form.
- **A spanning forest**, grown breadth first in O(n + m), has one tree per component and n − c edges; the path between two vertices of a tree is read as signed edges, which is how a loop of a complex is written (Hatcher 2002, §1.A).
- **A canonical form** of a multigraph with labelled vertices and edges is an encoding equal for two graphs exactly when they are isomorphic. Colour refinement (Weisfeiler and Leman 1968) partitions the vertices by label and then, repeatedly, by the multiset of (edge label, neighbour's colour) until the partition is stable. Every isomorphism respects that partition, so only the orderings within its cells are tried: the least encoding over them is the form, and the number that attain it is the order of the automorphism group (individualisation by refinement in its simplest form; McKay and Piperno 2014, §2). The search is exact and exponential only in the cells that refinement cannot split. Beyond a limit, the stable partition's encoding is returned, marked inexact: isomorphic graphs share it, and non-isomorphic ones rarely do.

## 3. Content

| Header | Content |
|---|---|
| `graph.hpp` | `Multigraph{vertices, edges}` with `incidence`; `Step{edge, sign}`; `spanning_forest(G, first)` → `SpanningForest` with `path(u, v, G)`, `components`; `distances(G, source)` |
| `canonical.hpp` | `LabelledGraph{vertex, edges}`; `canonical_form(G, limit)` → `CanonicalForm{code, automorphisms, exact, order}` |

## 4. Checks

`tests/check_combinatorics.cpp`: on 200 random multigraphs, n − c tree edges and forest paths that are walks between their ends; distances around a hexagon.

`tests/check_canonical.cpp` (15 checks): automorphism groups known by hand: the cycles C_3 to C_8 (2n), K_4 (24), the star K_{1,3} (6), the Petersen graph (120, where refinement splits nothing); labels and a doubled edge breaking symmetry; forty random relabellings of the Petersen graph sharing its code; the hexagon and two triangles, alike to colour refinement, told apart.

## 5. Dependencies

None.

## 6. References

- Hatcher, A. (2002), *Algebraic Topology*, Cambridge University Press, §1.A.
- Massey, W. S. (1991), *A Basic Course in Algebraic Topology*, Springer, ch. VI.
- Patterson, E., Lynch, O. and Fairbanks, J. (2022), Categorical data structures for technical computing, *Compositionality* 4, 5.
- Weisfeiler, B. and Leman, A. (1968), The reduction of a graph to canonical form and the algebra which appears therein, *Nauchno-Technicheskaya Informatsia* 2(9), 12–16.
- McKay, B. D. and Piperno, A. (2014), Practical graph isomorphism, II, *J. Symbolic Comput.* 60, 94–112.

## 7. Open Questions

1. **C-sets** as the common representation, once a second schema (semi-simplicial sets, Petri nets) is written.
2. **Graph isomorphism at scale.** The canonical form tries every ordering within the stable cells. Individualising one vertex at a time and refining again, with the automorphisms found pruning the search tree (McKay and Piperno 2014), would make highly symmetric graphs cheap; the symmetry factors of Feynman graphs and the types of a history's bubbles would use it.
