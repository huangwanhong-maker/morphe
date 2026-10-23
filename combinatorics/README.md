# combinatorics: Graphs and Combinatorial Structures

*Status: draft (2026-10-23): multigraphs, spanning forests, paths and distances.* Enumeration and generating functions, designs, graphs and hypergraphs as C-sets, flows, matchings, colourings, isomorphism by refinement, matroids, exact combinatorial optimisation and the enumeration of Feynman graphs with their automorphisms are planned ([decision 0003](../development/decisions/0003-logic-and-combinatorics.md), [decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md)).

## 1. Purpose

Combinatorial structures and their algorithms. The first is the graph, as the 1-skeleton of a complex: its components and spanning forests, through which the fundamental group of a complex is presented.

## 2. Concepts

- **A multigraph** has vertices and edges (tail, head), loops and parallel edges allowed: the plain C-set over the schema E ⇉ V ([idea of 2026-10-21](../development/ideas/2026-10-21-c-sets-as-common-data-model.md)), held as arrays until a second schema calls for the generic form.
- **A spanning forest**, grown breadth first in O(n + m), has one tree per component and n − c edges; the path between two vertices of a tree is read as signed edges, which is how a loop of a complex is written (Hatcher 2002, §1.A).

## 3. Content

| Header | Content |
|---|---|
| `graph.hpp` | `Multigraph{vertices, edges}` with `incidence`; `Step{edge, sign}`; `spanning_forest(G, first)` → `SpanningForest` with `path(u, v, G)`, `components`; `distances(G, source)` |

## 4. Checks

`tests/check_combinatorics.cpp`: on 200 random multigraphs, n − c tree edges and forest paths that are walks between their ends; distances around a hexagon.

## 5. Dependencies

None.

## 6. References

- Hatcher, A. (2002), *Algebraic Topology*, Cambridge University Press, §1.A.
- Massey, W. S. (1991), *A Basic Course in Algebraic Topology*, Springer, ch. VI.
- Patterson, E., Lynch, O. and Fairbanks, J. (2022), Categorical data structures for technical computing, *Compositionality* 4, 5.

## 7. Open Questions

1. **C-sets** as the common representation, once a second schema (semi-simplicial sets, Petri nets) is written.
2. **Graph isomorphism and automorphisms** by partition refinement, for the symmetry factors of Feynman graphs.
