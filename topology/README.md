# topology: Complexes and Their Invariants

*Status: draft (2026-10-23): CW and simplicial complexes, a catalogue of constructions, homology, cohomology and fundamental groups, surfaces, manifolds, collapses, dual 2-complexes and coloured graphs; persistent homology over ℤ/2, the cycles over ℤ/2, Atkin's Q-analysis and the worldline complex of a record of encounters (2026-10-24).* Cubical complexes, simplicial sets, persistence over other fields, cup products, discrete Morse theory, the recognition of 3-manifolds and knots are planned ([decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md), §2.2).

## 1. Purpose

Spaces built from cells, and what can be computed of them exactly. The 2-complex is the first object: it carries the fundamental group and the first homology groups; its 2-cells are where the curvature of a lattice gauge field is read; and the Feynman graphs of group field theory are 2-complexes, the dual 2-skeletons of the pseudomanifolds they encode. physis builds its gauge fields, spin foams and discrete spacetimes on these complexes.

## 2. Concepts

- **A CW complex, combinatorially.** Each 1-cell has a tail and a head; each 2-cell a base vertex and an attaching word, a closed edge path; each higher cell its cellular boundary as integer incidences (Hatcher 2002, §2.2). Up to homotopy in dimension 2, and for homology in every dimension, this is all that matters. The cellular chain complex is read from it, and `validate` checks that the words close and that ∂∂ = 0.
- **Simplicial complexes**, given by facets and holding every face, oriented by vertex order; links, stars, skeleta, the barycentric subdivision, joins, cones and suspensions (Spanier 1966; Rourke and Sanderson 1972). Each is also a CW complex, a triangle attached by [v₀v₁][v₁v₂][v₀v₂]⁻¹.
- **Constructions.** Polygons glued by words ("abAB" the torus, "aa" the projective plane, "abAb" the Klein bottle, "aaA" the dunce hat), with corners identified by union-find; spheres, real and complex projective spaces and lens spaces by their minimal cell structures; presentation complexes; products of CW complexes by the Leibniz rule, the 2-cells keeping attaching words, and with them lattice tori of any dimension; the standard small triangulations.
- **The fundamental group** on the edges outside a spanning tree of the 1-skeleton, with a relator from each 2-cell and a loop at the base for each generator (Hatcher, Proposition 1.26); its simplification and abelianisation by morphe `algebra`.
- **Surfaces** are recognised exactly (every edge on two sides or one, every vertex link a circle or an arc) and classified by orientability, Euler characteristic and boundary circles (Massey 1991, ch. 1), with a consistent orientation of the 2-cells when one exists.
- **Manifolds** among simplicial complexes: pseudomanifolds, and combinatorial manifolds by their vertex links, exactly in dimension at most 3, by homology above.
- **Collapses** of 2-complexes: free edges of 2-cells and free vertices of edges (Whitehead 1939). In dimension 2 the greedy order decides collapsibility; the dunce hat is contractible and not collapsible (Zeeman 1964).
- **The dual 2-complex** of a triangulated pseudomanifold: a vertex per top simplex, an edge per interior facet, a 2-cell per cycle around a codimension-2 simplex: the 2-complex of spin foams and lattice gauge theory on a triangulation (Reisenberger and Rovelli 1997; Oriti 2001).
- **Cycles over ℤ/2.** The k-cycles of a complex modulo 2, Z_k(X; ℤ/2) = ker ∂_k, with a basis by elimination over GF(2) and the support of the space: the cells that lie in some cycle. In dimension 2 they are the spin foams of a ℤ₂ gauge theory on X, which exist on non-orientable surfaces too (Reisenberger and Rovelli 1997; physis `gauge`).
- **The worldline complex.** A record of encounters (events in time, the participants present at each, pairs with the events at which they meet) makes a 2-complex with no further choice. Its vertices are the events, its edges each participant's segments between consecutive events, and its 2-cells the lenses, each pair between consecutive meetings, attached along a[p, q] · b[p, q]⁻¹. Its refinement opens each event into its presences, with a rung for each pair meeting there and a triangle for each three pairs that join three participants. It is the construction of the infrastructure's history complex and of its economics domain's spacetime of transactions.
- **Q-analysis.** A relation between two finite sets, read row by row, is a family of simplices (Dowker 1952). Two are q-near when they share at least q + 1 vertices, and q-connected through a chain of q-near ones; the numbers Q_q of q-connected components among the simplices of dimension at least q form the structure vector of the family (Atkin 1974, ch. 3).
- **Coloured graphs**, (D+1)-regular and properly edge-coloured, encode D-dimensional pseudomanifolds (Pezzana 1974; Ferri, Gagliardi and Grasselli 1986): k-cells are the components of the subgraphs on D − k colours. They are the Feynman graphs of coloured tensor models and group field theories (Gurau 2011); their jackets, degree and melonic family govern the large-N expansion (Bonzom, Gurau, Riello and Rivasseau 2011).

## 3. Content

| Header | Content |
|---|---|
| `cw_complex.hpp` | `Face{base, boundary}`, `Path` of `Step`; `CWComplex{name, vertices, edges, faces, higher, edge_names}` with `dimension`, `cells`, `counts`, `euler_characteristic`, `walk`, `one_skeleton`, `chain_complex`, `validate` |
| `simplicial.hpp` | `SimplicialComplex::from_facets`, `find`, `facets`, `is_pure`, `cofaces`, `chain_complex`, `cw`, `link`, `star`, `skeleton`; `renumbered`, `barycentric_subdivision`, `join`, `cone`, `suspension` |
| `constructions.hpp` | `parse_polygons`, `polygon_complex`; `orientable_surface`, `nonorientable_surface`, `klein_bottle`, `dunce_hat`, `sphere`, `real_projective_space`, `complex_projective_space`, `lens_space`, `presentation_complex`, `circle`, `interval`, `product`, `torus_lattice`; `simplex`, `simplex_boundary`, `torus_7`, `projective_plane_6`, `octahedron` |
| `invariants.hpp` | `homology`, `betti`; `fundamental_group(X, base)` → `FundamentalGroup{presentation, generator_edges, loops, simplified, simplified_loops, abelianisation}` |
| `cycles.hpp` | `cycles_mod2(X, k)` → `Cycles{dimension, basis, support}` |
| `q_analysis.hpp` | `q_analysis(family)` → the structure vector Q |
| `worldlines.hpp` | `Encounters{participants, present, pairs}`; `worldline_complex(R)` → `WorldlineComplex{cw, worldline, segments, segments_of, lenses, lenses_of, path}`; `refined_worldlines(R, W)` → `RefinedWorldlines{cw, presences, presence_of, temporal, rung_edge, rung_of_edge, triangles}` |
| `persistence.hpp` | `Filtration{vertex, edge, face, higher}` (a stage for each cell, faces no later than cofaces); `persistent_homology(X, F, zero_length)` → `PersistenceDiagram{intervals}`, each `PersistenceInterval{dimension, birth, creator, death, killer}`, and `betti(dim, stage)`: the standard reduction of the boundary matrix over ℤ/2 in the order of the filtration (Edelsbrunner, Letscher and Zomorodian 2002; Zomorodian and Carlsson 2005; Edelsbrunner and Harer 2010, §VII.1) |
| `surfaces.hpp` | `classify_surface(X)` → `SurfaceType{surface, closed, orientable, euler, genus, crosscaps, boundary_components, name, normal_form, reason, orientation}` |
| `manifolds.hpp` | `manifold_type(K)` → `ManifoldType{dimension, pure, pseudomanifold, manifold, exact, closed, reason}` |
| `collapse.hpp` | `collapse(X)` → `Collapse{attempted, collapsible, steps, moves, remaining, …}` |
| `dual.hpp` | `dual_two_complex(K)` → `DualTwoComplex{complex, vertex_simplex, edge_simplex, face_simplex}` |
| `coloured_graph.hpp` | `ColouredGraph{colours, p, sigma}` with `components`, `faces`, `face_count`, `bubbles`, `jacket_genera`, `degree`, `two_complex`, `dual_complex`; `elementary_melon`, `insert_melon`, `random_coloured_graph` |
| `bubbles.hpp` | bubbles, the interactions of invariant theories as bipartite r-coloured graphs: (φ̄φ)^k, (φ̄·φ)^k, tr(MM†)^k, the quadratic bubble, quartic melons, the wheel K₃,₃, the necklace; their automorphism groups; the gluing of bubbles by a propagator matching into (r+1)-coloured Feynman graphs (physis decision 0006) |

## 4. Checks

`tests/check_topology.cpp`: the homology and classification of eleven complexes given by words (sphere, torus, genus 2, projective plane, Klein bottle, three crosscaps, Möbius band, annulus, disc, two tori on a common circle, two disjoint tori); π₁ of the torus with its loops; the H₂ generator of a lattice torus is its fundamental class Σ o_f f for the orientation found; ℝP⁵, ℂP³, S⁴ and L(5, q) with its π₁; ℝP² × S¹ against Künneth; lattice 3- and 4-tori (Betti numbers 1, 4, 6, 4, 1) and a 32 × 32 torus with its generators; the presentation complex of the (2, 3, 5) triangle group, whose π₁ has 121 homomorphisms into A₅; the 7-vertex torus, the 6-vertex projective plane, a barycentric subdivision, ∂Δ⁴ and a suspension recognised as 3-manifolds, a cone as a ball, two spheres on a point as a pseudomanifold only; collapses (a grid square and a disc collapse, the dunce hat and the torus do not); the dual 2-complexes of ∂Δ⁴ and of the 7-vertex torus; coloured graphs: the melon and melonic insertions encode S³, Gurau's face formula on random graphs, the surfaces of 3-coloured graphs have genus equal to their degree, and rooted melonic 4-coloured graphs are counted by the Fuss–Catalan numbers 1, 4, 22, 140.

`tests/check_worldlines.cpp` (12 checks): three participants meeting pairwise at two events make the theta-sphere (Betti 1 0 1) and a refined prism, also a sphere; a circle of encounters met once each is a graph with one loop; paths run along worldlines in time and not against it; on sixty random records both complexes are valid, with the cells the construction makes and words of the right lengths; faulty records are refused.

`tests/check_cycles.cpp` (16 checks): on the torus, the surface of genus 2, the projective plane, the Klein bottle and a 4 × 3 lattice torus, one 2-cycle over ℤ/2, the whole surface, whose basis is closed modulo 2 and whose number is the Betti number over ℤ/2 of persistence; a disc beside a torus outside every 2-cycle; the 1-cycles of a graph, E − V + c of them; the Q-vectors of two families worked by hand, unchanged by relabelling, with Q_0 the number of components of the complex they generate.

## 5. Dependencies

morphe `algebra` (and through it `arithmetic`) and `combinatorics`.

## 6. References

- Hatcher, A. (2002), *Algebraic Topology*, Cambridge University Press.
- Massey, W. S. (1991), *A Basic Course in Algebraic Topology*, Springer.
- Spanier, E. H. (1966), *Algebraic Topology*, McGraw-Hill.
- Munkres, J. R. (1984), *Elements of Algebraic Topology*, Addison-Wesley.
- Rourke, C. P. and Sanderson, B. J. (1972), *Introduction to Piecewise-Linear Topology*, Springer.
- Whitehead, J. H. C. (1939), Simplicial spaces, nuclei and m-groups, *Proc. London Math. Soc.* 45, 243–327.
- Zeeman, E. C. (1964), On the dunce hat, *Topology* 2, 341–358.
- Cohen, M. M. (1973), *A Course in Simple-Homotopy Theory*, Springer.
- Gallier, J. and Xu, D. (2013), *A Guide to the Classification Theorem for Compact Surfaces*, Springer.
- Pezzana, M. (1974), Sulla struttura topologica delle varietà compatte, *Atti Sem. Mat. Fis. Univ. Modena* 23, 269–277.
- Ferri, M., Gagliardi, C. and Grasselli, L. (1986), A graph-theoretical representation of PL-manifolds, *Aequationes Math.* 31, 121–141.
- Gurau, R. (2011), The 1/N expansion of colored tensor models, *Ann. Henri Poincaré* 12, 829–847.
- Bonzom, V., Gurau, R., Riello, A. and Rivasseau, V. (2011), Critical behavior of colored tensor models in the large N limit, *Nucl. Phys. B* 853, 174–195.
- Reisenberger, M. and Rovelli, C. (1997), "Sum over surfaces" form of loop quantum gravity, *Phys. Rev. D* 56, 3490–3508.
- Oriti, D. (2001), Spacetime geometry from algebra: spin foam models for non-perturbative quantum gravity, *Rep. Prog. Phys.* 64, 1489–1544.
- Dowker, C. H. (1952), Homology groups of relations, *Ann. of Math.* 56, 84–95.
- Atkin, R. H. (1974), *Mathematical Structure in Human Affairs*, Heinemann, ch. 3.

## 7. Open Questions

1. **Cup products and cohomology operations,** for the intersection forms and characteristic classes that gauge theory and quantum gravity ask for.
2. **Collapses in dimension 3 and above,** and discrete Morse theory, which need regular complexes or more data than incidences.
3. **Recognition of the 3-sphere** (Rubinstein and Thompson), to make manifold recognition exact in dimension 4.
4. **Non-orientable coloured graphs** (non-bipartite), and the coloured graphs of a given manifold (crystallisations of minimal size).
5. **Products and joins of simplicial complexes, and quotients by free actions,** which physis's `quantum_gravity` checks now build for themselves (staircase products ∂Δ³ × Cₙ for S² × S¹, joins of cycles and their quotients for lens spaces); and a test that a coloured graph encodes a manifold (all its bubbles spheres).
