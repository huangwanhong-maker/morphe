# 0004. The First Slice: Exact Integers, Linear Algebra over ℤ, and Cell Complexes

- **Status:** Accepted, provisionally (the data model in §2.4 answers open decision 2 for complexes only)
- **Date:** 2026-10-23
- **Decided by:** the user asked for morphe and physis to be written in parallel, with their structure planned and only what is needed now implemented: first, the analysis and visualization of 2-complexes and CW complexes, for spacetime models, field theory (QFT and group field theory) and the study of symmetry and gauge. Claude designed this slice and wrote this record; physis's part is its decision 0005, and the composition is The Science Infrastructures' [ADR 0039](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0039-complexes-symmetry-and-fields.md).

## 1. Context

morphe had a plan and no code. Its open decision 1 offered two first slices; the first, arithmetic → algebra → topology, is the one now needed. A 2-complex is the arena of much of what the user studies:

- a CW complex of dimension 2 carries the fundamental group (by its attaching words) and the first two homology groups;
- the Feynman graphs of group field theory are 2-complexes (stranded graphs, dual to triangulations), and their amplitudes are spin-foam amplitudes on them;
- a lattice gauge field is a function from the 1-cells to a group, and its curvature is read on the 2-cells;
- a discrete spacetime is a simplicial or cellular complex with a geometry on it.

## 2. Decision

### 2.1 Four Modules Written Now, Header-Only, Status `draft`

| Module | Written now | Planned later in the module |
|---|---|---|
| `arithmetic` | `Integer` (arbitrary precision, held in a machine word while it fits); `Rational`; residues modulo a prime | algebraic numbers, polynomial rings, factorisation, lattices of ℤⁿ |
| `algebra` | exact matrices over ℤ, dense and sparse; the Smith normal form with its transforms; ranks over ℚ and over 𝔽_p; chain complexes, their homology and cohomology with generators; finitely presented groups (words, free and cyclic reduction, Tietze simplification, abelianisation); finite groups by their tables (cyclic, dihedral, symmetric, alternating, quaternion, products), their conjugacy classes, and the homomorphisms from a presented group into one | Hermite normal form, permutation groups by Schreier–Sims, character tables (Dixon–Schneider), coset enumeration (Todd–Coxeter), rewriting (Knuth–Bendix), Lie algebras over exact rings, universal algebra |
| `combinatorics` | multigraphs, their components, spanning forests and paths in them | as decision 0003 plans: C-sets, flows, matchings, isomorphism, Feynman-graph enumeration with automorphisms |
| `topology` | CW complexes of any dimension, with attaching words for their 2-cells; simplicial complexes (faces, links, stars, skeleta, subdivision, joins); a catalogue of constructions; cellular homology and cohomology with generators; fundamental groups as presentations; the classification of closed surfaces; elementary collapses; the dual 2-complex of a triangulated pseudomanifold; coloured graphs (crystallisations) and the pseudomanifolds they encode | cubical complexes, simplicial sets, persistence, cup products, discrete Morse theory, 3-manifold recognition, knots |

### 2.2 The Plan of Topology

`topology` is organised by what it holds and what it computes:

1. **Complexes:** simplicial (abstract, ordered), CW (cells with attaching data), cubical, Δ-complexes and simplicial sets; products, joins, cones, suspensions, quotients, subdivisions; duals.
2. **Invariants:** Euler characteristic; homology and cohomology over ℤ, ℚ and 𝔽_p, with generators; the fundamental group; later the cup product, torsion invariants (Reidemeister) and persistent homology of filtrations.
3. **Homotopy:** elementary collapses and expansions (simple homotopy), discrete Morse theory, homotopy equivalence of complexes.
4. **Manifolds:** recognition of combinatorial manifolds and pseudomanifolds; surfaces classified; 3-manifolds by normal surfaces and recognition of the sphere (research); coloured graphs as encodings of manifolds.
5. **Knots and links** (planned): diagrams, Reidemeister moves, polynomial invariants.

### 2.3 The Plan of Algebraic Geometry

`algebraic_geometry` stays `planned`, with this order of work when a domain needs it:

1. **Polynomials** in `arithmetic`: sparse multivariate polynomials over ℤ, ℚ and 𝔽_p, monomial orders.
2. **Ideals:** Gröbner bases (Buchberger, then F4), ideal membership, elimination, saturation and quotients.
3. **Invariants:** dimension and degree, Hilbert functions and series, free resolutions (Schreyer).
4. **Varieties:** affine and projective varieties, their points over finite fields, morphisms; toric varieties from polytopes and fans (with `combinatorics`).
5. **Real algebraic geometry:** Sturm sequences, real root isolation, cylindrical algebraic decomposition.
6. **Tropical geometry** and, for research, schemes and sheaf cohomology as far as computation reaches.

### 2.4 Data Model for Complexes

A CW complex is held as explicit arrays, not as a generic C-set (open decision 2): its 0-cells are counted; each 1-cell has a tail and a head; each 2-cell has a base vertex and an attaching word, a closed edge path; each higher cell has its cellular boundary as integer incidences. The attaching words are what a C-set over a plain schema cannot hold (as the idea of 2026-10-21 foresaw), and they are what the fundamental group and the holonomy of a gauge field need. The idea stays a seed, for graphs and meshes.

### 2.5 Exact Arithmetic of Our Own

Open decision 4 is answered for now: morphe writes its own integers (§10.3 of the infrastructure's map: control and learning), checked against identities and against machine arithmetic where it does not overflow. FLINT and GMP stay unlinked; a back end can sit behind `Integer` later.

### 2.6 Checks without calx

morphe depends on no project library (invariant 1), so its checks cannot use calx's `verification`. A header of its own, `morphe/check.hpp` at the library's root, prints the same report (`ok`, `FAIL`, `skip`, the count), so the infrastructure's runner adds morphe's checks to its total. Its build is a CMake project of its own, like calx's.

### 2.7 Algorithms with Their Limits

Following invariant 4:

- the Smith normal form terminates on every integer matrix; on a sparse boundary matrix it first eliminates unit pivots (an algebraic reduction of the chain complex, Kaczynski, Mischaikow and Mrozek 2004, ch. 4) and then reduces the dense remainder, whose size governs its cost;
- Tietze simplification always terminates and returns an isomorphic presentation, but it is a heuristic: the isomorphism problem for finitely presented groups is undecidable, and no normal form is claimed;
- counting homomorphisms into a finite group is exact, by backtracking over the generators, and exponential in their number in the worst case;
- recognition of surfaces is exact; recognition of manifolds of dimension 3 is limited to the local conditions (links of vertices are spheres by their homology and Euler characteristic), which characterise manifolds only with the Poincaré conjecture and a check of simple connectivity left open.

## 3. Alternatives Considered

- **The other slice first** (categories → equivalence → formal): it serves research in formal languages, which no work presently asks for.
- **Generic C-sets now:** the abstraction would come before its second instance; the complexes need attaching words in any case.
- **FLINT or GMP:** both LGPL, needing a decision of their own; the integers needed now are small, and writing them is how the method is learnt.
- **Checks on calx's `verification`:** it would break morphe's independence for its tests alone.

## 4. Consequences

- morphe has its first C++: four modules, a root `CMakeLists.txt`, and its own check report. The infrastructure's root `CMakeLists.txt` adds morphe before physis.
- physis builds its gauge fields, spin foams and discrete spacetimes on these complexes (physis decision 0005).
- morphe's README, architecture and CLAUDE.md list what is written; the glossary (`docs/glossary.md`) is created with this slice and fixes the order of composition and the side of actions.
