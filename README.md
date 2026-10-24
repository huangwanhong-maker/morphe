# morphe

*Mathematical structures and their equivalences.* The name comes from the Greek *μορφή*, "form", the root of "morphism".

## 1. Purpose

morphe computes with mathematical structures, whether algebraic, topological, categorical or formal. It also computes with the morphisms and equivalences between them. It serves two roles:

- **A foundation for computation elsewhere:**
  - exact arithmetic;
  - groups and representations, for symmetry;
  - cell complexes and homology, for the meshes of physical models.
- **An instrument of research in its own right:**
  - algebraic topology, algebraic geometry and category theory;
  - formal languages and automata;
  - notions of structural sameness, such as isomorphism, homotopy equivalence, bisimulation and bimorphism.

morphe is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures), but it stands on its own.

## 2. Planned Modules

| Module | Content | Status |
|---|---|---|
| `logic` | Propositional, first-order and modal logic; set theory as a foundation; computability; the simply and dependently typed λ-calculi, proof terms and normalisation (Curry–Howard–Lambek); satisfiability | planned |
| [`arithmetic`](arithmetic/) | ℤ, ℚ, finite fields, algebraic numbers, polynomial rings; computational number theory | draft: integers, rationals, residues |
| `categories` | Categories, functors and natural transformations; limits and colimits; monoidal structure; algebras and coalgebras of functors | planned |
| `equivalence` | Generic algorithms for sameness: canonical forms, partition refinement for simulation and bisimulation, morphisms between coalgebras | planned |
| [`algebra`](algebra/) | Groups (finite, permutation, finitely presented), rings, modules, representations; orders and lattices; universal algebra. Homological algebra: chain complexes, Smith normal form. Exact linear algebra over ℤ, ℚ and finite fields: elimination, Hermite and Smith normal forms, kernels, determinants without rounding. | draft: Smith normal form, homology with generators, presentations, finite groups, determinants |
| [`combinatorics`](combinatorics/) | Enumeration, generating functions, designs; graphs and hypergraphs as C-sets (paths, cuts, flows, matchings, colourings, isomorphism by refinement); matroids; exact combinatorial optimisation | draft: multigraphs, spanning forests, canonical forms and automorphism groups, Eulerian circuits, the matrix-tree and BEST theorems |
| [`topology`](topology/) | Simplicial, cubical and CW complexes. Homology and cohomology. Persistent homology. Fundamental groups. Homotopy equivalence. Surfaces and combinatorial manifolds; dual complexes; coloured graphs. | draft: CW and simplicial complexes, homology, π₁, surfaces, collapses, duals, coloured graphs, persistence and cycles over ℤ/2, Q-analysis, worldline complexes |
| `algebraic_geometry` | Polynomial ideals and Gröbner bases; elimination; dimension, Hilbert series and resolutions; varieties over ℚ and finite fields; toric varieties; real algebraic geometry; tropical geometry (the plan of [decision 0004](development/decisions/0004-the-first-slice-cell-complexes.md), §2.3) | planned |
| [`formal`](formal/) | Words and trees; automata, grammars and transducers; rewriting systems; transition systems. Bimorphisms in the sense of Arnold and Dauchet. | draft: words, Lyndon words, suffix arrays, the BWT, FM-indexes, Aho–Corasick; semirings and the alignment of words over them; the polytope algebra ([decision 0005](development/decisions/0005-words-their-indexes-and-their-alignment.md)) |

## 3. Principles

These are set out in full in [docs/architecture.md](docs/architecture.md).

- **Structures with morphisms.** Every structure comes with its morphisms, and every notion of equivalence names the category it lives in.
- **A categorical spine.** Category theory supplies the shared constructions.
- **Two directions of description.** Algebras describe how things are built; coalgebras describe how they behave. Bisimulation unifies equivalence for automata, processes and dynamical systems.
- **Exact by default.** Every algorithm also states when it terminates and how costly it is.

## 4. Relation to Sibling Libraries

| Library | Scope | Relation to morphe |
|---|---|---|
| [calx](https://github.com/huangwanhong-maker/calx) | Retargetable numerical computing | Independent of morphe |
| [physis](https://github.com/huangwanhong-maker/physis) | Geometry and theoretical physics | Depends on morphe |

morphe depends on no other library of the project.

## 5. Status

**Draft (2026-10-23).** The first slice is written ([decision 0004](development/decisions/0004-the-first-slice-cell-complexes.md)): exact integers and rationals, linear algebra over ℤ with the Smith normal form, the homology of chain complexes with its generators, presented and finite groups, and the topology of cell complexes: CW and simplicial complexes, their homology and fundamental groups, the classification of surfaces, collapses, dual 2-complexes and coloured graphs. Four modules, header-only, with 98 checks; the bubbles of invariant field theories and their Feynman graphs since physis decision 0006. The rest is planned; [decision 0003](development/decisions/0003-logic-and-combinatorics.md) added `logic` and `combinatorics` to the plan, and the place of each module in the whole of the project's mathematics is the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

Build and check with CMake 3.20 or later; morphe needs no other library:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 6. License

morphe is released under the MIT License. See [`LICENSE`](LICENSE).
