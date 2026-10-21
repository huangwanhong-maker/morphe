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

morphe is one of three libraries developed for the [Integrative Astronomy, Geography, Oceanography Laboratory](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory), but it stands on its own.

## 2. Planned Modules

| Module | Content | Status |
|---|---|---|
| `arithmetic` | ℤ, ℚ, finite fields, polynomial rings | planned |
| `categories` | Categories, functors and natural transformations; limits and colimits; monoidal structure; algebras and coalgebras of functors | planned |
| `equivalence` | Generic algorithms for sameness: canonical forms, partition refinement for simulation and bisimulation, morphisms between coalgebras | planned |
| `algebra` | Groups (finite, permutation, finitely presented), rings, modules, representations and lattices. Homological algebra: chain complexes, Smith normal form. | planned |
| `topology` | Simplicial, cubical and CW complexes. Homology and cohomology. Persistent homology. Fundamental groups. Homotopy equivalence. | planned |
| `algebraic_geometry` | Gröbner bases and varieties | planned |
| `formal` | Words and trees; automata, grammars and transducers; rewriting systems; transition systems. Bimorphisms in the sense of Arnold and Dauchet. | planned |

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

**Planned (2026-10-21).** Only design documents exist; there is no code yet.

## 6. License

morphe is released under the MIT License. See [`LICENSE`](LICENSE).
