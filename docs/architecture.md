# Architecture

*Status: draft, 2026-10-21.* This document describes how morphe is structured. The decisions behind it are recorded in [development/decisions/](../development/decisions/).

## 1. Principles

### 1.1 Structures Come with Their Morphisms

Every structure in morphe is defined together with its morphisms, so every structure lives in a category. Each notion of sameness is therefore relative to a stated category and is never left implicit. Examples are isomorphism, homotopy equivalence, bisimilarity and Morita equivalence.

### 1.2 Category Theory as the Spine

Generic constructions are written once in `categories/` and reused. They are not reimplemented for each structure. Examples are products, quotients, limits and colimits, and functors between kinds of structure.

### 1.3 Algebra and Coalgebra

- **Algebras of a functor describe construction.** Examples are syntax, terms, and groups presented by generators and relations.
- **Coalgebras describe behaviour.** Examples are automata, transition systems, streams and discrete dynamical systems.
- **Bisimulation is the canonical equivalence for coalgebras** (Rutten, universal coalgebra). One theory therefore covers the equivalence of automata, of processes and of dynamical systems.
- **This is a bridge to physics.** A discrete dynamical system is a coalgebra, so the same machinery that compares automata also compares models of physical evolution.

### 1.4 Exact by Default

Computations are exact: integers, rationals, finite fields, algebraic numbers and polynomials. Floating-point numbers appear only where a structure is built from real-valued data, such as a filtration in persistent homology. They never appear in algebraic identities.

### 1.5 Decidability and Complexity Stated

Many questions about structures are either undecidable in general or hard:

- the word problem for finitely presented groups is undecidable in general;
- graph isomorphism is hard.

Every algorithm therefore states three things:

- the class of inputs on which it terminates;
- its complexity;
- what it returns when it cannot decide.

### 1.6 Independence

morphe depends on no other project library. In particular it does not depend on calx, so it stays light enough to use on its own. Work on group-equivariant grammars, for example, can use it directly.

## 2. Modules and Their Dependencies

| Module | Depends on |
|---|---|
| `logic` | — |
| `arithmetic` | — |
| `categories` | — |
| `equivalence` | `categories` |
| `algebra` | `arithmetic`, `categories`, `equivalence` |
| `combinatorics` | `arithmetic`, `algebra`, `categories`, `equivalence` |
| `topology` | `algebra`, `combinatorics`, `categories`, `equivalence` |
| `algebraic_geometry` | `arithmetic`, `algebra` |
| `formal` | `arithmetic`, `categories`, `equivalence`, `logic` |

This is the planned order. Written so far (decision 0004): `arithmetic`; `algebra` on `arithmetic` alone, since `categories` and `equivalence` are not yet written; `combinatorics` on nothing until decision 0005 joined it to `algebra` and `arithmetic` (the matrix-tree and BEST theorems count by determinants over ℤ); `topology` on `algebra` and `combinatorics`; `formal` on `arithmetic` (decision 0005). When `categories` is written, the generic constructions that the first modules write for themselves (products of complexes, for instance) move there.

- **Generic algorithms.** Generic equivalence algorithms live in `equivalence/`. Examples are canonical forms, partition refinement and coalgebra morphisms.
- **Structure-specific algorithms.** An equivalence specific to one kind of structure lives in that structure's module and uses `equivalence/`. An example is homotopy equivalence of complexes.

## 3. Data Model (Provisional)

Combinatorial structures may share one representation: C-sets, that is, functors from a small category to **Set**. Graphs, semi-simplicial sets, cell complexes, Petri nets and meshes would then all be instances of a single generic implementation ([idea](../development/ideas/2026-10-21-c-sets-as-common-data-model.md)).

For now, CW complexes are explicit arrays (decision 0004, §2.4): vertices counted, each edge a (tail, head), each 2-cell a base vertex and an attaching word, each higher cell its incidences. The attaching words carry what a C-set over a plain schema cannot, the fundamental group and the holonomy of a gauge field.

## 4. Prior Art

Study these systems, but do not link to any GPL-licensed one.

| System | Subject | Licence |
|---|---|---|
| GAP | Groups | GPL |
| Macaulay2, Singular | Commutative algebra, algebraic geometry | GPL |
| SageMath | Broad mathematics | GPL |
| GUDHI | Topological data analysis | MIT (parts GPL) |
| Catlab.jl (AlgebraicJulia) | Applied category theory, C-sets | MIT |
| mCRL2 | Process equivalence, bisimulation | Boost |
| Lean 4 Mathlib | Formalised mathematics | Apache 2.0 |
| FLINT, GMP | Exact arithmetic | LGPL |

## 5. Module Layout

```
<module>/
├── README.md                 purpose, status, dependencies, references, open questions
├── CMakeLists.txt            the target morphe::<module>, and its checks
├── include/morphe/<module>/  public headers
├── src/                      implementation (when a module has code worth compiling once)
├── tests/                    checks, printing the shared report through morphe/check.hpp
└── bindings/python/          Python bindings (optional)
```
