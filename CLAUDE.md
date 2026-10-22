# CLAUDE.md

This file guides Claude Code in the morphe repository.

morphe is one of three libraries developed for the [Integrative Astronomy, Geography, Oceanography Laboratory](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory). The laboratory includes it as the git submodule `lib/morphe`. It is a project in its own right, and must make sense without the laboratory.

## 1. Purpose and Scope

morphe computes with mathematical structures and with the morphisms and equivalences between them. Its scope covers:

- exact arithmetic;
- algebra;
- category theory;
- topology;
- algebraic geometry;
- formal structures: automata, grammars, transducers, rewriting systems and transition systems;
- notions of sameness: isomorphism, homotopy equivalence, simulation and bisimulation, and bimorphism.

It is both a foundation for other code and an instrument of research.

Out of scope:

| Topic | Where it belongs |
|---|---|
| Floating-point numerical methods and devices | calx |
| Geometry as physics uses it, units and physical laws | physis |
| Anything about particular bodies, datasets or observations | the laboratory |

## 2. Invariants

1. **No project dependencies.** morphe depends on no other project library. In particular it does not depend on calx, which keeps it light enough to use on its own.
2. **Structures come with morphisms.** Every structure is defined together with its morphisms. Every notion of equivalence names the category in which it holds.
3. **Exact by default.** Floating-point numbers appear only where a structure is built from real-valued data, such as a filtration in persistent homology. They never appear in algebraic identities.
4. **Behaviour is documented.** Every algorithm states:
   - the class of inputs on which it terminates;
   - its complexity;
   - what it returns when it cannot decide. For example, the word problem for finitely presented groups is undecidable in general.
5. **Sources are cited.** Every definition, theorem and algorithm cites its source (book, paper, theorem number).

## 3. Design Philosophy

morphe shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status.** Every module has one of these statuses:
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **Generic constructions are written once.** Products, quotients, limits, colimits and functors between kinds of structure live in `categories/`. Do not write them again for each structure.
- **Abstract late.** Prefer small, composable units. Introduce an abstraction only after its pattern has appeared twice.

## 4. Repository Layout

```
morphe/
├── docs/          architecture
└── development/   decisions · ideas
```

Module directories (`arithmetic/`, `categories/`, `equivalence/`, `algebra/`, `topology/`, `algebraic_geometry/`, `formal/`) are created together with their first content. Until then, each is listed as `planned` in the README.

## 5. Working Procedures

### 5.1 Module Contract

Each module has:

- a `README.md` giving its purpose, status, dependencies, references and open questions;
- public headers under `include/morphe/<module>/`;
- its implementation in `src/`;
- unit tests in `tests/`;
- optional Python bindings in `bindings/python/`.

Module dependencies follow [architecture §2](docs/architecture.md).

### 5.2 Decisions and Ideas

- **Decisions** are numbered records in `development/decisions/`. A record is never rewritten; a later record supersedes it.
- **Ideas** are files named `development/ideas/YYYY-MM-DD-slug.md`. Their status is `seed`, `promoted` or `abandoned`.

### 5.3 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.4 Shared Conventions

calx, morphe, physis and the laboratory keep the following identical:

- the writing style;
- the status vocabulary;
- the conventions for decisions and ideas.

Change any of these in every repository of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures) together: the three libraries, the laboratories and the umbrella.

### 5.5 Work Inside the Laboratory

When editing morphe from the laboratory's `lib/morphe`:

1. Commit in morphe first.
2. Then commit the updated pointer in the laboratory.

## 6. Mathematical Conventions

- **Glossary.** Notation follows standard references and is recorded in `docs/glossary.md`, which is created with the first module.
- **Fixed once, stated in the glossary:**
  - the order of composition, whether applicative (g ∘ f applies f first) or diagrammatic;
  - whether group actions are left or right actions.
- **Strict or up to isomorphism.** State which one each construction is.
- **Mathematical precision.** Write with the precision of a mathematics text: state hypotheses, distinguish definitions from results, and give or cite proofs.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Asides are labelled *Rationale.*, *Remark.* or *Caveat.*
  - Everything is in English.

## 7. Technology

- **Language: C++.** This follows laboratory ADR 0003 and decision 0002 here.
- **Provisional choices:** C++20, CMake, and nanobind for Python bindings.
- **Exact-arithmetic back ends.** Candidates are FLINT and GMP. Both are LGPL, so adopting either needs an explicit decision. Writing our own is the alternative.
- **Licences.** morphe is MIT-licensed. Systems such as GAP, Macaulay2, Singular and Sage are GPL: study them, but do not link to them.

## 8. Open Decisions

1. **First vertical slice.** Two candidates:
   - arithmetic → algebra (Smith normal form) → topology (homology of complexes). physis needs this path for discrete exterior calculus.
   - categories → equivalence → formal (partition refinement, bisimulation), as a research slice.
2. **Data model.** Should C-sets be the common representation? See the [idea](development/ideas/2026-10-21-c-sets-as-common-data-model.md).
3. **Formal proofs.** Should key theorems be proved in Lean 4 with Mathlib, alongside the C++ code?
4. **Exact-arithmetic back end.** FLINT, GMP, or our own.
5. **Front end.** Python bindings only, or also a small language of morphe's own for presenting structures.

## 9. Current State

As of 2026-10-22, only design documents exist; there is no code. The repository is published at `https://github.com/huangwanhong-maker/morphe`.
