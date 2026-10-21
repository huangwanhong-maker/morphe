---
title: C-Sets as the Common Data Model
status: seed
date: 2026-10-21
related: [categories, topology, formal, equivalence]
---

# C-Sets as the Common Data Model

## 1. Idea

Represent combinatorial structures uniformly as C-sets: functors X: C → **Set** from a small category C (the schema) to sets.

- **Graphs.** A graph is a C-set whose schema has two objects, E and V, and two arrows, src and tgt: E → V.
- **Other structures.** Semi-simplicial sets, cell complexes, Petri nets, wiring diagrams and meshes are C-sets over other schemas.
- **Attributes.** Attributed C-sets add data such as coordinates or labels (Patterson, Lynch and Fairbanks, 2022; Catlab.jl).

## 2. Consequences

- **Shared machinery.** One generic implementation serves every schema. It provides:
  - morphisms;
  - limits and colimits, including gluing by pushout;
  - isomorphism testing.
- **Meshes for physis.** Meshes become attributed C-sets, so discrete exterior calculus can run on them. CombinatorialSpaces.jl already does this.
- **Rewriting.** Rewriting of structures becomes double-pushout rewriting on C-sets, which connects to `formal/`.

## 3. Questions

- **Performance.** How does a generic C-set representation in C++ compare with hand-written mesh data structures?
- **Coverage.** Does every structure we need fit? CW complexes, for example, carry incidence degrees and may need attributes beyond plain C-sets.
