# 0003. Logic and Type Theory, and Combinatorics and Graphs, in the Plan

- **Status:** Accepted
- **Date:** 2026-10-23
- **Decided by:** the user asked that the project's mathematics be designed as a whole and comprehensively (topology, geometry, category and type theory, differential equations, linear algebra, statistics, numerical computation), and written only as it is needed. Claude designed morphe's part and wrote this record. The whole map is The Science Infrastructures' [ADR 0035](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0035-the-mathematical-foundations.md) and [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), checked against every class of the Mathematics Subject Classification 2020.

## 1. Context

morphe's plan (decision 0002) covered categories, equivalence, arithmetic, algebra, topology, algebraic geometry and formal languages. Measured against the Mathematics Subject Classification, two areas of exact mathematics had no home: logic and the foundations (MSC 03), with type theory, which the user named; and combinatorics (MSC 05), with graph theory, which networks of every science need. Exact linear algebra (MSC 15, over ℤ, ℚ and finite fields) was implicit in homological algebra and is now named.

## 2. Decision

1. **`logic`**: propositional, first-order and modal logic; set theory as a foundation; computability; the simply and dependently typed λ-calculi, proof terms and their normalisation (the Curry–Howard–Lambek correspondence with `categories`); satisfiability. It depends on nothing; `formal` may use it.
2. **`combinatorics`**: enumeration, generating functions and designs; graphs and hypergraphs as C-sets (the data model of the idea of 2026-10-21); paths, cuts, flows, matchings, colourings, isomorphism by partition refinement; matroids; exact combinatorial optimisation. It depends on `categories` and `equivalence`; `topology` may use it.
3. **`algebra`** names exact linear algebra (elimination, Hermite and Smith normal forms, kernels and images, determinants without rounding) and orders, lattices and universal algebra. Its numerical counterpart is calx's `linalg`.
4. No code is written now: nothing the first simulation needs is exact. Modules are written when a domain or a sibling needs them, first `arithmetic`, `algebra` and `topology` toward discrete exterior calculus in physis.

## 3. Alternatives Considered

- **Type theory inside `formal`** (with rewriting and λ-terms). Type theory is logic and the foundation of proof, with its own semantics in categories; `formal` keeps words, trees, automata, grammars and rewriting, and uses `logic`.
- **Graphs inside `topology`** (as one-dimensional complexes). Graph algorithms (flows, matchings, isomorphism) are combinatorics, and a graph is the simplest C-set; `topology` uses `combinatorics`, not the reverse.

## 4. Consequences

- morphe's plan has nine modules; its README, architecture and CLAUDE.md list them with their order.
- The C++ conventions of the infrastructure's ADR 0035 apply when the first module is written.
