# 0006. Representations of Finite Groups, and Grammars with Symmetry

- **Status:** Accepted
- **Date:** 2026-10-25
- **Decided by:**
  - **The request.** The user gave Claude their draft, *Towards a Group Context-Sensitive Grammar for Genome Analysis* (2026), and asked it to implement it first, then to survey advanced genome-analysis algorithms and implement them.
  - **The guidance.** They added that the draft's formalism is preliminary: where it cannot be reproduced as written, the implementation should follow the idea's direction and be correct.
  - **What Claude did.** It placed the general mathematics here and wrote the action, grammar and parser headers and this record. Agents wrote the representations, the cyclotomic numbers, the r-index and the positional BWT to the contracts given them.
  - **The composition.** The infrastructure's [ADR 0047](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0047-group-grammars-and-advanced-genome-analysis.md).

## 1. Context

The draft makes a finite group G the primitive of a grammar of genomes. It derives the alphabet from G's group algebra ℂ[G] and its Wedderburn decomposition, parses with rules taken up to G, and obtains concrete letters by symmetry-breaking maps. Everything in it that names no object in nature is mathematics, and belongs here:
- groups acting on sets and words;
- characters, group algebras and Fourier transforms;
- grammars and their parsers;
- invariant grammars, and maps between alphabets.

morphe's plan already gave grammars to `formal` and representations to `algebra`; neither was written. Its architecture anticipated this use (§1.6: "work on group-equivariant grammars, for example, can use it directly").

The draft's own statements needed care. Some claims hold as written. Some hold only when restated: the bound on the number of rule orbits runs the other way, and the per-position orbit scan decides a larger language. A few do not hold: the group algebra as an algebra does not determine the group. The infrastructure's note `domains/genomics/theory/group-grammars.md` gives the precise forms, with proofs, and the amendments. The code here implements the precise forms.

## 2. Decision

1. **`algebra` gains the representation theory of finite groups.** Everything is exact.
   - **Actions** (`action.hpp`): actions on finite sets, with a sign under which elements act on words by anti-automorphisms (reversal); orbits, stabilisers; Burnside's count; canonical forms of words; the orbits of Σⁿ by Burnside's lemma.
   - **Cyclotomic numbers** (`arithmetic/cyclotomic.hpp`): ℚ(ζ_m) exactly, with Φ_m, inverses, Galois automorphisms and lifting.
   - **Exact linear algebra** over any exact field (`field_matrix.hpp`).
   - **Character tables** of any finite group by the Dixon–Schneider method (`characters.hpp`).
   - **Group algebras** K[G] (`group_algebra.hpp`): central idempotents, Wedderburn blocks measured by rank, isotypic components, rational blocks, and ℚ[C_n] ≅ ⊕_{d|n} ℚ(ζ_d).
   - **Representations of S_n** in Young's seminormal form, with Murnaghan–Nakayama characters, hooks, the Fourier transform on S_n and its matrix units (`symmetric_representations.hpp`).
   - **Abelian transforms**: the exact Fourier transform of finite abelian groups, and the Walsh–Hadamard transform over ℤ (`abelian_fourier.hpp`).
2. **`formal` gains grammars** (`grammar.hpp`, `earley.hpp`, `inside_outside.hpp`).
   - Phrase-structure grammars, context-free and monotone.
   - Earley's parser with ε-rules and parse trees.
   - The inside and outside algorithms over any semiring of `semiring.hpp`, on a binary form that is a bijection on derivations, with expected rule counts and Viterbi trees.
   - The search of a monotone grammar's forms.
3. **`formal` gains grammars with symmetry** (`equivariant.hpp`, `breaking.hpp`), in the precise form.
   - **Invariant grammars** under actions with signs, and their rule orbits checked against Burnside's count.
   - **An equivariant Earley parser** that binds each rule's group element as a coset narrowed at each scan; it decides exactly the language. The draft's orbit-by-position scan is kept as an option, so that the larger language it decides, π⁻¹(L(Γ/G)), can be shown.
   - **Shared spans**: spans of one class share their values.
   - **The orbit-canonical search** of monotone grammars.
   - **Letter-to-letter maps**: image grammars, the breakings of a torsor, and coarse-to-fine masks that keep every complete parse exactly.
4. **`formal` gains two indexes** for collections of words: the run-length FM-index and the r-index (`r_index.hpp`), and the positional BWT (`pbwt.hpp`).
5. **`combinatorics` gains increasing subsequences** (`increasing.hpp`), for the chaining of anchors in whole-genome comparison.
6. **`formal` now depends on `algebra`**, for groups and their actions. The order stays acyclic, since `algebra` depends on `arithmetic` alone.

## 3. Alternatives Considered

- **The draft's formalism as written.** Rejected where it fails: a parser built on its scan rule accepts words outside the language. The checks show this on S → aa | bb, and on random invariant grammars.
- **Numerical characters and group Fourier transforms**, in calx. Rejected for the identities: invariant 3 keeps floating point out of algebraic identities. The exact tables, idempotents and transforms are here, and a domain evaluates them numerically where its data are real (the genomics domain's isotypic channels).
- **Characters of named families only** (abelian groups and S_n). Rejected in favour of a general method: the draft's open question Q2 asks for groups discovered from data, which can be any group.
- **Grammars in a module of their own.** Rejected: morphe's plan gives automata, grammars and transducers to `formal`, beside the words they act on.

## 4. Consequences

- The genomics domain's `grammar` module, and the infrastructure's Genome Workbench, use these headers. The domain writes no general method of its own.
- The checks:
  - representations: cyclotomic, characters, group algebra, symmetric, abelian Fourier, actions;
  - grammars: 31;
  - indexes: the r-index and the PBWT.
- **Open.**
  - Weighted automata and transducers.
  - Grammars beyond monotone membership.
  - Representations of large groups, whose character tables need memory for the multiplication table.
  - Irreducible representations as matrices beyond S_n and abelian groups.
  - The r-index built in O(r) space, and compressed PBWT columns.
