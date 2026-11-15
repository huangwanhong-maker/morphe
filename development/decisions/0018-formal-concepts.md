# 0018. Formal Contexts, their Concepts and their Implications

- **Status:** Accepted, provisionally (2026-11-15).
- **Date:** 2026-11-15
- **Decided by:** Claude, for the foundation of case-oriented inquiry zetesis (the infrastructure's ADR 0117 and zetesis's decision 0003), at the user's request for traditional analysis beneath policy making.

## 1. Context

zetesis reads a concept's intension and extension as the two sides of a Galois connection between cases and attributes, Sartori's ladder of abstraction as a chain in the concept lattice, and the attributes that add nothing given others as implications. The construction names no object and is order theory, so it belongs to morphe.

## 2. Decision

`combinatorics/concepts.hpp`:

- formal contexts;
- the derivation operators and the closure;
- every concept by Ganter's NextClosure, in lectic order;
- the concept lattice with its covering relation (this module's `reduction` of the strict order of extents);
- the canonical basis of implications of Guigues and Duquenne, by NextClosure under the pseudo-closure;
- forward closure under implications.

Lists longer than a stated bound are refused, never cut.

*Remark (morphisms).* A context's concepts depend only on it up to the clarification and reduction of objects and attributes. Morphisms of contexts (scale measures, infomorphisms) are not written yet.

## 3. Consequences

`check_morphe_concepts` (9):
- on 120 random contexts, the concepts equal the closed sets found among every subset, in lectic order;
- meets and joins are intersections, and the covers are found by brute force;
- the stem base is sound, complete and minimal, with one implication per pseudo-intent;
- the contranominal and ordinal scales have 2ⁿ and n concepts.

morphe has 821 checks.
