# 0021. Normal Forms of Traces

- **Status:** Accepted, provisionally (2026-11-20).
- **Date:** 2026-11-20
- **Decided by:** Claude, for the infrastructure's framework for allocating heterogeneous resources by evolution, in which the order of grants matters except where grants commute (its design `development/ideas/2026-11-20-evolutionary-allocation-of-heterogeneous-resources.md`, Stage 2), and heuresis's `sequential`, which evaluates each order once up to exchanges of commuting grants.

## 1. Context

`formal/traces.hpp` (decision 0008) gives a record's causal order, its Foata steps and the number of its linearisations, for events that declare what they read and write. A search over orders of grants needs the same for words over an alphabet with an independence relation, the trace monoid of Mazurkiewicz, and above all a canonical form: one word for each trace, so that two orders equal up to exchanges of commuting grants are recognised as one and evaluated once.

## 2. Decision

In `formal/traces.hpp`:

- `Independence`, letters 0…k−1 with a symmetric, irreflexive relation, built from the commuting pairs or from each letter's footprint, and checked whenever it is used;
- `causal_order(word, I)`, sharing the closure with the causal order of a record;
- `lexicographic_normal_form`, the least word of the trace, by Anisimov and Knuth's greedy choice of the least available letter;
- `foata_normal_form`, the steps written in increasing order of their letters;
- `equivalent`, by the lexicographic normal form, and `projections_agree`, by the projection lemma of Cori and Perrin, written apart so that each checks the other.

The checks hold all three to the traces found by exchanging adjacent independent letters breadth-first, on every word up to length five over forty random independence relations, and each trace's size to `linearisations`.

## 3. Alternatives Considered

- **The Foata form as the only key.** Equally canonical; but the lexicographic form is itself a word of the trace, which a search can use as the representative it evaluates.
- **Leaving words to heuresis.** The normal forms name no search; they are a property of the trace monoid, and morphe already holds the causal order they rest on.

## 4. Consequences

- 11 more checks in `check_morphe_traces` (19).
- heuresis's `sequential` (its decision 0009) keys its cache of evaluations by the lexicographic normal form.
- Not yet: the trace monoid's Möbius function and generating series, and normal forms for traces with an infinite or partial alphabet.
