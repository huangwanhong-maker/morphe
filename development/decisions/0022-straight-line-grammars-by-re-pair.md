# 0022. Straight-Line Grammars by Re-Pair

- **Status:** Accepted, provisionally (2026-11-26).
- **Date:** 2026-11-26
- **Decided by:** Claude, for syntaxis's grammars of foams (its decision 0004), at the infrastructure's request to learn a generative grammar of a city's event spacetime: the lives of its buildings and streets, read along their worldlines, are words, and their repeated patterns are rules.

## 1. Context

A straight-line grammar derives one word from each of its nonterminals; built greedily by Re-Pair it is a compressed form of a set of words whose rules are their most repeated pairs, nested. syntaxis's induction compresses geometric graphs by their own patterns; words are morphe's, so their grammar-based compression belongs here.

## 2. Decision

`formal/include/morphe/formal/repair.hpp`: `StraightLine` (the first new letter, the rules, the words left, the input's length; `derives`, `expand`, `length`, `size`) and `re_pair(words, max_rules)`: Larsson and Moffat's Re-Pair over a set of words, occurrences counted without overlap (⌊k/2⌋ of aa in a run of k), pairs never spanning two words, ties to the least pair, stopping when no pair occurs twice.

## 3. Alternatives Considered

- **Sequitur** (Nevill-Manning and Witten 1997). It reads online and keeps digram uniqueness and rule utility as invariants; Re-Pair's offline greed gives smaller grammars in practice and an order of rules by frequency, which is what a reader of the rules wants first.
- **The smallest grammar.** Finding it is NP-hard (Charikar et al. 2005); Re-Pair is the usual approximation.

## 4. Consequences

`check_morphe_formal_repair`, 11 checks: abababab to BB with B → AA → abab; a run of five equal letters; pairs never spanning words; ties to the least pair; on 300 random sets of words, exact expansion, digram uniqueness, rules made of earlier letters, a grammar never larger than its input; the bound on rules.
