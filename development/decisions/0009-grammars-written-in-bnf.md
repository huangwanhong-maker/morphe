# 0009. Grammars Written in BNF over Bytes, and the Diagnostics of Earley's Chart

- **Status:** Accepted, provisionally (left-recursive lists are read back by recursion, to a depth of the derivation's height; the grammar's ambiguity is not decided)
- **Date:** 2026-11-01
- **Decided by:**
  - **The request.** The user asked for the international systems that describe abnormal genomes: HGVS for sequence variants, ISCN for chromosomes. Both are notations defined by their standards, so both are languages.
  - **The rule.** General methods are written in the foundations.
  - **What Claude did.** It wrote the header, extended Earley's parser, and wrote the check and this record.
  - **The composition.** The infrastructure's [ADR 0064](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0064-abnormal-genomes-their-description-and-expression.md), whose genomics modules `nomenclature` and `karyotype` write HGVS's and ISCN's grammars in this notation.

## 1. Context

`formal` held grammars and Earley's parser over integer alphabets (decision 0006). A notation people write is a language over bytes, and its standard often gives its grammar in an extended Backus–Naur form. Writing such a grammar as numbered rules by hand would be error-prone.

A reader of a rejected description also needs to know two things: how far the text is right, and what could follow. Earley's chart holds that exactly, but the parser did not report it.

## 2. Decision

`formal/bnf.hpp`:

| Name | What it is |
|---|---|
| `compile_bnf(text)` → `Bnf` | rules `name ::= alternatives ;`; literals, classes of bytes, groups, `?`, `*`, `+`; `#` comments; hidden rules (`_name`); errors as `BnfError` with line and column, unproductive rules among them |
| `Bnf` | the compiled `Grammar` and its Earley parser; `atoms`, `atom(byte)`; `parse(input)`, `parse(input, rule)`, `accepts` |
| `ParseNode` | the named tree, with each node's span and text; `child`, `all`, `find`, `find_all` |
| `BnfParse` | `accepted`, `tree`, `furthest`, `expected` (literals from where the text left them, classes as written), `error` |

The terminals are the atoms of the Boolean algebra of byte sets that the grammar's classes and literal letters generate. A class of one atom is one terminal.

`earley.hpp`:

- `parse` takes a start other than nonterminal 0.
- `EarleyResult` reports the furthest viable prefix and the terminals and items awaiting there.

*Rationale.*
- **Atoms.** Without atoms, each class `[0-9]` would be ten rules predicted at every position. With them it is usually a single letter. The partition is the coarsest that makes the compiled grammar equivalent to the text grammar.
- **Diagnostics from the chart.** S_j is non-empty exactly when w₁…w_j is a viable prefix, which is a prefix of a word of the language when every nonterminal is productive. `compile_bnf` rejects unproductive rules, so the furthest prefix and its expected letters are exact rather than heuristic.

## 3. Alternatives Considered

- **A parser generator (LL or LR) or regular expressions.**
  - **Rejected.** The standards' grammars are written for readers, not for a parser class, and LL(1) or LALR(1) would force rewriting them.
  - **Earley instead.** Earley's parser takes them as written, in O(n²) on unambiguous grammars. The notations are short: a description of at most a few hundred bytes.
- **Hand-written recursive-descent parsers in the domains.**
  - **Rejected.** They would hide the grammar in the code, and each domain would write its own.

## 4. Consequences

- **The checks** (`check_bnf.cpp`):
  - **Atoms and trees.**
    - atoms on overlapping classes;
    - trees evaluated on 400 random expressions;
    - hidden rules spliced into their parents;
    - a start other than the first rule.
  - **Agreement with references.**
    - 80 random regular grammars agree with `std::regex` on every word up to length 6;
    - JSON's numbers agree with their regular expression.
  - **Diagnostics.**
    - on 60 finite languages enumerated completely, the furthest viable prefix and the expected letters are exact;
    - the messages are checked;
    - errors in grammar texts are located.
  - **Scale.** A list of 5000 items has a linear chart.
- **The genomics domain** writes HGVS's grammar (`nomenclature`) and ISCN's (`karyotype`) in this notation.
- **Open.**
  - Read the tree back iteratively, for lists of tens of thousands of items.
  - Count parse trees in the chart (a semiring over the forest), so that a grammar's unambiguity can be checked on its examples.
  - Case-insensitive literals.
