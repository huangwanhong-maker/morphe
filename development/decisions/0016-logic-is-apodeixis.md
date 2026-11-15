# 0016. Logic Is apodeixis's

- **Status:** Accepted, provisionally (2026-11-14).
- **Decided by:** the user, who asked for a foundation of formal logical reasoning with derivations and their verification, and chose two foundations, apodeixis for logic and tekmerion for evidence (the infrastructure's ADR 0112), and Claude.

## 1. Context

morphe planned a module `logic` from decision 0003: propositional, first-order and modal logic; set theory as a foundation; computability; the typed λ-calculi with proof terms and normalisation; satisfiability. Nothing of it was written. The infrastructure has now founded apodeixis, on morphe and calx, for logic and verification: formulas and sequents read by morphe's grammars in BNF, calculi declared as rule schemata and derivations checked against them, decision procedures with countermodels, satisfiability with checked refutations, argumentation, and derivability relative to an interpretation of a norm.

## 2. Decision

1. **Logic is apodeixis's.** morphe's planned `logic` row points there. apodeixis stands on morphe and uses its public headers (grammars in BNF with Earley's parser, directed graphs, acyclic relations and the closure of rules).
2. **Nothing moves.** morphe keeps grammars, parsers, rewriting and the closure of rules over several premises (decision 0015), which apodeixis uses and does not write again.
3. **What remains open.** Set theory as a foundation, computability and the typed λ-calculi are planned in apodeixis; if a part of them is one day wanted beneath morphe's own modules, it is proposed here, and morphe still depends on nothing in the programme.

## 3. Alternatives Considered

- **Write `logic` here, as planned.** Set aside at the user's word, for the reasons of decision 0013: the body is large, of another kind than morphe's structures, and needs morphe's grammars beneath it.
- **Drop the row.** It would hide where logic is; it is kept, pointing to apodeixis.

## 4. Consequences

- `formal`'s planned use of `logic` (docs/architecture.md, §2) is read as a use of apodeixis by those who need both; morphe's own modules need no logic.
