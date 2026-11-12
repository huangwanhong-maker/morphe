# 0013. Categories as Objects of Study Are kategoria's

- **Status:** Accepted, provisionally (2026-11-12).
- **Decided by:** the user, who asked for a foundation of category theory standing on its own (the infrastructure's ADR 0100), and Claude.

## 1. Context

morphe planned a module `categories` from its first day: categories, functors and natural transformations; limits and colimits; monoidal structure; algebras and coalgebras of functors. Its architecture gives that module two jobs at once. One is a *spine* beneath morphe's other modules: generic constructions (products, quotients, functors between kinds of structure) written once and not again for each structure. The other is category theory as a subject: computing with categories as objects of study.

Nothing of either was written. The infrastructure has now founded kategoria for the second.

## 2. Decision

1. **Categories as objects of study are kategoria's:** finite categories and those presented by generators and equations, functors, natural transformations, equivalences, adjunctions, monads, limits and colimits, set-valued functors and Kan extensions, monoidal theories and their diagrams, and the deciding of claims about them. kategoria stands on morphe and uses its public headers (canonical forms of labelled graphs, grammars in BNF with Earley's parser, permutation groups).
2. **morphe's `categories` is narrowed to the spine:** the generic constructions that two or more of morphe's own modules share. It stays planned, and is written when a second module needs a construction the first already has. It does not depend on kategoria: morphe depends on nothing in the programme.
3. **No structure moves.** Groups, graphs, complexes, words, grammars and rewriting stay in morphe with their algorithms. kategoria reads them as categories where it needs them.

## 3. Alternatives Considered

- **Write all of it here, as planned.** One foundation fewer. Set aside at the user's word; and the subject, with a language of claims and its session, is a body of its own that needs morphe beneath it.
- **Drop `categories` from morphe.** The spine is still wanted: a product of two structures should be written once. It is kept, narrowed.
- **Let morphe depend on kategoria.** A cycle, since kategoria needs morphe's algorithms.

## 4. Consequences

- The README's table and the map of mathematics say where category theory is.
- The planned order of morphe's modules is unchanged: `categories`, as the spine, still stands alone at the bottom.
- When morphe's structures are wanted as categories (the category of a group's actions, the nerve of a category as a simplicial set), the construction is kategoria's, on morphe's types.
