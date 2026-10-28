# 0007. Graph Rewriting by Double Pushouts, with Relabelling

- **Status:** Accepted, provisionally (matches are injective; application conditions are predicates on the host and the match, not graphs; classes of graphs rest on `combinatorics/canonical.hpp`, exact component by component)
- **Date:** 2026-10-29
- **Decided by:**
  - **The request.** The user asked for a chemistry workbench whose first plug-in executes their manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking* (v1.4), with a serious verification of its formalism.
  - **The rule.** General methods are written in the foundations.
  - **What Claude did.** It wrote the header, its check and this record.
  - **The composition.** The infrastructure's [ADR 0051](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0051-the-chemistry-workbench.md).

## 1. Context

The manuscript's first formalism is a typed graph grammar. Its productions replace nonterminal sites with ports (hyperedge replacement). Its rules are spans L ← K → R applied by double pushout (Ehrig, Pfender and Schneider 1973). It asks that a symmetry group acting on labels carry derivations to derivations.

morphe's `formal` held words, their indexes and their alignment (decision 0005), but nothing on graphs.

Verifying the formalism turned up a gap. Its acid–base rule keeps every atom and changes two oxygens' charges, which a span of label-preserving monomorphisms cannot do. The fix is relabelling in the sense of Habel and Plump (2002): K partially labelled, kept items relabelled by the double pushout.

## 2. Decision

`formal/graph_rewriting.hpp`, in the namespace `morphe::formal::rewriting` (beside the words' `Rule`):

| Name | What it is |
|---|---|
| `Graph` | vertices and edges with labels, loops and parallel edges, directed or not |
| `Rule` | L, R and K's correspondences of vertices and edges; `"*"` on a kept vertex as a label schema; an optional application condition |
| `for_each_match`, `matches`, `count_matches`, `automorphisms` | the injective matches |
| `dangling_ok`, `apply` | the derived graph and the trace of the host's vertices |
| `derivations` | every direct derivation |
| `code`, `isomorphic`, `rule_code` | classes of graphs (component by component) and of rules |
| `Relabelling`, `act`, `compose`, `closed_under`, `equivariant` | the action of label permutations |
| `explore` → `DerivationSpace` | classes by depth, arcs with their multiplicities of matches, derivations of every length counted |
| `orbits` | orbits of classes under generated relabellings |

*Rationale.*
- **Hyperedge replacement** needs no separate machinery. It is the special case whose L is a nonterminal vertex with its port edges, K its attachment vertices, and R the replacement.
- **Equivariance is tested** as an equality of multisets of derived classes, the content of the manuscript's Proposition 4.1. **Closure is tested** separately, because the proposition holds only under it.
- **Classes are encoded per component.** A graph's class is the sorted list of its components' classes, so that graphs of many identical parts (free molecules, released water) keep an exact canonical form.

`formal` now links morphe's `combinatorics`, for the canonical forms.

## 3. Checks

`formal/tests/check_graph_rewriting.cpp` (20 checks):
- |Aut| of cycles (directed and not), K₅, paths and stars;
- matches as |Aut L| times the occurrences, in complete graphs and on parallel edges;
- the dangling condition; relabelling with the trace of the vertices;
- one step of string rewriting (`grammar.hpp`, where it is present) and of graph rewriting on the words as paths giving the same forms, on 300 random forms of a monotone grammar for aⁿbⁿcⁿ;
- the derivations of chains from n labelled units, n!(n − 1)! in all, (n − 1)! for each, one orbit of Sₙ, equivariant at every derivable graph (n = 3, 4);
- twelve identical triangles with ten isolated vertices keeping one exact class under renumbering;
- a rule set closed under a reflection is equivariant, one not closed is not;
- rule classes ignoring names.

## 4. Alternatives Considered

- **Single-pushout rewriting** (Löwe 1993) deletes dangling edges implicitly. Chemistry wants every deleted bond named, so double pushouts with the dangling condition were kept.
- **Negative application conditions as graphs** (Habel, Heckel and Taentzer 1996) are planned. A predicate on the match serves the present conditions (two termini in different chains).

## 5. Consequences

The Chemistry Workbench's grammars run here (ADR 0051). A chemistry domain will state its reaction rules as `Rule`s and its symmetries as `Relabelling`s.
