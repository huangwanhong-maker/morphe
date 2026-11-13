# 0015. Acyclic Relations, and Forward Chaining

- **Status:** Accepted, provisionally (2026-11-13).
- **Decided by:** Claude, for the infrastructure's ADR 0110, in which the history domain reads the dependence of historical objects as an order; the user's rule that a general method is written in the foundation where it belongs.

## 1. Context

The infrastructure's history domain gained records in which one object's appearing rests on others. What is asked of such a record is order theory and a little logic: whether the relation has a cycle, and where; what lies before and after an element; the least relation with the same closure; the relation a map to classes induces; and, since an appearing has several premises together, what follows from a set of givens and what rests on one of them. None of it names history. morphe held Möbius functions of posets given whole, and graphs without direction; it had none of these.

## 2. Decision

1. **`combinatorics/acyclic.hpp`** holds, for a relation on {0, …, n−1} given as arcs:
   - a topological order by Kahn's algorithm, or a cycle as the place where there is none;
   - what is reached from an element forward or backward, with the least number of arcs, to a bound if one is given;
   - reachability, and whether two elements are comparable;
   - the covering relation of an acyclic relation (Aho, Garey and Ullman 1972);
   - the quotient by a map to classes, each arc between classes with the number of arcs it gathers.
2. **Forward chaining over rules** with several premises and one conclusion, in the same header: the least set holding the givens and closed under the rules, in time linear in their size (Dowling and Gallier 1984), with elements that may be withheld; and what rests on an element, those had with it and not without it.
3. **Indices, not names.** Elements are indices, as in the rest of `combinatorics`; a domain keeps its own names.
4. **A cycle is an answer.** A relation with a cycle is not an error: the cycle is returned. Only the covering relation, defined for acyclic relations, refuses one.

## 3. Alternatives Considered

- **In the history domain.** Dependence among packages, tasks, proofs and reactions is the same mathematics; the philosophy domain's "raised by" already needs the cycle.
- **On `Multigraph`.** It has no direction, and its incidence serves spanning forests. Arcs as pairs are the least that states a relation.
- **Horn clauses in `logic`.** That module is planned and will hold formulae; forward chaining on sets of indices is combinatorics, and can be reached from there when it is written.

## 4. Consequences

- 18 checks (`combinatorics/tests/check_acyclic.cpp`): orders and cycles on 200 random relations; cones against Warshall's closure; the covering relation with the relation's closure and no arc to spare; quotients' counts; forward chaining and what rests on an element against iteration to a fixed point on 300 random sets of rules.
- The history domain's `objects` and `dependence` stand on it.
- **Not built:** strong components (tyche's `networks` holds them, waiting to move here), lattices of down-sets, linear extensions counted.
