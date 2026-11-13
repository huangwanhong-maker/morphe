# 0014. Braid Groups, and Laurent Polynomials

- **Status:** Accepted, provisionally (2026-11-13).
- **Decided by:** the user, for whom the braid group is "sehr important" and who asked that it be built early, on their own draft on trajectories and braids; and Claude. The infrastructure's ADR 0103.

## 1. Context

The infrastructure founded kinesis, a foundation for dynamical systems, with the reading of braids from records of trajectories among its first modules. A braid read from a record is a word in Artin's generators, and what is asked of it is algebra: whether two words are the same braid, what a braid does to the curves of a punctured disc, what its representation and the polynomial of its closure are. None of that names a trajectory.

morphe holds groups given by generators and relations and finite groups, and had no braid groups. Its presented groups are simplified by Tietze moves, which decide nothing in general; the braid groups have a word problem that is solved, by several methods.

## 2. Decision

1. **The braid groups are morphe's,** in `algebra/braid.hpp`: words in Artin's generators with their permutation, exponent sum and crossings; the word problem by Dehornoy's handle reduction, with Dehornoy's order; Artin's faithful action on the free group, as the independent check of the first; Dynnikov's coordinates of curve systems, on integers without bound; the unreduced Burau representation over ℤ[t, t⁻¹], with the Alexander polynomial of a closure as a first minor of I − ψ(β).
2. **Laurent polynomials are `arithmetic`'s,** in `arithmetic/laurent.hpp`: ℤ[t, t⁻¹] in one variable, exact, with the representative of a polynomial up to units.
3. **Conventions, fixed once.** A word is read from the left in time. σ_i is the counterclockwise half-exchange of the strands at positions i and i + 1, the left one passing below, so that a full counterclockwise turn of two strands is σ_1² and the exponent sum of a pure braid on two strands is twice its winding number. Sources that take σ_i clockwise (Thiffeault 2010; braidlab) are read with every letter inverted: Dynnikov's rules as Thiffeault states them are for his σ_i, which is this library's σ_i⁻¹.
4. **What is read from trajectories is kinesis's:** the braid of a record, what a record certifies under a bound on speed, the windings a record admits, and the entropy of a braid as a statement about a flow. kinesis uses this header.
5. **Bounds are reported.** Handle reduction ends for every word, by Dehornoy's theorem, and no polynomial bound on its steps is proved: a bound on steps is taken, and passing it is reported as such, never answered as an inequality. Artin's images can grow exponentially, and are computed within a bound on length.

## 3. Alternatives Considered

- **In kinesis, with the trajectories.** The word problem, the coordinates and the representation serve knots, mapping classes and configuration spaces as well; they are algebra and belong with the groups.
- **Through `presentation.hpp`.** A braid group is a presented group, and Tietze simplification does not solve its word problem. The braid groups have their own structure, and a type of their own.
- **Garside's normal form first.** It solves the word and conjugacy problems and gives a canonical word. Handle reduction is shorter to state, is fast in practice, and needed nothing else beneath it; Garside's form is planned, for conjugacy.
- **The Burau representation to decide equality.** It is faithful for three strands and not for five or more, and its faithfulness for four is open in the sources read. It is given as a representation and never used as a test of equality.
- **Polynomials in general.** `arithmetic` plans polynomial rings. One variable with negative exponents is what a braid needs today, and it is written so that a general ring can take its place.

## 4. Consequences

- 39 checks (`algebra/tests/check_braid.cpp`): the relations by handle reduction; handle reduction against Artin's action on words made trivial by inserted relators and on random words, for three to six strands; Garside's element; Dehornoy's order on a sample; Dynnikov's rules by the braid relations on random coordinates, and the growth rates of σ_1 σ_2⁻¹ and σ_1 σ_2 σ_3⁻¹ against log((3 + √5)/2) and the logarithm of the largest root of x⁴ − 2x³ − 2x + 1; the Burau matrices, and the Alexander polynomials of the unknot, the Hopf link, the trefoil, the figure-eight knot and 5_1.
- Not written: Garside's normal form, conjugacy, the Lawrence–Krammer representation, braids on surfaces, and monoids in which strands are created and removed.
- The planned order of morphe's modules is unchanged.
