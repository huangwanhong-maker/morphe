# 0002. Scope, Boundaries and Independence

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user, who proposed algebraic structures as a project of its own, for investigating topology, algebraic geometry, category theory, formal structures and equivalences such as bisimulation and bimorphism. Claude drew the boundaries ([laboratory ADR 0008](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/decisions/0008-independent-libraries.md)).

## 1. Context

The laboratory first planned a `mathematics/` layer inside its own core, covering algebra, geometry and topology. Mathematical structures are broader than their uses in physics, however. They are also objects of research for the user, in formal language theory and in group-theoretic approaches to sequences, among other areas.

## 2. Decision

1. **Scope.** morphe covers:
   - exact arithmetic;
   - algebra (including homological algebra and representations);
   - category theory;
   - topology;
   - algebraic geometry;
   - formal structures;
   - notions of equivalence.
2. **Excluded.**
   - Geometry as physics uses it belongs to physis. This includes differential geometry, rotations, the sphere and ellipsoid, and map projections.
   - Floating-point numerical methods belong to calx.
3. **Dependencies.** morphe depends on no other project library. physis depends on morphe.
4. **Exactness.** Computation is exact by default.
5. **Language.** morphe is written in C++, following laboratory ADR 0003.
6. **Organising principle (provisional).** Category theory forms the spine of the library. Algebras describe construction, and coalgebras describe behaviour.

## 3. Alternatives Considered

- **Mathematics inside physis.** Rejected. Structures are broader than physics, and research in formal languages and genomics can reuse them directly.
- **Depending on calx for fast linear algebra.** Rejected, to keep morphe light. Computations at GPU scale can go through optional adapters, each decided separately.
- **Including differential geometry.** Rejected. In physics it is numerical and tied to coordinates, so it belongs with physis. Algebraic geometry and topology stay here.

## 4. Consequences

- Exact sparse matrices over ℤ and over finite fields are implemented here. They are kept separate from calx's floating-point sparse matrices.
- physis builds structure-preserving discretisations on morphe's complexes.
