# 0010. Surfaces Given as Lists of Triangles: Their Complexes Read Directly, Their Components, and Their Orientation

- **Status:** Accepted, provisionally
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked for a Body Atlas, the anatomy of the whole body from BodyParts3D's meshes.
  - **The rule.** General methods are written in the foundations.
  - **What Claude did.** It wrote the header, its check and this record.
  - **The composition.** The infrastructure's [ADR 0066](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0066-the-anatomy-domain-and-the-body-atlas.md). Its anatomy domain tells, for each part of the body, whether its mesh is a closed surface, of what genus, and whether its triangles are consistently oriented, which it must be to bound a solid.

## 1. Context

morphe classifies surfaces exactly (`topology/surfaces.hpp`), on a CW complex. A list of triangles became one through `SimplicialComplex::from_facets(...).cw()`, which closes the list under subsets: it stores every face of every triangle in ordered sets. For the 6.7 million triangles of a whole body, some meshes with 200,000 triangles, that costs far more than the classification itself, which is linear.

That route also discards something the meshes carry. A simplicial complex holds each triangle by its sorted vertices, so it forgets the orientation the triangle's corners give. Yet whether that orientation is consistent is exactly what decides whether the surface bounds a solid by the divergence theorem.

`classify_surface` asks for a connected complex, and a mesh may have several components.

## 2. Decision

`topology/triangulated.hpp`:
- **`triangle_complex`.** The 2-complex of a list of triangles, built directly:
  - its edges are the corner pairs that occur, each once, from the lesser vertex to the greater, found by sorting the 3F pairs in O(F log F);
  - each triangle is attached along its corners in their own order, a → b → c → a.
- **`triangle_components`.** The triangles grouped by the components of the graph of their edges (`combinatorics::spanning_forest`).
- **`classify_triangles`.** Each component renumbered on its own vertices and classified by `classify_surface`. A component is *consistent* when it is an orientable surface and the orientation `classify_surface` returns gives every cell the same sign. That holds exactly when the triangles' own orientation traverses every interior edge once each way.

## 3. Alternatives Considered

- **Counting V, E, F and boundary edges in the domain.** Rejected: it would write the classification of surfaces again outside morphe. `classify_surface` also checks the links of vertices, which counts do not.
- **Making `from_facets` faster.** Rejected for now. The closure under subsets is what a simplicial complex is, and other callers need its faces. A surface given as triangles is a different input, and deserves a direct reading.

## 4. Consequences

- **Checks.** morphe gains `check_triangulated`, with 10 checks:
  - a tetrahedron's boundary, and the same with one triangle turned over: an orientable sphere whose given orientation is inconsistent;
  - the 3 × 3 torus;
  - a Möbius strip;
  - a fan, which is a disc;
  - two tetrahedra apart, and two sharing a vertex (a pinch point);
  - on every example, agreement with `from_facets(...).cw()` in cells, Euler characteristic and classification.
- **The anatomy domain** classifies BodyParts3D's 2,234 meshes with it in seconds. Of them, 2,201 are closed surfaces, consistently oriented.
