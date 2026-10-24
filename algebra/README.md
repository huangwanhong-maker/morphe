# algebra: Exact Linear Algebra, Homology, Presented and Finite Groups

*Status: draft (2026-10-23): integer matrices and the Smith normal form, chain complexes and their homology, finitely presented groups, finite groups.* Hermite normal forms, permutation groups by Schreier–Sims, character tables, coset enumeration, rewriting, Lie algebras over exact rings and universal algebra are planned ([decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md), §2.1).

## 1. Purpose

The algebra that topology and physics need exactly: linear algebra over ℤ, the homology of chain complexes with its generators, groups given by generators and relations, and finite groups with the homomorphisms into them. Its numerical counterpart, linear algebra in floating point, is calx's `linalg`; Lie groups acting on space are physis's `symmetry`.

## 2. Concepts

- **The Smith normal form** U A V = diag(d₁ | d₂ | … | d_r, 0, …) with U, V unimodular (Smith 1861; Newman 1972, Theorem II.9). The invariant factors give the abelian group a matrix presents, and the rank over ℚ. The algorithm brings a least entry to the pivot, clears its row and column by division with remainder, and restores divisibility (Cohen 1993, Algorithm 2.4.14); the transforms and their inverses are accumulated on request.
- **Homology in two phases.** First an algebraic reduction of the chain complex by unit pivots (Kaczynski, Mischaikow and Mrozek 2004, §4.4), which for a cell complex contracts a spanning tree and removes free faces, keeping for each surviving cell its lift to an original chain; then the Smith normal form of the small remainder, from whose transforms the generators are read (Munkres 1984, §11). Generators are cycles of the original complex: a basis of the free part, and a cycle of the stated order for each torsion factor. Cohomology follows by the universal coefficient theorem, Betti numbers over 𝔽_p from ranks modulo p.
- **Presentations and Tietze moves.** Words are sequences of letters ±i; `simplify` reduces relators freely and cyclically, removes repeats up to rotation and inversion, and eliminates a generator occurring once in a relator (Lyndon and Schupp 1977, §II.2). It terminates and preserves the group, and claims no normal form: the isomorphism problem is undecidable (Adian 1957; Rabin 1958).
- **Finite groups and homomorphisms.** Groups by their tables, built for the classical families and by closure of permutations; homomorphisms from a presented group are enumerated by backtracking, testing each relator as soon as its generators are assigned. For the fundamental group of a surface the count is a character sum (Frobenius 1896; Mednykh 1978), which the checks use.

## 3. Content

| Header | Content |
|---|---|
| `matrix.hpp` | `Matrix` (dense, row-major), `SparseMatrix` (by columns), `Vector` of `Entry{index, value}`, `axpy`, `coefficient`, products |
| `smith.hpp` | `smith(A, {U, V, U_inverse, V_inverse})` → `Smith{D, diagonal, U, V, U_inverse, V_inverse}`; `torsion`; `rank_mod(A, p)` |
| `chain_complex.hpp` | `ChainComplex{n, boundary}` with `validate`, `euler_characteristic`; `reduce`; `homology(C, generators)` → `HomologyGroup{rank, torsion, free_generators, torsion_generators}`; `betti(C, p)`; `cohomology(H)` |
| `presentation.hpp` | `Word`, `inverse`, `reduce`, `cyclic_reduce`, `canonical_cyclic`; `Presentation{generators, relators, names}`; `simplify` → `Simplified{presentation, kept}`; `abelianisation` → `AbelianGroup{rank, torsion}` |
| `determinant.hpp` | `determinant(A)` over ℤ by Bareiss's fraction-free elimination with row exchanges (Bareiss 1968), every intermediate a minor; `minor_matrix`, `minor`; `determinant(n, entries)` over ℚ by clearing the denominators of each row |
| `finite_group.hpp` | `FiniteGroup` (table, inverses, `evaluate`, `element_order`, `is_abelian`); `cyclic`, `dihedral`, `symmetric`, `alternating`, `quaternion`, `direct_product`, `permutation_group`; `conjugacy_classes`; `for_each_homomorphism`, `count_homomorphisms` |

## 4. Checks

`tests/check_algebra.cpp`: a textbook Smith form; U A V = D, unimodularity and the divisibility chain on 300 random matrices, d₁ and d₁d₂ against gcds of minors; ranks modulo primes; the homology of 40 complexes built from known elementary pieces and disguised by random unimodular changes of basis, with generators that are cycles, free generators independent modulo boundaries, and Betti numbers mod 2 by universal coefficients; ℝP² with its torsion generator and its cohomology; Tietze simplification and abelianisation (surface groups, the Klein bottle); orders and class numbers of seven groups up to A₅; homomorphisms from the groups of the torus, the genus-2 surface, the projective plane and the Klein bottle into them against the formulas of Frobenius and Mednykh, with Frobenius–Schur indicators for the non-orientable ones.

`tests/check_determinant.cpp`: Bareiss's elimination against the Leibniz expansion on random matrices with many zeros (so that pivots are exchanged); the sign of every permutation of five letters; multiplicativity on matrices whose determinants exceed 64 bits; Vandermonde's product; singular matrices; the Hilbert matrices over ℚ against det Hₙ = cₙ⁴ / c₂ₙ, cₙ = ∏_{i<n} i!; minors against the expansion along a row.

## 5. Dependencies

morphe `arithmetic`.

## 6. References

- Smith, H. J. S. (1861), On systems of linear indeterminate equations and congruences, *Phil. Trans. R. Soc.* 151, 293–326.
- Newman, M. (1972), *Integral Matrices*, Academic Press.
- Cohen, H. (1993), *A Course in Computational Algebraic Number Theory*, Springer, §2.4.
- Kaczynski, T., Mischaikow, K. and Mrozek, M. (2004), *Computational Homology*, Springer, ch. 3–4.
- Munkres, J. R. (1984), *Elements of Algebraic Topology*, Addison-Wesley, §11.
- Dumas, J.-G., Saunders, B. D. and Villard, G. (2001), On efficient sparse integer matrix Smith normal form computations, *J. Symb. Comput.* 32, 71–99.
- Lyndon, R. C. and Schupp, P. E. (1977), *Combinatorial Group Theory*, Springer.
- Mednykh, A. D. (1978), Determination of the number of nonequivalent coverings over a compact Riemann surface, *Soviet Math. Dokl.* 19, 318–320.
- James, G. and Liebeck, M. (2001), *Representations and Characters of Groups*, 2nd ed., Cambridge University Press.

## 7. Open Questions

1. **Large residues.** The dense phase is cubic; complexes whose reduction leaves thousands of cells would need modular Smith forms (Dumas, Saunders and Villard 2001).
2. **Characters of finite groups** computed rather than quoted (Dixon–Schneider), so that the formula of Frobenius and Mednykh and the gauge theories of finite groups can use any group.
3. **Coset enumeration and rewriting** (Todd–Coxeter, Knuth–Bendix) for presented groups that are finite or automatic.
