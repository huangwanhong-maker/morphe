# arithmetic: Exact Arithmetic

*Status: draft (2026-10-24): integers, rationals, residues and the cyclotomic fields.* Other algebraic numbers, polynomial rings, factorisation and lattices of ℤⁿ are planned ([decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md), §2.1).

*New (2026-10-24):* the cyclotomic fields ℚ(ζ_m), exactly, in which the characters of every finite group take their values (`algebra/characters.hpp`).

## 1. Purpose

The numbers on which every exact structure of morphe is built: the integers ℤ without bound, the rationals ℚ, the residues of ℤ/n and the prime fields 𝔽_p, and the cyclotomic fields ℚ(ζ_m). Nothing here rounds and nothing overflows; a result that cannot be computed exactly is refused with an exception, never approximated.

## 2. Concepts

- **Integers in two forms.** An `Integer` lives in one machine word while its value fits a signed 64-bit integer, and as a sign and a magnitude in base 2³² otherwise. The form is canonical (a value that fits a word is always held in one), so equality compares parts. Arithmetic on words checks for overflow and moves to limbs only when it must, so the small entries of boundary matrices cost what machine integers cost.
- **The classical algorithms** (Knuth 1997, §4.3.1): schoolbook addition, subtraction and multiplication, and division by Algorithm D with the normalisation and the correction of the trial quotient written as in Warren 2012. Division truncates toward zero, as C++ does; `floor_div` and `mod` give the Euclidean pair.
- **Rationals in lowest terms** with a positive denominator, reduced by one gcd per operation.
- **Residues** of word size, multiplied in 128 bits; `is_prime` is deterministic for every 64-bit integer (Miller–Rabin with the first twelve primes as bases, exact below 3.3 × 10²⁴; Sorenson and Webster 2017).
- **Cyclotomic fields.** *Definition.* Φ_m = Π (x − ζ) over the primitive m-th roots of unity ζ, and ℚ(ζ_m) = ℚ[x]/(Φ_m). *Results used* (Washington 1997, ch. 2; Lang 2002, ch. VI, §3): Φ_m ∈ ℤ[x] is monic of degree φ(m) and irreducible over ℚ; x^m − 1 = Π_{d | m} Φ_d; the automorphisms of ℚ(ζ_m) are the σ_k: ζ ↦ ζ^k, k ∈ (ℤ/m)^×, complex conjugation being σ_{−1}. A `Cyclotomic` is held as its coordinates in the basis 1, ζ, …, ζ^{φ(m)−1}, unique by irreducibility, so equality in one field is equality of coordinates; numbers of two fields meet in the field of the least common multiple, ζ_m = ζ_M^{M/m}. Φ_m is computed as (x^m − 1) / Π_{d | m, d < m} Φ_d, every division exact by a monic divisor, and cached; the inverse comes from the extended Euclidean algorithm on (Φ_m, a), which ends on a non-zero constant because Φ_m is irreducible and deg a < φ(m). Rational numbers are always held in ℚ itself, so that a table of rational characters costs what rationals cost. Costs: O(φ) for a sum, O(φ²) for a product, O(m φ) for an automorphism, O(m²) once for Φ_m. `to_complex()` is for display only; nothing exact is decided from it.

## 3. Content

| Header | Content |
|---|---|
| `integer.hpp` | `Integer`: parsing and decimal text, `+ - * / %`, `divmod`, comparisons, `bit_length`, `to_double`; `abs`, `mod`, `floor_div`, `exact_div`, `gcd`, `extended_gcd`, `pow`, `factorial`, `double_factorial`, `binomial`; `std::hash<Integer>` |
| `rational.hpp` | `Rational`: arithmetic, comparisons, `parse`, `str`, `to_double` accurate for parts of any size; `abs`, `pow` |
| `modular.hpp` | `mod_add`, `mod_sub`, `mod_mul`, `mod_pow`, `mod_inverse`, `residue`, `is_prime` |
| `cyclotomic.hpp` | `euler_phi`, `cyclotomic_polynomial`; `Cyclotomic`: `zeta(m, k)`, `from_powers`, `from_basis`, `+ − × /`, `inverse`, `power`, `==`, `conj`, `galois(k)`, `lift(M)`, `coefficients_in(M)`, `field`, `is_rational`, `to_rational`, `to_complex` (display only), `str` |

## 4. Checks

`tests/check_cyclotomic.cpp`: Φ_m monic of degree φ(m) and Π_{d | m} Φ_d = x^m − 1 for m ≤ 120; Φ_p, Φ_{2m}(x) = Φ_m(−x), Φ_m(1); 105 the least order with a coefficient outside {−1, 0, 1}, that coefficient −2; the field axioms and inverses on random numbers of fourteen fields; ζ_m of order exactly m and Σ_k ζ_m^k = 0; each σ_k a ring homomorphism with σ_l σ_k = σ_{kl}, conjugation σ_{−1}, norms and traces rational; lifting a homomorphism, and numbers of different fields compared and combined; √2, √−3 and √5 found in ℚ(ζ_8), ℚ(ζ_3) and ℚ(ζ_5), and the Gauss sums of ten primes with g² = (−1)^{(p−1)/2} p (Ireland and Rosen 1990, ch. 6).

`tests/check_arithmetic.cpp`: 40,000 random operations against 128-bit machine arithmetic; division beyond 128 bits by the identity a = q b + r with |r| < |b| on operands built from all-ones and all-zero limbs, the hard cases of Algorithm D; 2¹⁰⁰, 25!, C(100, 50) and the Catalan numbers by two routes; gcd(2ᵃ − 1, 2ᵇ − 1) = 2^gcd(a, b) − 1 and Bézout's identity; telescoping and harmonic sums of rationals; π(10⁵) = 9592 and inverses modulo a prime.

## 5. Dependencies

None: the C++20 standard library and GCC's or Clang's overflow builtins and 128-bit integers.

## 6. References

- Knuth, D. E. (1997), *The Art of Computer Programming*, vol. 2, *Seminumerical Algorithms*, 3rd ed., Addison-Wesley, §4.3.1 and §4.5.
- Warren, H. S. (2012), *Hacker's Delight*, 2nd ed., Addison-Wesley, §9-2.
- Cohen, H. (1993), *A Course in Computational Algebraic Number Theory*, Springer, ch. 1.
- Sorenson, J. and Webster, J. (2017), Strong pseudoprimes to twelve prime bases, *Math. Comp.* 86, 985–1003.
- Washington, L. C. (1997), *Introduction to Cyclotomic Fields*, 2nd ed., Springer, ch. 2.
- Lang, S. (2002), *Algebra*, revised 3rd ed., Springer, ch. VI, §3.
- Ireland, K. and Rosen, M. (1990), *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, ch. 6.

## 7. Open Questions

1. **Faster multiplication** (Karatsuba, then FFT) when integers of thousands of digits appear.
2. **A borrowed back end** (GMP or FLINT, both LGPL) behind `Integer`, which needs a decision of its own (CLAUDE.md, open decision 4).
3. **Polynomials** over ℤ, ℚ and 𝔽_p, the first step toward `algebraic_geometry` (decision 0004, §2.3). `cyclotomic.hpp` keeps the few operations on polynomials over ℚ that it needs (product, division with remainder, reduction modulo Φ_m) in its `detail` namespace until a polynomial module is decided.
4. **Faster cyclotomic arithmetic.** Coordinates are rationals one by one; a common denominator, and the smallest field of each number (its conductor) rather than the field it was formed in, would make large fields cheaper. Numbers whose value is rational are already reduced to ℚ.
