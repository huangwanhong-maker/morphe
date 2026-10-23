# arithmetic: Exact Arithmetic

*Status: draft (2026-10-23): integers, rationals and residues.* Algebraic numbers, polynomial rings, factorisation and lattices of ℤⁿ are planned ([decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md), §2.1).

## 1. Purpose

The numbers on which every exact structure of morphe is built: the integers ℤ without bound, the rationals ℚ, and the residues of ℤ/n and the prime fields 𝔽_p. Nothing here rounds and nothing overflows; a result that cannot be computed exactly is refused with an exception, never approximated.

## 2. Concepts

- **Integers in two forms.** An `Integer` lives in one machine word while its value fits a signed 64-bit integer, and as a sign and a magnitude in base 2³² otherwise. The form is canonical (a value that fits a word is always held in one), so equality compares parts. Arithmetic on words checks for overflow and moves to limbs only when it must, so the small entries of boundary matrices cost what machine integers cost.
- **The classical algorithms** (Knuth 1997, §4.3.1): schoolbook addition, subtraction and multiplication, and division by Algorithm D with the normalisation and the correction of the trial quotient written as in Warren 2012. Division truncates toward zero, as C++ does; `floor_div` and `mod` give the Euclidean pair.
- **Rationals in lowest terms** with a positive denominator, reduced by one gcd per operation.
- **Residues** of word size, multiplied in 128 bits; `is_prime` is deterministic for every 64-bit integer (Miller–Rabin with the first twelve primes as bases, exact below 3.3 × 10²⁴; Sorenson and Webster 2017).

## 3. Content

| Header | Content |
|---|---|
| `integer.hpp` | `Integer`: parsing and decimal text, `+ - * / %`, `divmod`, comparisons, `bit_length`, `to_double`; `abs`, `mod`, `floor_div`, `exact_div`, `gcd`, `extended_gcd`, `pow`, `factorial`, `double_factorial`, `binomial`; `std::hash<Integer>` |
| `rational.hpp` | `Rational`: arithmetic, comparisons, `parse`, `str`, `to_double` accurate for parts of any size; `abs`, `pow` |
| `modular.hpp` | `mod_add`, `mod_sub`, `mod_mul`, `mod_pow`, `mod_inverse`, `residue`, `is_prime` |

## 4. Checks

`tests/check_arithmetic.cpp`: 40,000 random operations against 128-bit machine arithmetic; division beyond 128 bits by the identity a = q b + r with |r| < |b| on operands built from all-ones and all-zero limbs, the hard cases of Algorithm D; 2¹⁰⁰, 25!, C(100, 50) and the Catalan numbers by two routes; gcd(2ᵃ − 1, 2ᵇ − 1) = 2^gcd(a, b) − 1 and Bézout's identity; telescoping and harmonic sums of rationals; π(10⁵) = 9592 and inverses modulo a prime.

## 5. Dependencies

None: the C++20 standard library and GCC's or Clang's overflow builtins and 128-bit integers.

## 6. References

- Knuth, D. E. (1997), *The Art of Computer Programming*, vol. 2, *Seminumerical Algorithms*, 3rd ed., Addison-Wesley, §4.3.1 and §4.5.
- Warren, H. S. (2012), *Hacker's Delight*, 2nd ed., Addison-Wesley, §9-2.
- Cohen, H. (1993), *A Course in Computational Algebraic Number Theory*, Springer, ch. 1.
- Sorenson, J. and Webster, J. (2017), Strong pseudoprimes to twelve prime bases, *Math. Comp.* 86, 985–1003.

## 7. Open Questions

1. **Faster multiplication** (Karatsuba, then FFT) when integers of thousands of digits appear.
2. **A borrowed back end** (GMP or FLINT, both LGPL) behind `Integer`, which needs a decision of its own (CLAUDE.md, open decision 4).
3. **Polynomials** over ℤ, ℚ and 𝔽_p, the first step toward `algebraic_geometry` (decision 0004, §2.3).
