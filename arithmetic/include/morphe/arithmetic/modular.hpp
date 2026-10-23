// Residues modulo a word-sized modulus: the rings ℤ/n and the prime fields 𝔽_p.
//
// Residues are unsigned 64-bit integers in [0, n); products are formed in 128 bits, so any
// modulus below 2^63 is allowed. is_prime is deterministic for every 64-bit integer: the
// Miller–Rabin test with the first twelve primes as bases is exact below 3.3 × 10^24
// (Sorenson and Webster 2017, Math. Comp. 86, 985–1003).

#pragma once

#include <cstdint>
#include <stdexcept>

#include "integer.hpp"

namespace morphe::arithmetic {

using Residue = std::uint64_t;

inline Residue mod_add(Residue a, Residue b, Residue n) { return a >= n - b ? a - (n - b) : a + b; }
inline Residue mod_sub(Residue a, Residue b, Residue n) { return a >= b ? a - b : a + (n - b); }
inline Residue mod_mul(Residue a, Residue b, Residue n) { return static_cast<Residue>(static_cast<unsigned __int128>(a) * b % n); }
inline Residue mod_pow(Residue a, std::uint64_t e, Residue n) {
  Residue r = 1 % n;
  a %= n;
  while (e) {
    if (e & 1) r = mod_mul(r, a, n);
    a = mod_mul(a, a, n);
    e >>= 1;
  }
  return r;
}
// The residue of an Integer modulo n.
inline Residue residue(const Integer& x, Residue n) {
  return static_cast<Residue>(mod(x, Integer(static_cast<unsigned long long>(n))).to_int64());
}
// The inverse of a modulo n; throws std::domain_error when gcd(a, n) != 1.
inline Residue mod_inverse(Residue a, Residue n) {
  const Bezout b = extended_gcd(Integer(static_cast<unsigned long long>(a % n)), Integer(static_cast<unsigned long long>(n)));
  if (!b.g.is_one()) throw std::domain_error("morphe::arithmetic::mod_inverse: not invertible");
  return residue(b.s, n);
}

inline bool is_prime(std::uint64_t n) {
  if (n < 2) return false;
  static constexpr std::uint64_t small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (auto p : small) {
    if (n % p == 0) return n == p;
  }
  std::uint64_t d = n - 1;
  int s = 0;
  while (!(d & 1)) { d >>= 1; s++; }
  for (auto a : small) {
    Residue x = mod_pow(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool composite = true;
    for (int r = 1; r < s; r++) {
      x = mod_mul(x, x, n);
      if (x == n - 1) { composite = false; break; }
    }
    if (composite) return false;
  }
  return true;
}

}  // namespace morphe::arithmetic
