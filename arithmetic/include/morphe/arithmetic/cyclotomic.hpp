// The cyclotomic fields ℚ(ζ_m), exactly.
//
// ℚ(ζ_m) = ℚ[x]/(Φ_m), where Φ_m is the m-th cyclotomic polynomial, the product of x − ζ over
// the primitive m-th roots of unity ζ. Φ_m has integer coefficients and degree φ(m), and it is
// irreducible over ℚ (Washington 1997, Introduction to Cyclotomic Fields, 2nd ed., ch. 2; Lang
// 2002, Algebra, rev. 3rd ed., ch. VI, §3), so 1, ζ, …, ζ^{φ(m)−1} is a basis of ℚ(ζ_m) over ℚ
// and every number has exactly one set of coordinates in it. A Cyclotomic is held as those
// coordinates, rationals, and equality of numbers held in the same field is equality of the
// coordinates.
//
// Φ_m is computed from x^m − 1 = Π_{d | m} Φ_d: x^m − 1 is divided by Φ_d for every proper
// divisor d of m, each division exact and by a monic polynomial, so over ℤ. The polynomials
// are kept in a cache shared by every number, guarded by a mutex; computing Φ_m costs
// O(m Σ_{d | m} φ(d)) = O(m²) operations on integers.
//
// Operations, with φ = φ(m):
//   - addition and subtraction coordinate by coordinate, O(φ);
//   - multiplication of the polynomials, then reduction modulo Φ_m from the top degree down,
//     O(φ²) plus O(φ · t) for the t non-zero coefficients of Φ_m;
//   - the inverse by the extended Euclidean algorithm in ℚ[x] on (Φ_m, a): it ends on a
//     non-zero constant c = s a + t Φ_m because Φ_m is irreducible and deg a < φ, so
//     a⁻¹ = s / c; O(φ²) operations on rationals whose size can grow, as in any Euclidean
//     algorithm over ℚ;
//   - the Galois automorphisms σ_k: ζ ↦ ζ^k for gcd(k, m) = 1, which are all the automorphisms
//     of ℚ(ζ_m) (the Galois group is (ℤ/m)^×, by the irreducibility of Φ_m), O(m φ); complex
//     conjugation is σ_{−1};
//   - lifting to ℚ(ζ_M) for m | M, through ζ_m = ζ_M^{M/m}, O(M φ(M)).
// Numbers of different fields are combined in the field of the least common multiple of
// their orders, to which both are lifted. A rational number is always held in ℚ = ℚ(ζ_1)
// (a number whose only non-zero coordinate is the constant one is rational, and conversely),
// so that rational values cost what rationals cost.
//
// to_complex() evaluates a number under the embedding ζ_m ↦ exp(2πi/m), in double precision,
// for display only: nothing exact is ever decided from it.

#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "integer.hpp"
#include "rational.hpp"

namespace morphe::arithmetic {

// Euler's totient φ(m) = |(ℤ/m)^×|, by trial division, O(√m).
inline std::uint64_t euler_phi(std::uint64_t m) {
  if (m == 0) throw std::invalid_argument("morphe::arithmetic::euler_phi: m = 0");
  std::uint64_t r = m, x = m;
  for (std::uint64_t q = 2; q * q <= x; q++)
    if (x % q == 0) {
      while (x % q == 0) x /= q;
      r -= r / q;
    }
  if (x > 1) r -= r / x;
  return r;
}

namespace detail {

// Polynomials over ℚ, the coefficient of x^i at index i, without trailing zeros.
using QPoly = std::vector<Rational>;

inline void trim_poly(QPoly& p) {
  while (!p.empty() && p.back().is_zero()) p.pop_back();
}

// P / D for a monic integer polynomial D that divides P; throws std::logic_error on a remainder.
inline std::vector<Integer> divide_monic_exactly(std::vector<Integer> P, const std::vector<Integer>& D) {
  const std::size_t dd = D.size() - 1;
  if (P.size() < D.size()) throw std::logic_error("morphe::arithmetic::cyclotomic: a divisor of higher degree");
  std::vector<Integer> Q(P.size() - dd, Integer(0));
  for (std::size_t t = P.size(); t-- > dd;) {
    const Integer c = P[t];
    if (c.is_zero()) continue;
    Q[t - dd] = c;
    for (std::size_t i = 0; i <= dd; i++) P[t - dd + i] -= c * D[i];
  }
  for (const auto& x : P)
    if (!x.is_zero()) throw std::logic_error("morphe::arithmetic::cyclotomic: an inexact division of x^m - 1");
  return Q;
}

// The field ℚ(ζ_m): its modulus Φ_m and the reduction modulo it.
struct CyclotomicField {
  unsigned m = 1;
  std::size_t phi = 1;
  std::vector<Integer> poly;  // Φ_m, low degree first, monic of degree φ(m)
  struct Term {
    std::size_t degree;
    Rational value;
    int unit;  // +1 or −1 when the value is a unit of ℤ, 0 otherwise
  };
  std::vector<Term> low;  // the non-zero coefficients of Φ_m below x^φ

  // p ← p mod Φ_m: degree below φ, trimmed. With Φ_m = x^φ + Σ a_i x^i, a term c x^t of degree
  // t ≥ φ is replaced by −c Σ a_i x^{t−φ+i}, from the top degree down.
  void reduce(QPoly& p) const {
    for (std::size_t t = p.size(); t-- > phi;) {
      if (p[t].is_zero()) continue;
      const Rational c = p[t];
      p[t] = Rational(0);
      for (const auto& term : low) {
        Rational& q = p[t - phi + term.degree];
        if (term.unit == 1) q -= c;
        else if (term.unit == -1) q += c;
        else q -= c * term.value;
      }
    }
    if (p.size() > phi) p.resize(phi);
    trim_poly(p);
  }
};

// The field ℚ(ζ_m), from a cache that lives as long as the program; Φ_m is computed on first use.
inline const CyclotomicField& cyclotomic_field(unsigned m) {
  if (m == 0) throw std::invalid_argument("morphe::arithmetic::cyclotomic: m = 0");
  static std::mutex mutex;
  static std::map<unsigned, std::unique_ptr<CyclotomicField>> cache;
  {
    std::lock_guard<std::mutex> lock(mutex);
    const auto it = cache.find(m);
    if (it != cache.end()) return *it->second;
  }
  // Φ_m = (x^m − 1) / Π_{d | m, d < m} Φ_d. The lock is not held while the divisors are computed.
  std::vector<Integer> P(static_cast<std::size_t>(m) + 1, Integer(0));
  P[0] = -1;
  P[m] = 1;
  for (unsigned d = 1; d < m; d++)
    if (m % d == 0) P = divide_monic_exactly(std::move(P), cyclotomic_field(d).poly);
  auto F = std::make_unique<CyclotomicField>();
  F->m = m;
  F->phi = P.size() - 1;
  if (F->phi != euler_phi(m)) throw std::logic_error("morphe::arithmetic::cyclotomic: deg Φ_m != φ(m)");
  for (std::size_t i = 0; i < F->phi; i++)
    if (!P[i].is_zero()) F->low.push_back({i, Rational(P[i]), P[i] == Integer(1) ? 1 : (P[i] == Integer(-1) ? -1 : 0)});
  F->poly = std::move(P);
  std::lock_guard<std::mutex> lock(mutex);
  const auto [it, inserted] = cache.emplace(m, std::move(F));
  return *it->second;
}

// a = q b + r in ℚ[x], deg r < deg b, for b non-zero.
inline void divmod_poly(QPoly a, const QPoly& b, QPoly& q, QPoly& r) {
  trim_poly(a);
  if (b.empty() || b.back().is_zero()) throw std::domain_error("morphe::arithmetic::cyclotomic: division by the zero polynomial");
  const std::size_t db = b.size() - 1;
  if (a.size() < b.size()) {
    q.clear();
    r = std::move(a);
    return;
  }
  q.assign(a.size() - db, Rational(0));
  const Rational lead = Rational(1) / b.back();
  for (std::size_t t = a.size(); t-- > db;) {
    if (a[t].is_zero()) continue;
    const Rational c = a[t] * lead;
    q[t - db] = c;
    for (std::size_t i = 0; i <= db; i++) a[t - db + i] -= c * b[i];
  }
  a.resize(db);
  trim_poly(a);
  trim_poly(q);
  r = std::move(a);
}

inline QPoly mul_poly(const QPoly& a, const QPoly& b) {
  if (a.empty() || b.empty()) return {};
  QPoly p(a.size() + b.size() - 1, Rational(0));
  for (std::size_t i = 0; i < a.size(); i++) {
    if (a[i].is_zero()) continue;
    for (std::size_t j = 0; j < b.size(); j++)
      if (!b[j].is_zero()) p[i + j] += a[i] * b[j];
  }
  trim_poly(p);
  return p;
}

inline QPoly sub_poly(const QPoly& a, const QPoly& b) {
  QPoly d(std::max(a.size(), b.size()), Rational(0));
  for (std::size_t i = 0; i < a.size(); i++) d[i] += a[i];
  for (std::size_t i = 0; i < b.size(); i++) d[i] -= b[i];
  trim_poly(d);
  return d;
}

// lcm(a, b), refused when it does not fit an unsigned.
inline unsigned common_order(unsigned a, unsigned b) {
  const std::uint64_t l = std::lcm(static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(b));
  if (l > 0x7fffffffull) throw std::overflow_error("morphe::arithmetic::Cyclotomic: the common field is too large");
  return static_cast<unsigned>(l);
}

}  // namespace detail

// The m-th cyclotomic polynomial Φ_m, coefficients low degree first; monic, of degree φ(m).
inline std::vector<Integer> cyclotomic_polynomial(unsigned m) { return detail::cyclotomic_field(m).poly; }

class Cyclotomic {
 public:
  Cyclotomic() = default;
  Cyclotomic(int v) : Cyclotomic(Rational(v)) {}
  Cyclotomic(long long v) : Cyclotomic(Rational(v)) {}
  Cyclotomic(const Integer& v) : Cyclotomic(Rational(v)) {}
  Cyclotomic(Rational q) {
    if (!q.is_zero()) c_.push_back(std::move(q));
  }

  // ζ_m^k, for any integer k.
  static Cyclotomic zeta(unsigned m, long long k = 1) {
    const long long mm = static_cast<long long>(m);
    const std::size_t e = static_cast<std::size_t>(((k % mm) + mm) % mm);
    std::vector<Rational> c(e + 1, Rational(0));
    c[e] = 1;
    return from_powers(m, c);
  }
  // Σ_k c[k] ζ_m^k, the exponents taken modulo m.
  static Cyclotomic from_powers(unsigned m, const std::vector<Rational>& c) {
    if (m == 1) {
      Rational s = 0;
      for (const auto& x : c) s += x;
      return Cyclotomic(s);
    }
    const auto& F = detail::cyclotomic_field(m);
    detail::QPoly p(std::min<std::size_t>(c.size(), m), Rational(0));
    for (std::size_t k = 0; k < c.size(); k++)
      if (!c[k].is_zero()) p[k % m] += c[k];
    F.reduce(p);
    return Cyclotomic(m, &F, std::move(p));
  }
  // The number with coordinates c in the basis 1, ζ_m, …, ζ_m^{φ(m)−1}.
  static Cyclotomic from_basis(unsigned m, std::vector<Rational> c) {
    if (c.size() > euler_phi(m)) throw std::invalid_argument("morphe::arithmetic::Cyclotomic::from_basis: more than φ(m) coordinates");
    if (m == 1) return Cyclotomic(c.empty() ? Rational(0) : c[0]);
    return Cyclotomic(m, &detail::cyclotomic_field(m), std::move(c));
  }

  // The order m of the field ℚ(ζ_m) in which the number is held (1 for a rational).
  unsigned field() const { return m_; }
  // Its coordinates in the power basis of that field, without trailing zeros.
  const std::vector<Rational>& coefficients() const { return c_; }
  bool is_zero() const { return c_.empty(); }
  bool is_rational() const { return m_ == 1; }
  Rational to_rational() const {
    if (m_ != 1) throw std::domain_error("morphe::arithmetic::Cyclotomic: " + str() + " is not rational");
    return c_.empty() ? Rational(0) : c_[0];
  }

  // The same number held in ℚ(ζ_M), for m | M (a rational stays in ℚ).
  Cyclotomic lift(unsigned M) const {
    if (M == 0 || M % m_ != 0) throw std::invalid_argument("morphe::arithmetic::Cyclotomic::lift: the order of the field must divide M");
    if (m_ == 1 || M == m_) return *this;
    const auto& F = detail::cyclotomic_field(M);
    const std::size_t t = M / m_;
    detail::QPoly p((c_.size() - 1) * t + 1, Rational(0));
    for (std::size_t i = 0; i < c_.size(); i++) p[i * t] = c_[i];
    F.reduce(p);
    return Cyclotomic(M, &F, std::move(p));
  }
  // The coordinates in the power basis of ℚ(ζ_M), m | M, all φ(M) of them: the canonical form there.
  std::vector<Rational> coefficients_in(unsigned M) const {
    std::vector<Rational> c = lift(M).c_;
    c.resize(euler_phi(M), Rational(0));
    return c;
  }

  // σ_k: ζ_m ↦ ζ_m^k, gcd(k, m) = 1. It is the restriction of σ_k of every ℚ(ζ_M) ⊇ ℚ(ζ_m).
  Cyclotomic galois(long long k) const {
    if (m_ == 1) return *this;
    const long long mm = static_cast<long long>(m_);
    const std::uint64_t kk = static_cast<std::uint64_t>(((k % mm) + mm) % mm);
    if (std::gcd(kk, static_cast<std::uint64_t>(m_)) != 1)
      throw std::invalid_argument("morphe::arithmetic::Cyclotomic::galois: k is not prime to the order of the field");
    detail::QPoly p(m_, Rational(0));
    for (std::size_t i = 0; i < c_.size(); i++) p[(i * kk) % m_] += c_[i];
    f_->reduce(p);
    return Cyclotomic(m_, f_, std::move(p));
  }
  // Complex conjugation, ζ ↦ ζ⁻¹.
  Cyclotomic conj() const { return galois(-1); }

  Cyclotomic inverse() const {
    if (is_zero()) throw std::domain_error("morphe::arithmetic::Cyclotomic: the inverse of zero");
    if (m_ == 1) return Cyclotomic(Rational(1) / c_[0]);
    // The extended Euclidean algorithm on (Φ_m, a), keeping r_i ≡ s_i a (mod Φ_m).
    detail::QPoly r0, r1 = c_, s0, s1{Rational(1)}, q, r;
    for (const auto& x : f_->poly) r0.push_back(Rational(x));
    while (!r1.empty()) {
      detail::divmod_poly(r0, r1, q, r);
      detail::QPoly s2 = detail::sub_poly(s0, detail::mul_poly(q, s1));
      r0 = std::move(r1);
      r1 = std::move(r);
      s0 = std::move(s1);
      s1 = std::move(s2);
    }
    if (r0.size() != 1) throw std::logic_error("morphe::arithmetic::Cyclotomic: Φ_m has a common factor with a number");
    const Rational c = Rational(1) / r0[0];
    for (auto& x : s0) x *= c;
    f_->reduce(s0);
    return Cyclotomic(m_, f_, std::move(s0));
  }
  // a^e for any integer e (a non-zero when e < 0).
  Cyclotomic power(long long e) const {
    Cyclotomic base = e < 0 ? inverse() : *this, r = 1;
    unsigned long long k = e < 0 ? static_cast<unsigned long long>(-(e + 1)) + 1 : static_cast<unsigned long long>(e);
    while (k) {
      if (k & 1) r *= base;
      k >>= 1;
      if (k) base *= base;
    }
    return r;
  }

  // For display only: the value under ζ_m ↦ exp(2πi/m), in double precision.
  std::complex<double> to_complex() const {
    std::complex<double> z = 0;
    for (std::size_t i = 0; i < c_.size(); i++)
      z += c_[i].to_double() * std::polar(1.0, 2 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(m_));
    return z;
  }
  // "1 - 2*ζ12^3 + 1/2*ζ12": the coordinates in the power basis of the field it is held in.
  std::string str() const {
    if (c_.empty()) return "0";
    std::string s;
    for (std::size_t i = 0; i < c_.size(); i++) {
      if (c_[i].is_zero()) continue;
      const bool neg = c_[i].sign() < 0;
      const Rational a = neg ? -c_[i] : c_[i];
      std::string term;
      if (i == 0) term = a.str();
      else {
        term = (a == Rational(1) ? std::string() : a.str() + "*") + "ζ" + std::to_string(m_);
        if (i > 1) term += "^" + std::to_string(i);
      }
      if (s.empty()) s = (neg ? "-" : "") + term;
      else s += (neg ? " - " : " + ") + term;
    }
    return s;
  }

  Cyclotomic operator-() const {
    Cyclotomic x = *this;
    for (auto& c : x.c_) c = -c;
    return x;
  }
  friend Cyclotomic operator+(const Cyclotomic& a, const Cyclotomic& b) {
    if (a.m_ == b.m_ || a.m_ == 1 || b.m_ == 1) {
      const Cyclotomic& host = a.m_ == 1 ? b : a;
      detail::QPoly s(std::max(a.c_.size(), b.c_.size()), Rational(0));
      for (std::size_t i = 0; i < a.c_.size(); i++) s[i] += a.c_[i];
      for (std::size_t i = 0; i < b.c_.size(); i++) s[i] += b.c_[i];
      return Cyclotomic(host.m_, host.f_, std::move(s));
    }
    const unsigned M = detail::common_order(a.m_, b.m_);
    return a.lift(M) + b.lift(M);
  }
  friend Cyclotomic operator-(const Cyclotomic& a, const Cyclotomic& b) { return a + (-b); }
  friend Cyclotomic operator*(const Cyclotomic& a, const Cyclotomic& b) {
    if (a.is_zero() || b.is_zero()) return Cyclotomic();
    if (a.m_ == 1 || b.m_ == 1) {
      const Cyclotomic& x = a.m_ == 1 ? b : a;
      const Rational& q = (a.m_ == 1 ? a : b).c_[0];
      detail::QPoly p = x.c_;
      for (auto& c : p) c *= q;
      return Cyclotomic(x.m_, x.f_, std::move(p));
    }
    if (a.m_ != b.m_) {
      const unsigned M = detail::common_order(a.m_, b.m_);
      return a.lift(M) * b.lift(M);
    }
    detail::QPoly p = detail::mul_poly(a.c_, b.c_);
    a.f_->reduce(p);
    return Cyclotomic(a.m_, a.f_, std::move(p));
  }
  friend Cyclotomic operator/(const Cyclotomic& a, const Cyclotomic& b) { return a * b.inverse(); }
  Cyclotomic& operator+=(const Cyclotomic& b) { return *this = *this + b; }
  Cyclotomic& operator-=(const Cyclotomic& b) { return *this = *this - b; }
  Cyclotomic& operator*=(const Cyclotomic& b) { return *this = *this * b; }
  Cyclotomic& operator/=(const Cyclotomic& b) { return *this = *this / b; }

  // Equality of numbers: of the coordinates when both are held in one field, and otherwise in
  // the field of the least common multiple (a rational is never equal to an irrational).
  friend bool operator==(const Cyclotomic& a, const Cyclotomic& b) {
    if (a.m_ == b.m_) return a.c_ == b.c_;
    if (a.m_ == 1 || b.m_ == 1) return false;
    const unsigned M = detail::common_order(a.m_, b.m_);
    return a.lift(M).c_ == b.lift(M).c_;
  }

 private:
  // From coordinates in ℚ(ζ_m); trims, and holds a rational result in ℚ.
  Cyclotomic(unsigned m, const detail::CyclotomicField* f, detail::QPoly c) : m_(m), f_(f), c_(std::move(c)) {
    detail::trim_poly(c_);
    if (c_.size() <= 1) m_ = 1, f_ = nullptr;
  }

  unsigned m_ = 1;
  const detail::CyclotomicField* f_ = nullptr;  // the field ℚ(ζ_m) from the cache; null for ℚ
  std::vector<Rational> c_;
};

}  // namespace morphe::arithmetic
