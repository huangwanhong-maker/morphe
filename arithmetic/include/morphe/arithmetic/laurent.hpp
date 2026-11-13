// Laurent polynomials ℤ[t, t⁻¹] in one variable, exact.
//
// A polynomial is a finite sum Σ c_k tᵏ over k ∈ ℤ with coefficients in ℤ (integers without bound), kept as the
// map from the exponents with a coefficient that is not zero. The ring is commutative with units ±tᵏ. normalised()
// picks the representative of a polynomial up to units: its lowest exponent 0 and its lowest coefficient positive,
// which is how an Alexander polynomial is written.
#pragma once

#include <morphe/arithmetic/integer.hpp>

#include <complex>
#include <map>
#include <string>
#include <utility>

namespace morphe::arithmetic {

class Laurent {
 public:
  Laurent() = default;
  Laurent(const Integer& c) {  // a constant
    if (!c.is_zero()) terms_[0] = c;
  }
  Laurent(long long c) : Laurent(Integer(c)) {}
  // c tᵏ
  static Laurent monomial(const Integer& c, long k) {
    Laurent p;
    if (!c.is_zero()) p.terms_[k] = c;
    return p;
  }
  static Laurent t(long k = 1) { return monomial(Integer(1), k); }

  bool is_zero() const { return terms_.empty(); }
  const std::map<long, Integer>& terms() const { return terms_; }
  Integer coefficient(long k) const {
    const auto it = terms_.find(k);
    return it == terms_.end() ? Integer(0) : it->second;
  }
  long lowest() const { return terms_.empty() ? 0 : terms_.begin()->first; }
  long highest() const { return terms_.empty() ? 0 : terms_.rbegin()->first; }

  Laurent operator-() const {
    Laurent p;
    for (const auto& [k, c] : terms_) p.terms_[k] = -c;
    return p;
  }
  friend Laurent operator+(Laurent a, const Laurent& b) {
    for (const auto& [k, c] : b.terms_) a.add(k, c);
    return a;
  }
  friend Laurent operator-(Laurent a, const Laurent& b) {
    for (const auto& [k, c] : b.terms_) a.add(k, -c);
    return a;
  }
  friend Laurent operator*(const Laurent& a, const Laurent& b) {
    Laurent p;
    for (const auto& [k, c] : a.terms_)
      for (const auto& [l, d] : b.terms_) p.add(k + l, c * d);
    return p;
  }
  friend bool operator==(const Laurent& a, const Laurent& b) { return a.terms_ == b.terms_; }

  // The polynomial up to units: divided by ±t^lowest so that it begins with a positive constant.
  Laurent normalised() const {
    if (is_zero()) return *this;
    const long shift = lowest();
    const bool negative = terms_.begin()->second.sign() < 0;
    Laurent p;
    for (const auto& [k, c] : terms_) p.terms_[k - shift] = negative ? -c : c;
    return p;
  }

  // Its value at a complex number that is not zero.
  std::complex<double> at(std::complex<double> z) const {
    std::complex<double> v = 0;
    for (const auto& [k, c] : terms_) v += c.to_double() * std::pow(z, static_cast<double>(k));
    return v;
  }

  // As it is written: "1 - 3t + t^2", "t^-1 + 2".
  std::string str() const {
    if (is_zero()) return "0";
    std::string out;
    for (const auto& [k, c] : terms_) {
      const bool negative = c.sign() < 0;
      const Integer magnitude = negative ? -c : c;
      out += out.empty() ? (negative ? "-" : "") : (negative ? " - " : " + ");
      const std::string power = k == 0 ? "" : k == 1 ? "t" : "t^" + std::to_string(k);
      if (k == 0 || !magnitude.is_one()) out += magnitude.str();
      out += power;
    }
    return out;
  }

 private:
  void add(long k, const Integer& c) {
    if (c.is_zero()) return;
    auto [it, fresh] = terms_.emplace(k, c);
    if (!fresh) {
      it->second = it->second + c;
      if (it->second.is_zero()) terms_.erase(it);
    }
  }
  std::map<long, Integer> terms_;
};

}  // namespace morphe::arithmetic
