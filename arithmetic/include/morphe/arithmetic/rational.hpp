// The rationals ℚ, exactly.
//
// A Rational is a pair of Integers (numerator, denominator) in lowest terms with a positive
// denominator, so that equality is equality of the parts. Each operation reduces its result
// by one gcd (Knuth 1997, §4.5.1, with the cross-cancellation of Henrici left for later).

#pragma once

#include <compare>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "integer.hpp"

namespace morphe::arithmetic {

class Rational {
 public:
  Rational() : num_(0), den_(1) {}
  Rational(Integer n) : num_(std::move(n)), den_(1) {}
  Rational(long long n) : num_(n), den_(1) {}
  Rational(int n) : num_(n), den_(1) {}
  Rational(Integer n, Integer d) : num_(std::move(n)), den_(std::move(d)) { normalise(); }
  // "p/q" or "p".
  static Rational parse(std::string_view text) {
    const auto slash = text.find('/');
    if (slash == std::string_view::npos) return Rational(Integer::parse(text));
    return Rational(Integer::parse(text.substr(0, slash)), Integer::parse(text.substr(slash + 1)));
  }

  const Integer& numerator() const { return num_; }
  const Integer& denominator() const { return den_; }
  int sign() const { return num_.sign(); }
  bool is_zero() const { return num_.is_zero(); }
  bool is_integer() const { return den_.is_one(); }
  std::string str() const { return is_integer() ? num_.str() : num_.str() + "/" + den_.str(); }
  // The nearest double, to within a few units in the last place, for any size of the parts.
  double to_double() const {
    if (num_.fits_int64() && den_.fits_int64() && num_.bit_length() <= 53 && den_.bit_length() <= 53)
      return num_.to_double() / den_.to_double();
    // Scale so that the quotient has 64 significant bits, then divide exactly.
    const long long shift = 64 - static_cast<long long>(num_.bit_length()) + static_cast<long long>(den_.bit_length());
    const Integer n = shift > 0 ? num_ * pow(Integer(2), static_cast<unsigned long long>(shift)) : num_;
    const Integer d = shift < 0 ? den_ * pow(Integer(2), static_cast<unsigned long long>(-shift)) : den_;
    return std::ldexp((n / d).to_double(), static_cast<int>(-shift));
  }

  Rational operator-() const { return Rational(-num_, den_, true); }
  friend Rational operator+(const Rational& a, const Rational& b) { return Rational(a.num_ * b.den_ + b.num_ * a.den_, a.den_ * b.den_); }
  friend Rational operator-(const Rational& a, const Rational& b) { return Rational(a.num_ * b.den_ - b.num_ * a.den_, a.den_ * b.den_); }
  friend Rational operator*(const Rational& a, const Rational& b) { return Rational(a.num_ * b.num_, a.den_ * b.den_); }
  friend Rational operator/(const Rational& a, const Rational& b) {
    if (b.is_zero()) throw std::domain_error("morphe::arithmetic::Rational: division by zero");
    return Rational(a.num_ * b.den_, a.den_ * b.num_);
  }
  Rational& operator+=(const Rational& b) { return *this = *this + b; }
  Rational& operator-=(const Rational& b) { return *this = *this - b; }
  Rational& operator*=(const Rational& b) { return *this = *this * b; }
  Rational& operator/=(const Rational& b) { return *this = *this / b; }

  friend bool operator==(const Rational& a, const Rational& b) { return a.num_ == b.num_ && a.den_ == b.den_; }
  friend std::strong_ordering operator<=>(const Rational& a, const Rational& b) { return a.num_ * b.den_ <=> b.num_ * a.den_; }

 private:
  Rational(Integer n, Integer d, bool) : num_(std::move(n)), den_(std::move(d)) {}
  void normalise() {
    if (den_.is_zero()) throw std::domain_error("morphe::arithmetic::Rational: zero denominator");
    if (den_.sign() < 0) { num_ = -num_; den_ = -den_; }
    const Integer g = gcd(num_, den_);
    if (!g.is_one() && !g.is_zero()) { num_ = exact_div(num_, g); den_ = exact_div(den_, g); }
    if (num_.is_zero()) den_ = 1;
  }

  Integer num_, den_;
};

inline Rational abs(const Rational& a) { return a.sign() < 0 ? -a : a; }

// a^e for an integer exponent e (a nonzero when e < 0).
inline Rational pow(const Rational& a, long long e) {
  if (e >= 0) return Rational(pow(a.numerator(), static_cast<unsigned long long>(e)), pow(a.denominator(), static_cast<unsigned long long>(e)));
  if (a.is_zero()) throw std::domain_error("morphe::arithmetic::pow: zero to a negative power");
  return Rational(pow(a.denominator(), static_cast<unsigned long long>(-e)), pow(a.numerator(), static_cast<unsigned long long>(-e)));
}

}  // namespace morphe::arithmetic
