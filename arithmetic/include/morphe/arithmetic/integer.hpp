// The integers ℤ, exactly and without bound.
//
// An Integer is held in one machine word while its value fits in a signed 64-bit integer,
// and otherwise as a sign and a magnitude in base 2^32 (least significant limb first). The
// representation is canonical: a value that fits a word is always held in one, so equality
// is equality of the parts. Arithmetic on words checks for overflow and moves to limbs only
// when a result needs them, which keeps the common case of homology (small entries) at the
// speed of machine arithmetic.
//
// The algorithms are the classical ones (Knuth 1997, The Art of Computer Programming 2,
// 3rd ed., §4.3.1): schoolbook addition, subtraction and multiplication, O(nm) on operands
// of n and m limbs; division by Algorithm D, O(nm), with the normalisation and the
// correction of the trial quotient as Warren 2012 (Hacker's Delight, 2nd ed., §9-2) writes
// it. Faster multiplication (Karatsuba, FFT) waits until a computation needs it.
//
// Division follows C++: a / b truncates toward zero and a % b has the sign of a. floor_div
// and mod give the Euclidean pair, 0 <= mod(a, b) < |b|. Division by zero throws
// std::domain_error. The word fast paths use GCC's and Clang's overflow builtins.

#pragma once

#include <algorithm>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace morphe::arithmetic {

class Integer {
 public:
  using Limbs = std::vector<std::uint32_t>;

  Integer() = default;
  Integer(long long v) : small_(v) {}
  Integer(long v) : small_(v) {}
  Integer(int v) : small_(v) {}
  Integer(unsigned v) : small_(v) {}
  Integer(unsigned long long v) { *this = make(false, from_u64(v)); }
  Integer(unsigned long v) : Integer(static_cast<unsigned long long>(v)) {}
  // From decimal digits, with an optional sign: "-12345678901234567890".
  explicit Integer(std::string_view text) { *this = parse(text); }

  static Integer parse(std::string_view text) {
    std::size_t i = 0;
    bool neg = false;
    if (i < text.size() && (text[i] == '-' || text[i] == '+')) neg = text[i++] == '-';
    if (i == text.size()) throw std::invalid_argument("morphe::arithmetic::Integer: no digits in \"" + std::string(text) + "\"");
    Limbs m;
    for (; i < text.size(); i++) {
      const char c = text[i];
      if (c < '0' || c > '9') throw std::invalid_argument("morphe::arithmetic::Integer: not a decimal integer: \"" + std::string(text) + "\"");
      mul_small(m, 10, static_cast<std::uint32_t>(c - '0'));
    }
    return make(neg, std::move(m));
  }

  int sign() const { return big_.empty() ? (small_ > 0) - (small_ < 0) : (neg_ ? -1 : 1); }
  bool is_zero() const { return big_.empty() && small_ == 0; }
  bool is_one() const { return big_.empty() && small_ == 1; }
  // Whether the value fits in a signed 64-bit integer.
  bool fits_int64() const { return big_.empty(); }
  long long to_int64() const {
    if (!big_.empty()) throw std::overflow_error("morphe::arithmetic::Integer: " + str() + " does not fit in 64 bits");
    return small_;
  }
  // The nearest double, to within one unit in the last place.
  double to_double() const {
    if (big_.empty()) return static_cast<double>(small_);
    const std::size_t n = big_.size();
    // The top three limbs carry at least 65 significant bits.
    double v = 0;
    for (std::size_t k = 0; k < 3 && k < n; k++) v += std::ldexp(static_cast<double>(big_[n - 1 - k]), static_cast<int>(32 * (n - 1 - k)));
    return neg_ ? -v : v;
  }
  // The number of bits of |value|; 0 for zero.
  std::size_t bit_length() const {
    const Limbs m = magnitude();
    if (m.empty()) return 0;
    std::size_t b = 32 * (m.size() - 1);
    for (std::uint32_t top = m.back(); top; top >>= 1) b++;
    return b;
  }
  std::string str() const {
    if (big_.empty()) return std::to_string(small_);
    Limbs m = big_;
    std::string digits;
    while (!m.empty()) {
      std::uint32_t r = div_small(m, 1000000000u);
      for (int k = 0; k < 9; k++) {
        digits.push_back(static_cast<char>('0' + r % 10));
        r /= 10;
        if (m.empty() && r == 0) break;
      }
    }
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    if (neg_) digits.push_back('-');
    std::reverse(digits.begin(), digits.end());
    return digits;
  }

  Integer operator-() const {
    if (big_.empty() && small_ != std::numeric_limits<long long>::min()) return Integer(-small_);
    return make(sign() > 0, magnitude());
  }
  Integer& operator+=(const Integer& b) { return *this = *this + b; }
  Integer& operator-=(const Integer& b) { return *this = *this - b; }
  Integer& operator*=(const Integer& b) { return *this = *this * b; }
  Integer& operator/=(const Integer& b) { return *this = *this / b; }
  Integer& operator%=(const Integer& b) { return *this = *this % b; }

  friend Integer operator+(const Integer& a, const Integer& b) {
    long long r;
    if (a.big_.empty() && b.big_.empty() && !__builtin_add_overflow(a.small_, b.small_, &r)) return Integer(r);
    return signed_add(a.sign() < 0, a.magnitude(), b.sign() < 0, b.magnitude());
  }
  friend Integer operator-(const Integer& a, const Integer& b) {
    long long r;
    if (a.big_.empty() && b.big_.empty() && !__builtin_sub_overflow(a.small_, b.small_, &r)) return Integer(r);
    return signed_add(a.sign() < 0, a.magnitude(), b.sign() > 0, b.magnitude());
  }
  friend Integer operator*(const Integer& a, const Integer& b) {
    long long r;
    if (a.big_.empty() && b.big_.empty() && !__builtin_mul_overflow(a.small_, b.small_, &r)) return Integer(r);
    return make((a.sign() < 0) != (b.sign() < 0), mul_mag(a.magnitude(), b.magnitude()));
  }
  // Truncated toward zero.
  friend Integer operator/(const Integer& a, const Integer& b) {
    Integer q, r;
    divmod(a, b, q, r);
    return q;
  }
  // With the sign of a.
  friend Integer operator%(const Integer& a, const Integer& b) {
    Integer q, r;
    divmod(a, b, q, r);
    return r;
  }
  // a = q b + r, q truncated toward zero, r with the sign of a and |r| < |b|.
  friend void divmod(const Integer& a, const Integer& b, Integer& q, Integer& r) {
    if (b.is_zero()) throw std::domain_error("morphe::arithmetic::Integer: division by zero");
    if (a.big_.empty() && b.big_.empty() && !(a.small_ == std::numeric_limits<long long>::min() && b.small_ == -1)) {
      q = Integer(a.small_ / b.small_);
      r = Integer(a.small_ % b.small_);
      return;
    }
    Limbs qm, rm;
    divmod_mag(a.magnitude(), b.magnitude(), qm, rm);
    q = make((a.sign() < 0) != (b.sign() < 0), std::move(qm));
    r = make(a.sign() < 0, std::move(rm));
  }

  friend bool operator==(const Integer& a, const Integer& b) {
    return a.small_ == b.small_ && a.neg_ == b.neg_ && a.big_ == b.big_;
  }
  friend std::strong_ordering operator<=>(const Integer& a, const Integer& b) {
    if (a.big_.empty() && b.big_.empty()) return a.small_ <=> b.small_;
    const int sa = a.sign(), sb = b.sign();
    if (sa != sb) return sa <=> sb;
    const int c = cmp_mag(a.magnitude(), b.magnitude());
    return (sa < 0 ? -c : c) <=> 0;
  }

  // |value| as limbs, least significant first, with no leading zero limb.
  Limbs magnitude() const {
    if (!big_.empty()) return big_;
    const unsigned long long m = small_ < 0 ? static_cast<unsigned long long>(-(small_ + 1)) + 1 : static_cast<unsigned long long>(small_);
    return from_u64(m);
  }
  // The value with sign neg and magnitude m (any leading zeros), in canonical form.
  static Integer make(bool neg, Limbs m) {
    trim(m);
    Integer x;
    if (m.size() <= 2) {
      const unsigned long long v = m.empty() ? 0 : (m.size() == 1 ? m[0] : (static_cast<unsigned long long>(m[1]) << 32 | m[0]));
      constexpr unsigned long long top = 1ull << 63;
      if (!neg && v < top) { x.small_ = static_cast<long long>(v); return x; }
      if (neg && v <= top) { x.small_ = v == top ? std::numeric_limits<long long>::min() : -static_cast<long long>(v); return x; }
    }
    x.small_ = 0;
    x.neg_ = neg;
    x.big_ = std::move(m);
    return x;
  }

 private:
  static void trim(Limbs& m) {
    while (!m.empty() && m.back() == 0) m.pop_back();
  }
  static Limbs from_u64(unsigned long long v) {
    Limbs m;
    if (v) m.push_back(static_cast<std::uint32_t>(v));
    if (v >> 32) m.push_back(static_cast<std::uint32_t>(v >> 32));
    return m;
  }
  static int cmp_mag(const Limbs& a, const Limbs& b) {
    if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
    for (std::size_t i = a.size(); i-- > 0;)
      if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
  }
  static Limbs add_mag(const Limbs& a, const Limbs& b) {
    const Limbs& x = a.size() >= b.size() ? a : b;
    const Limbs& y = a.size() >= b.size() ? b : a;
    Limbs s(x.size() + 1);
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
      carry += static_cast<std::uint64_t>(x[i]) + (i < y.size() ? y[i] : 0);
      s[i] = static_cast<std::uint32_t>(carry);
      carry >>= 32;
    }
    s[x.size()] = static_cast<std::uint32_t>(carry);
    trim(s);
    return s;
  }
  // a - b for |a| >= |b|.
  static Limbs sub_mag(const Limbs& a, const Limbs& b) {
    Limbs d(a.size());
    std::int64_t borrow = 0;
    for (std::size_t i = 0; i < a.size(); i++) {
      std::int64_t t = static_cast<std::int64_t>(a[i]) - (i < b.size() ? b[i] : 0) - borrow;
      borrow = t < 0;
      if (t < 0) t += std::int64_t(1) << 32;
      d[i] = static_cast<std::uint32_t>(t);
    }
    trim(d);
    return d;
  }
  static Integer signed_add(bool na, const Limbs& a, bool nb, const Limbs& b) {
    if (na == nb) return make(na, add_mag(a, b));
    const int c = cmp_mag(a, b);
    if (c == 0) return Integer(0);
    return c > 0 ? make(na, sub_mag(a, b)) : make(nb, sub_mag(b, a));
  }
  static Limbs mul_mag(const Limbs& a, const Limbs& b) {
    if (a.empty() || b.empty()) return {};
    Limbs p(a.size() + b.size(), 0);
    for (std::size_t i = 0; i < a.size(); i++) {
      std::uint64_t carry = 0;
      for (std::size_t j = 0; j < b.size(); j++) {
        carry += static_cast<std::uint64_t>(a[i]) * b[j] + p[i + j];
        p[i + j] = static_cast<std::uint32_t>(carry);
        carry >>= 32;
      }
      p[i + b.size()] = static_cast<std::uint32_t>(carry);
    }
    trim(p);
    return p;
  }
  // m <- m * f + add.
  static void mul_small(Limbs& m, std::uint32_t f, std::uint32_t add) {
    std::uint64_t carry = add;
    for (auto& x : m) {
      carry += static_cast<std::uint64_t>(x) * f;
      x = static_cast<std::uint32_t>(carry);
      carry >>= 32;
    }
    if (carry) m.push_back(static_cast<std::uint32_t>(carry));
  }
  // m <- m / d; returns m mod d.
  static std::uint32_t div_small(Limbs& m, std::uint32_t d) {
    std::uint64_t r = 0;
    for (std::size_t i = m.size(); i-- > 0;) {
      const std::uint64_t cur = r << 32 | m[i];
      m[i] = static_cast<std::uint32_t>(cur / d);
      r = cur % d;
    }
    trim(m);
    return static_cast<std::uint32_t>(r);
  }
  static int leading_zeros(std::uint32_t x) {
    int n = 0;
    while (!(x & 0x80000000u)) { x <<= 1; n++; }
    return n;
  }
  // Knuth's Algorithm D: u = q v + r, 0 <= r < v, for v nonzero.
  static void divmod_mag(const Limbs& u, const Limbs& v, Limbs& q, Limbs& r) {
    if (cmp_mag(u, v) < 0) { q.clear(); r = u; return; }
    if (v.size() == 1) {
      q = u;
      const std::uint32_t rem = div_small(q, v[0]);
      r = rem ? Limbs{rem} : Limbs{};
      return;
    }
    const std::size_t n = v.size(), m = u.size() - v.size();
    const int s = leading_zeros(v.back());
    // Normalise so that the top limb of the divisor has its high bit set.
    Limbs vn(n), un(u.size() + 1);
    for (std::size_t i = n - 1; i > 0; i--) vn[i] = s ? (v[i] << s) | (v[i - 1] >> (32 - s)) : v[i];
    vn[0] = v[0] << s;
    un[u.size()] = s ? u.back() >> (32 - s) : 0;
    for (std::size_t i = u.size() - 1; i > 0; i--) un[i] = s ? (u[i] << s) | (u[i - 1] >> (32 - s)) : u[i];
    un[0] = u[0] << s;
    q.assign(m + 1, 0);
    constexpr std::uint64_t b = 1ull << 32;
    for (std::size_t j = m + 1; j-- > 0;) {
      // The trial quotient from the top two limbs, corrected at most twice.
      const std::uint64_t num = static_cast<std::uint64_t>(un[j + n]) << 32 | un[j + n - 1];
      std::uint64_t qhat = num / vn[n - 1], rhat = num % vn[n - 1];
      while (qhat >= b || qhat * vn[n - 2] > (rhat << 32 | un[j + n - 2])) {
        qhat--;
        rhat += vn[n - 1];
        if (rhat >= b) break;
      }
      // Multiply and subtract.
      std::int64_t borrow = 0;
      std::uint64_t carry = 0;
      for (std::size_t i = 0; i < n; i++) {
        const std::uint64_t p = qhat * vn[i] + carry;
        carry = p >> 32;
        const std::int64_t t = static_cast<std::int64_t>(un[i + j]) - borrow - static_cast<std::int64_t>(p & 0xffffffffu);
        un[i + j] = static_cast<std::uint32_t>(t);
        borrow = t < 0 ? 1 : 0;
      }
      const std::int64_t t = static_cast<std::int64_t>(un[j + n]) - borrow - static_cast<std::int64_t>(carry);
      un[j + n] = static_cast<std::uint32_t>(t);
      if (t < 0) {
        // The trial quotient was one too large: add the divisor back.
        qhat--;
        std::uint64_t c = 0;
        for (std::size_t i = 0; i < n; i++) {
          c += static_cast<std::uint64_t>(un[i + j]) + vn[i];
          un[i + j] = static_cast<std::uint32_t>(c);
          c >>= 32;
        }
        un[j + n] = static_cast<std::uint32_t>(un[j + n] + c);
      }
      q[j] = static_cast<std::uint32_t>(qhat);
    }
    trim(q);
    // Unnormalise the remainder.
    r.assign(n, 0);
    for (std::size_t i = 0; i < n; i++) r[i] = s ? (un[i] >> s) | (un[i + 1] << (32 - s)) : un[i];
    trim(r);
  }

  long long small_ = 0;
  bool neg_ = false;
  Limbs big_;
};

inline Integer abs(const Integer& a) { return a.sign() < 0 ? -a : a; }

// The Euclidean pair: a = floor_div(a, b) b + mod(a, b), with 0 <= mod(a, b) < |b|.
inline Integer mod(const Integer& a, const Integer& b) {
  Integer r = a % b;
  if (r.sign() < 0) r += abs(b);
  return r;
}
inline Integer floor_div(const Integer& a, const Integer& b) { return (a - mod(a, b)) / b; }

// a / b when b divides a; throws std::domain_error otherwise.
inline Integer exact_div(const Integer& a, const Integer& b) {
  Integer q, r;
  divmod(a, b, q, r);
  if (!r.is_zero()) throw std::domain_error("morphe::arithmetic::exact_div: " + b.str() + " does not divide " + a.str());
  return q;
}

// The greatest common divisor, nonnegative; gcd(0, 0) = 0. Euclid's algorithm.
inline Integer gcd(Integer a, Integer b) {
  a = abs(a);
  b = abs(b);
  while (!b.is_zero()) {
    Integer r = a % b;
    a = std::move(b);
    b = std::move(r);
  }
  return a;
}

// g = gcd(a, b) = s a + t b, by the extended Euclidean algorithm.
struct Bezout {
  Integer g, s, t;
};
inline Bezout extended_gcd(const Integer& a, const Integer& b) {
  Integer r0 = a, r1 = b, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
  while (!r1.is_zero()) {
    Integer q = r0 / r1;
    Integer r2 = r0 - q * r1, s2 = s0 - q * s1, t2 = t0 - q * t1;
    r0 = std::move(r1), r1 = std::move(r2);
    s0 = std::move(s1), s1 = std::move(s2);
    t0 = std::move(t1), t1 = std::move(t2);
  }
  if (r0.sign() < 0) return {-r0, -s0, -t0};
  return {r0, s0, t0};
}

// a^e by repeated squaring.
inline Integer pow(Integer a, unsigned long long e) {
  Integer r = 1;
  while (e) {
    if (e & 1) r *= a;
    e >>= 1;
    if (e) a *= a;
  }
  return r;
}

// n!, the double factorial n!! (with (-1)!! = 0!! = 1), and the binomial coefficient.
inline Integer factorial(unsigned long long n) {
  Integer r = 1;
  for (unsigned long long k = 2; k <= n; k++) r *= Integer(k);
  return r;
}
inline Integer double_factorial(long long n) {
  Integer r = 1;
  for (long long k = n; k > 1; k -= 2) r *= Integer(k);
  return r;
}
inline Integer binomial(long long n, long long k) {
  if (k < 0 || n < 0 || k > n) return 0;
  k = std::min(k, n - k);
  Integer r = 1;
  for (long long i = 1; i <= k; i++) r = exact_div(r * Integer(n - k + i), Integer(i));
  return r;
}

}  // namespace morphe::arithmetic

template <>
struct std::hash<morphe::arithmetic::Integer> {
  std::size_t operator()(const morphe::arithmetic::Integer& x) const noexcept {
    if (x.fits_int64()) return std::hash<long long>{}(x.to_int64());
    std::size_t h = x.sign() < 0 ? 0x9e3779b97f4a7c15ull : 0;
    for (auto limb : x.magnitude()) h = (h ^ limb) * 0x100000001b3ull;
    return h;
  }
};
