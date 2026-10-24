// Semirings, the exact ones that alignment and weighted automata use, and the homomorphisms
// between them.
//
// A semiring (K, ⊕, ⊗, 0, 1) is a set with two associative operations, ⊕ commutative with
// identity 0, ⊗ with identity 1, ⊗ distributing over ⊕ on both sides, and 0 annihilating:
// 0 ⊗ a = a ⊗ 0 = 0 (Golan 1999, Semirings and their Applications, Kluwer, ch. 1; Mohri 2002,
// Semiring frameworks and algorithms for shortest-distance problems, J. Automata, Languages and
// Combinatorics 7, 321–350, §2). It is idempotent when a ⊕ a = a, and selective when
// a ⊕ b ∈ {a, b}: then ⊕ chooses, a sum over paths is attained by one path, and that path can
// be read back (traceback). A semiring here is a type S with
//
//   S::value_type, S::zero(), S::one(), S::plus(a, b), S::times(a, b), S::equal(a, b),
//   S::idempotent, S::selective, S::commutative
//
// checked by the concept Semiring. The semirings of this header:
//
//   Boolean          ({0, 1}, ∨, ∧): whether a path exists. Selective.
//   MaxPlus<T>       (T ∪ {−∞}, max, +): the optimal score (Viterbi). Selective.
//   MinPlus<T>       (T ∪ {+∞}, min, +): the least cost, edit distance. Selective.
//   Arithmetic<T>    (T, +, ×): with T = Integer, Counting, the number of paths over ℕ, exact
//                    without bound; with T = double, the probability semiring of a pair HMM's
//                    forward algorithm (exact only up to rounding: the one inexact semiring here,
//                    for structures built from real-valued data, CLAUDE.md invariant 3).
//   MaxPlusCount     the tropical semiring with multiplicities: pairs (s, c), s ∈ ℤ ∪ {−∞} a
//                    score and c ∈ ℕ the number of paths attaining it, with
//                    (s, c) ⊕ (t, d) = (s, c) if s > t, (t, d) if t > s, (s, c + d) if s = t,
//                    (s, c) ⊗ (t, d) = (s + t, c d), 0 = (−∞, 0), 1 = (0, 1); the carrier is the
//                    pairs with c ≥ 1 and the zero. It is the image of the Laurent polynomials
//                    ℕ[t, t⁻¹] under the leading term, p ↦ (deg p, its leading coefficient),
//                    which is a homomorphism because ℕ has neither negatives nor zero divisors
//                    (Pachter and Sturmfels 2005, Algebraic Statistics
//                    for Computational Biology, Cambridge, §2.1, the "tropical semiring" and its
//                    use for counting optimal alignments). It is not idempotent:
//                    (s, c) ⊕ (s, c) = (s, 2c).
//   Product<A, B>    the componentwise product; idempotent and commutative when both are,
//                    never claimed selective.
//
// For integral T the infinities of MaxPlus and MinPlus are the type's lowest and highest
// values, reserved: a sum that would reach or pass them throws std::overflow_error, so the
// arithmetic is exact or refused, never wrapped. For floating T they are ±infinity.
//
// Homomorphisms. A homomorphism h : S → R preserves both operations and both identities. A
// morphism here is a function object with source and target types (the concept SemiringMorphism;
// words.hpp's Morphism is a morphism of free monoids), checked on samples by
// violated_morphism_law. The maps of this header, and what they are:
//
//   ForgetCount      MaxPlusCount → MaxPlus<int64>, (s, c) ↦ s: a homomorphism (the projection
//                    to the leading degree).
//   EmbedCount       Counting → MaxPlusCount, n ↦ (0, n) for n ≥ 1, 0 ↦ (−∞, 0): an injective
//                    homomorphism, so ℕ is the sub-semiring of scores 0.
//   Negate<T>        MaxPlus<T> → MinPlus<T>, x ↦ −x: an isomorphism (max(a, b) = −min(−a, −b)).
//   Finite<T>        MaxPlus<T> → Boolean, x ↦ [x ≠ −∞]: a homomorphism (a path exists).
//   Support          Counting → Boolean, n ↦ [n ≠ 0]: a homomorphism, because ℕ has no
//                    negatives (n + m = 0 only when n = m = 0); on ℤ it would not be.
//   First, Second    the projections of a Product: homomorphisms.
//   not a homomorphism: (s, c) ↦ c from MaxPlusCount to Counting is multiplicative but not
//                    additive ((1, 1) ⊕ (0, 1) = (1, 1), while 1 + 1 ≠ 1); the checks exhibit it.
//
// Because a sum over the paths of a graph is a polynomial expression in its edge weights, a
// homomorphism carries the sum in S to the sum in R of the images of the weights: counts of
// optimal paths project to optimal scores, and the polytope algebra (polytope.hpp) evaluates to
// the tropical semiring at every parameter. This is why one generic dynamic programme serves
// every scoring.
//
// Complexity: the operations are O(1), except on Integer (O(limbs), arithmetic/integer.hpp).
// violated_law is O(k³) in the number of samples, violated_morphism_law O(k²).

#pragma once

#include <morphe/arithmetic/integer.hpp>

#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace morphe::formal {

template <class S>
concept Semiring = requires(const typename S::value_type& a, const typename S::value_type& b) {
  { S::zero() } -> std::convertible_to<typename S::value_type>;
  { S::one() } -> std::convertible_to<typename S::value_type>;
  { S::plus(a, b) } -> std::convertible_to<typename S::value_type>;
  { S::times(a, b) } -> std::convertible_to<typename S::value_type>;
  { S::equal(a, b) } -> std::convertible_to<bool>;
  { S::idempotent } -> std::convertible_to<bool>;
  { S::selective } -> std::convertible_to<bool>;
  { S::commutative } -> std::convertible_to<bool>;
};

// A semiring in which ⊕ chooses one of its arguments: sums are attained, and traceback defined.
template <class S>
concept SelectiveSemiring = Semiring<S> && S::selective;

// ---- Boolean ------------------------------------------------------------------------------------
struct Boolean {
  using value_type = bool;
  static constexpr bool idempotent = true, selective = true, commutative = true;
  static bool zero() { return false; }
  static bool one() { return true; }
  static bool plus(bool a, bool b) { return a || b; }
  static bool times(bool a, bool b) { return a && b; }
  static bool equal(bool a, bool b) { return a == b; }
};

// ---- tropical semirings -----------------------------------------------------------------------------
namespace detail {

// a + b for integers, refusing a result at or beyond the reserved infinities.
template <class T>
T exact_add(T a, T b) {
  T r;
  if (__builtin_add_overflow(a, b, &r) || r == std::numeric_limits<T>::min() || r == std::numeric_limits<T>::max())
    throw std::overflow_error("morphe::formal: a tropical sum beyond the range of its integers");
  return r;
}

}  // namespace detail

template <class T = std::int64_t>
struct MaxPlus {
  static_assert(std::is_arithmetic_v<T> && std::is_signed_v<T>, "MaxPlus over a signed arithmetic type");
  using value_type = T;
  static constexpr bool idempotent = true, selective = true, commutative = true;
  static T zero() {
    if constexpr (std::is_floating_point_v<T>) return -std::numeric_limits<T>::infinity();
    else return std::numeric_limits<T>::min();
  }
  static T one() { return T(0); }
  static T plus(T a, T b) { return a < b ? b : a; }
  static T times(T a, T b) {
    if constexpr (std::is_floating_point_v<T>) return a + b;
    else return (a == zero() || b == zero()) ? zero() : detail::exact_add(a, b);
  }
  static bool equal(T a, T b) { return a == b; }
  static bool is_zero(T a) { return a == zero(); }
};

template <class T = std::int64_t>
struct MinPlus {
  static_assert(std::is_arithmetic_v<T> && std::is_signed_v<T>, "MinPlus over a signed arithmetic type");
  using value_type = T;
  static constexpr bool idempotent = true, selective = true, commutative = true;
  static T zero() {
    if constexpr (std::is_floating_point_v<T>) return std::numeric_limits<T>::infinity();
    else return std::numeric_limits<T>::max();
  }
  static T one() { return T(0); }
  static T plus(T a, T b) { return b < a ? b : a; }
  static T times(T a, T b) {
    if constexpr (std::is_floating_point_v<T>) return a + b;
    else return (a == zero() || b == zero()) ? zero() : detail::exact_add(a, b);
  }
  static bool equal(T a, T b) { return a == b; }
  static bool is_zero(T a) { return a == zero(); }
};

// ---- (T, +, ×) --------------------------------------------------------------------------------------
template <class T>
struct Arithmetic {
  using value_type = T;
  static constexpr bool idempotent = false, selective = false, commutative = true;
  static T zero() { return T(0); }
  static T one() { return T(1); }
  static T plus(const T& a, const T& b) { return a + b; }
  static T times(const T& a, const T& b) { return a * b; }
  static bool equal(const T& a, const T& b) { return a == b; }
};

// The number of paths, exactly: (ℕ, +, ×) on morphe's integers.
using Counting = Arithmetic<arithmetic::Integer>;

// ---- the tropical semiring with multiplicities ---------------------------------------------------------
struct ScoreCount {
  std::int64_t score = std::numeric_limits<std::int64_t>::min();  // −∞ for the zero
  arithmetic::Integer count = 0;                                  // the number of paths attaining it
  friend bool operator==(const ScoreCount&, const ScoreCount&) = default;
  std::string str() const {
    return "(" + (score == std::numeric_limits<std::int64_t>::min() ? std::string("-inf") : std::to_string(score)) + ", " + count.str() + ")";
  }
};

struct MaxPlusCount {
  using value_type = ScoreCount;
  static constexpr bool idempotent = false, selective = false, commutative = true;
  static ScoreCount zero() { return {}; }
  static ScoreCount one() { return {0, 1}; }
  // An element (s, c) with c ≥ 1 and s finite: the weight of one path of score s.
  static ScoreCount single(std::int64_t s) { return {s, 1}; }
  static bool is_zero(const ScoreCount& a) { return a.count.is_zero(); }
  static ScoreCount plus(const ScoreCount& a, const ScoreCount& b) {
    if (is_zero(a)) return b;
    if (is_zero(b)) return a;
    if (a.score > b.score) return a;
    if (b.score > a.score) return b;
    return {a.score, a.count + b.count};
  }
  static ScoreCount times(const ScoreCount& a, const ScoreCount& b) {
    if (is_zero(a) || is_zero(b)) return zero();
    return {detail::exact_add(a.score, b.score), a.count * b.count};
  }
  static bool equal(const ScoreCount& a, const ScoreCount& b) { return a == b; }
};

// ---- products -------------------------------------------------------------------------------------------
template <Semiring A, Semiring B>
struct Product {
  using value_type = std::pair<typename A::value_type, typename B::value_type>;
  static constexpr bool idempotent = A::idempotent && B::idempotent;
  static constexpr bool selective = false;
  static constexpr bool commutative = A::commutative && B::commutative;
  static value_type zero() { return {A::zero(), B::zero()}; }
  static value_type one() { return {A::one(), B::one()}; }
  static value_type plus(const value_type& a, const value_type& b) { return {A::plus(a.first, b.first), B::plus(a.second, b.second)}; }
  static value_type times(const value_type& a, const value_type& b) { return {A::times(a.first, b.first), B::times(a.second, b.second)}; }
  static bool equal(const value_type& a, const value_type& b) { return A::equal(a.first, b.first) && B::equal(a.second, b.second); }
};

// ---- homomorphisms ----------------------------------------------------------------------------------------
// A morphism of semirings: a function object from Source's values to Target's.
template <class H>
concept SemiringMorphism = Semiring<typename H::source> && Semiring<typename H::target> && requires(const H& h, const typename H::source::value_type& a) {
  { h(a) } -> std::convertible_to<typename H::target::value_type>;
};

struct ForgetCount {
  using source = MaxPlusCount;
  using target = MaxPlus<std::int64_t>;
  std::int64_t operator()(const ScoreCount& a) const { return MaxPlusCount::is_zero(a) ? target::zero() : a.score; }
};

struct EmbedCount {
  using source = Counting;
  using target = MaxPlusCount;
  ScoreCount operator()(const arithmetic::Integer& n) const { return n.is_zero() ? MaxPlusCount::zero() : ScoreCount{0, n}; }
};

template <class T = std::int64_t>
struct Negate {
  using source = MaxPlus<T>;
  using target = MinPlus<T>;
  T operator()(T a) const {
    if (a == source::zero()) return target::zero();
    return -a;
  }
};

template <class T = std::int64_t>
struct Finite {
  using source = MaxPlus<T>;
  using target = Boolean;
  bool operator()(T a) const { return a != source::zero(); }
};

struct Support {
  using source = Counting;
  using target = Boolean;
  bool operator()(const arithmetic::Integer& n) const { return !n.is_zero(); }
};

template <Semiring A, Semiring B>
struct First {
  using source = Product<A, B>;
  using target = A;
  typename A::value_type operator()(const typename source::value_type& a) const { return a.first; }
};
template <Semiring A, Semiring B>
struct Second {
  using source = Product<A, B>;
  using target = B;
  typename B::value_type operator()(const typename source::value_type& a) const { return a.second; }
};

// ---- the laws, on samples --------------------------------------------------------------------------------
// The first law of a semiring that fails on the samples (all triples of them), or nothing:
// associativity and commutativity of ⊕, its identity 0; associativity of ⊗, its identity 1 on
// both sides; both distributive laws; 0 annihilating on both sides; and the properties S claims
// (idempotent, selective, commutative ⊗).
template <Semiring S>
std::optional<std::string> violated_law(const std::vector<typename S::value_type>& samples) {
  const auto eq = [](const auto& a, const auto& b) { return S::equal(a, b); };
  const auto z = S::zero(), u = S::one();
  for (const auto& a : samples) {
    if (!eq(S::plus(a, z), a) || !eq(S::plus(z, a), a)) return "0 is not the identity of ⊕";
    if (!eq(S::times(a, u), a) || !eq(S::times(u, a), a)) return "1 is not the identity of ⊗";
    if (!eq(S::times(a, z), z) || !eq(S::times(z, a), z)) return "0 does not annihilate";
    if (S::idempotent && !eq(S::plus(a, a), a)) return "⊕ is not idempotent";
    for (const auto& b : samples) {
      if (!eq(S::plus(a, b), S::plus(b, a))) return "⊕ is not commutative";
      if (S::commutative && !eq(S::times(a, b), S::times(b, a))) return "⊗ is not commutative";
      if (S::selective && !eq(S::plus(a, b), a) && !eq(S::plus(a, b), b)) return "⊕ is not selective";
      for (const auto& c : samples) {
        if (!eq(S::plus(S::plus(a, b), c), S::plus(a, S::plus(b, c)))) return "⊕ is not associative";
        if (!eq(S::times(S::times(a, b), c), S::times(a, S::times(b, c)))) return "⊗ is not associative";
        if (!eq(S::times(a, S::plus(b, c)), S::plus(S::times(a, b), S::times(a, c)))) return "⊗ does not distribute over ⊕ on the left";
        if (!eq(S::times(S::plus(a, b), c), S::plus(S::times(a, c), S::times(b, c)))) return "⊗ does not distribute over ⊕ on the right";
      }
    }
  }
  return std::nullopt;
}

// The first law of a homomorphism h : S → R that fails on the samples, or nothing:
// h(a ⊕ b) = h(a) ⊕ h(b), h(a ⊗ b) = h(a) ⊗ h(b), h(0) = 0, h(1) = 1.
template <SemiringMorphism H>
std::optional<std::string> violated_morphism_law(const H& h, const std::vector<typename H::source::value_type>& samples) {
  using S = typename H::source;
  using R = typename H::target;
  if (!R::equal(h(S::zero()), R::zero())) return "h(0) ≠ 0";
  if (!R::equal(h(S::one()), R::one())) return "h(1) ≠ 1";
  for (const auto& a : samples)
    for (const auto& b : samples) {
      if (!R::equal(h(S::plus(a, b)), R::plus(h(a), h(b)))) return "h(a ⊕ b) ≠ h(a) ⊕ h(b)";
      if (!R::equal(h(S::times(a, b)), R::times(h(a), h(b)))) return "h(a ⊗ b) ≠ h(a) ⊗ h(b)";
    }
  return std::nullopt;
}

}  // namespace morphe::formal
