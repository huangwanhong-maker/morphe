// Checks of morphe/formal/semiring.hpp: the semiring laws on samples of Boolean, MaxPlus,
// MinPlus, Counting, the tropical semiring with multiplicities and a product; that the law
// checker refuses (max, ×) over ℤ, which is not a semiring; the homomorphisms the header names,
// each law on samples, and the one map it names as not a homomorphism refused; exact counts
// beyond 64 bits; and tropical sums that would reach the reserved infinities refused.

#include <morphe/check.hpp>
#include <morphe/formal/semiring.hpp>

#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace mf = morphe::formal;
using morphe::arithmetic::Integer;

namespace {

std::string verdict(const std::optional<std::string>& v) { return v ? *v : "every law holds"; }

// (ℤ, max, ×): ⊗ does not distribute over max when a factor is negative.
struct MaxTimes {
  using value_type = std::int64_t;
  static constexpr bool idempotent = true, selective = true, commutative = true;
  static std::int64_t zero() { return std::numeric_limits<std::int64_t>::min(); }
  static std::int64_t one() { return 1; }
  static std::int64_t plus(std::int64_t a, std::int64_t b) { return a < b ? b : a; }
  static std::int64_t times(std::int64_t a, std::int64_t b) { return (a == zero() || b == zero()) ? zero() : a * b; }
  static bool equal(std::int64_t a, std::int64_t b) { return a == b; }
};

// (s, c) ↦ c, which the header says is not a homomorphism.
struct CountOnly {
  using source = mf::MaxPlusCount;
  using target = mf::Counting;
  Integer operator()(const mf::ScoreCount& a) const { return a.count; }
};

}  // namespace

int main() {
  morphe::Suite s;
  const auto ninf = mf::MaxPlus<>::zero(), pinf = mf::MinPlus<>::zero();
  const Integer big = morphe::arithmetic::pow(Integer(2), 70) + Integer(3);

  // ---- the laws ------------------------------------------------------------------------------------
  {
    const auto v = mf::violated_law<mf::Boolean>({false, true});
    s.require("semiring: Boolean ({0, 1}, ∨, ∧) satisfies every law, selective", !v, verdict(v));
  }
  {
    const auto v = mf::violated_law<mf::MaxPlus<>>({ninf, -7, -1, 0, 3, 12, 1000000});
    s.require("semiring: (ℤ ∪ {−∞}, max, +) satisfies every law, selective, on 7 samples", !v, verdict(v));
  }
  {
    const auto v = mf::violated_law<mf::MinPlus<>>({pinf, -7, -1, 0, 3, 12, 1000000});
    s.require("semiring: (ℤ ∪ {+∞}, min, +) satisfies every law, selective, on 7 samples", !v, verdict(v));
  }
  {
    const auto v = mf::violated_law<mf::MaxPlus<double>>({-std::numeric_limits<double>::infinity(), -2.5, 0.0, 0.25, 7.0});
    s.require("semiring: (ℝ ∪ {−∞}, max, +) on dyadic samples (exact in binary) satisfies every law", !v, verdict(v));
  }
  {
    const auto v = mf::violated_law<mf::Counting>({Integer(0), Integer(1), Integer(2), Integer(7), big});
    s.require("semiring: Counting (ℕ, +, ×) on morphe's integers satisfies every law, 2^70 + 3 among the samples", !v, verdict(v));
  }
  const std::vector<mf::ScoreCount> counted = {mf::MaxPlusCount::zero(), {0, 1}, {3, 2}, {3, 5}, {-4, 1}, {12, big}, {-4, 3}};
  {
    const auto v = mf::violated_law<mf::MaxPlusCount>(counted);
    s.require("semiring: the tropical semiring with multiplicities satisfies every law on 7 samples", !v, verdict(v));
    s.require("semiring: it is not idempotent: (3, 2) ⊕ (3, 2) = (3, 4)", mf::MaxPlusCount::plus({3, 2}, {3, 2}) == mf::ScoreCount{3, 4});
  }
  {
    using P = mf::Product<mf::MaxPlus<>, mf::Counting>;
    const auto v = mf::violated_law<P>({P::zero(), P::one(), {3, Integer(2)}, {ninf, Integer(5)}, {-2, big}});
    s.require("semiring: the product MaxPlus × Counting satisfies every law, componentwise", !v, verdict(v));
  }
  {
    const auto v = mf::violated_law<MaxTimes>({-3, -1, 0, 2, 5});
    s.require("semiring: the law checker refuses (ℤ, max, ×), where max(a, b) × c ≠ max(ac, bc) for c < 0", v.has_value(), verdict(v));
  }

  // ---- homomorphisms ------------------------------------------------------------------------------------
  {
    const auto v = mf::violated_morphism_law(mf::ForgetCount{}, counted);
    s.require("homomorphism: forgetting multiplicities, (s, c) ↦ s, preserves ⊕, ⊗, 0 and 1", !v, verdict(v));
  }
  {
    const auto v = mf::violated_morphism_law(mf::EmbedCount{}, {Integer(0), Integer(1), Integer(4), big});
    s.require("homomorphism: n ↦ (0, n) embeds Counting in the tropical semiring with multiplicities", !v, verdict(v));
  }
  {
    const auto v = mf::violated_morphism_law(mf::Negate<>{}, {ninf, -7, 0, 3, 12});
    const bool inverse = mf::Negate<>{}(mf::Negate<>{}(5)) == 5 && mf::Negate<>{}(ninf) == pinf;
    s.require("homomorphism: x ↦ −x is an isomorphism (max, +) → (min, +), −∞ ↦ +∞", !v && inverse, verdict(v));
  }
  {
    const auto v = mf::violated_morphism_law(mf::Finite<>{}, {ninf, -7, 0, 3});
    s.require("homomorphism: x ↦ [x ≠ −∞] maps (max, +) to Boolean", !v, verdict(v));
  }
  {
    const auto v = mf::violated_morphism_law(mf::Support{}, {Integer(0), Integer(1), Integer(9), big});
    s.require("homomorphism: n ↦ [n ≠ 0] maps Counting over ℕ to Boolean", !v, verdict(v));
  }
  {
    using P = mf::Product<mf::MaxPlus<>, mf::Counting>;
    const std::vector<P::value_type> samples = {P::zero(), P::one(), {3, Integer(2)}, {ninf, Integer(5)}, {-2, big}};
    const auto v1 = mf::violated_morphism_law(mf::First<mf::MaxPlus<>, mf::Counting>{}, samples);
    const auto v2 = mf::violated_morphism_law(mf::Second<mf::MaxPlus<>, mf::Counting>{}, samples);
    s.require("homomorphism: both projections of a product", !v1 && !v2, verdict(v1 ? v1 : v2));
  }
  {
    const auto v = mf::violated_morphism_law(CountOnly{}, counted);
    s.require("not a homomorphism: (s, c) ↦ c fails additivity, as the header says ((1, 1) ⊕ (0, 1) = (1, 1))", v.has_value() && *v == "h(a ⊕ b) ≠ h(a) ⊕ h(b)", verdict(v));
  }

  // ---- exactness ---------------------------------------------------------------------------------------
  {
    Integer x = 1;
    for (int k = 0; k < 100; k++) x = mf::Counting::times(x, Integer(2));
    s.equal("exact: 2^100 in Counting is 1267650600228229401496703205376", x.str(), std::string("1267650600228229401496703205376"));
    bool refused = false;
    try {
      mf::MaxPlus<>::times(std::numeric_limits<std::int64_t>::max() - 1, 5);
    } catch (const std::overflow_error&) { refused = true; }
    bool refused_low = false;
    try {
      mf::MinPlus<>::times(std::numeric_limits<std::int64_t>::min() + 2, -5);
    } catch (const std::overflow_error&) { refused_low = true; }
    s.require("exact: a tropical sum past the 64-bit range is refused, not wrapped", refused && refused_low);
    s.require("exact: −∞ absorbs in ⊗ and is the identity of ⊕", mf::MaxPlus<>::times(ninf, 5) == ninf && mf::MaxPlus<>::plus(ninf, -9) == -9);
  }
  return s.done();
}
