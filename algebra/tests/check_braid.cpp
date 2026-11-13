// Checks of morphe/algebra/braid.hpp and morphe/arithmetic/laurent.hpp.
//
// - The relations of Artin's presentation are the identity by handle reduction; σ_1 σ_2 σ_1⁻¹ and σ_1 σ_2 σ_2⁻¹…
//   are told apart; σ_1 σ_2 is not σ_2 σ_1.
// - Two methods agree. On words made trivial by inserting relators at random, and on random words, handle
//   reduction and Artin's faithful action on the free group give the same answer, for 3 to 6 strands.
// - Garside's element: Δ² is central, Δ σ_i Δ⁻¹ = σ_{n−i}, and (σ_1 ⋯ σ_{n−1})ⁿ = Δ².
// - What a word shows: permutations, exponent sums, the crossings of labelled strands (σ_1² winds once), strands
//   forgotten.
// - Dehornoy's order: 1 < σ_1, σ_2 < σ_1, σ_1⁻¹ < 1, and it is a total order on a sample.
// - Dynnikov's coordinates: the relations hold on random coordinates; the full twist Δ² fixes every curve system;
//   the growth rate of σ_1 σ_2⁻¹ is log((3 + √5)/2) and that of σ_1 σ_2 σ_3⁻¹ is the logarithm of the largest root
//   of x⁴ − 2x³ − 2x + 1 (Ko, Los and Song 2002, the least dilatation on four strands), and both agree with the
//   growth of the words of Artin's action.
// - Burau: the relations hold exactly, rows sum to 1, the determinant is (−t)^e; for three strands the matrix is
//   the identity exactly for the trivial braids of a sample; Alexander polynomials of the unknot, the Hopf link,
//   the trefoil, the figure-eight knot and the knot 5_1.

#include <morphe/algebra/braid.hpp>
#include <morphe/check.hpp>

#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace ma = morphe::algebra;
using ma::Braid;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Laurent;

namespace {

Braid garside(std::size_t n) {  // Δ_n = (σ_1 ⋯ σ_{n−1})(σ_1 ⋯ σ_{n−2}) ⋯ σ_1
  Braid d(n);
  for (std::size_t top = n - 1; top >= 1; top--)
    for (std::size_t i = 1; i <= top; i++) d.word.push_back(static_cast<int>(i));
  return d;
}

std::vector<int> random_word(std::mt19937_64& rng, std::size_t n, std::size_t length) {
  std::vector<int> w;
  for (std::size_t k = 0; k < length; k++) w.push_back((rng() % 2 ? 1 : -1) * static_cast<int>(1 + rng() % (n - 1)));
  return w;
}

// A word of the identity: relators and cancelling pairs inserted at random places.
std::vector<int> trivial_word(std::mt19937_64& rng, std::size_t n, std::size_t insertions) {
  std::vector<int> w;
  for (std::size_t k = 0; k < insertions; k++) {
    std::vector<int> piece;
    const int i = static_cast<int>(1 + rng() % (n - 1));
    const unsigned kind = rng() % 3;
    if (kind == 0) {
      piece = {i, -i};
    } else if (kind == 1 && i + 1 <= static_cast<int>(n) - 1) {
      piece = {i, i + 1, i, -(i + 1), -i, -(i + 1)};
    } else if (n >= 4) {
      int j = static_cast<int>(1 + rng() % (n - 1));
      if (std::abs(i - j) < 2) j = i >= 3 ? i - 2 : i + 2;
      if (j > static_cast<int>(n) - 1) j = i - 2;
      if (j < 1) continue;
      piece = {i, j, -i, -j};
    } else {
      piece = {-i, i};
    }
    if (rng() % 2) {  // a relator read backwards is a relator
      std::reverse(piece.begin(), piece.end());
      for (int& l : piece) l = -l;
    }
    const std::size_t at = w.empty() ? 0 : rng() % (w.size() + 1);
    w.insert(w.begin() + static_cast<std::ptrdiff_t>(at), piece.begin(), piece.end());
  }
  return w;
}

bool identity_images(const std::vector<std::vector<int>>& images) {
  for (std::size_t j = 0; j < images.size(); j++)
    if (images[j] != std::vector<int>{static_cast<int>(j) + 1}) return false;
  return true;
}

double growth_by_dynnikov(const Braid& b, int iterations) {
  ma::Dynnikov u = ma::generating_loops(b.strands);
  double before = u.norm().to_double(), rate = 0;
  for (int k = 0; k < iterations; k++) {
    u = ma::act(b, u);
    const double now = u.norm().to_double();
    rate = std::log(now / before), before = now;
  }
  return rate;
}

double growth_by_free_group(const Braid& b, int iterations) {
  double before = 1, rate = 0;
  for (int k = 1; k <= iterations; k++) {
    const auto images = ma::artin_images(ma::power(b, k), 50000000);
    if (!images) break;
    double longest = 0;
    for (const auto& w : *images) longest = std::max(longest, static_cast<double>(w.size()));
    rate = std::log(longest / before), before = longest;
  }
  return rate;
}

}  // namespace

int main() {
  morphe::Suite s;

  // ---- Laurent polynomials ----
  {
    const Laurent t = Laurent::t(), one(1);
    s.equal("(1 - t)(1 + t) is 1 - t^2", ((one - t) * (one + t)).str(), std::string("1 - t^2"));
    s.equal("t t^-1 is 1", (t * Laurent::t(-1)).str(), std::string("1"));
    s.equal("a polynomial up to units: -t^-1 + 3 - t is 1 - 3t + t^2", (-Laurent::t(-1) + Laurent(3) - t).normalised().str(), std::string("1 - 3t + t^2"));
    s.require("p - p is 0, and 0 has no terms", (t - t).is_zero() && (t - t) == Laurent());
  }

  // ---- the relations, and words told apart ----
  const auto is_trivial = [&](std::size_t n, std::vector<int> w) { return ma::trivial(Braid(n, std::move(w))).value(); };
  s.require("far generators commute: s1 s3 s1^-1 s3^-1 is 1 in B_4", is_trivial(4, {1, 3, -1, -3}));
  s.require("the braid relation: s1 s2 s1 s2^-1 s1^-1 s2^-1 is 1", is_trivial(3, {1, 2, 1, -2, -1, -2}));
  s.require("s1 s2 s1^-1 is not 1, nor is s1 s2 s1^-1 s2^-1", !is_trivial(3, {1, 2, -1}) && !is_trivial(3, {1, 2, -1, -2}));
  s.require("s1 s2 is not s2 s1, and s1 s2 s1 is s2 s1 s2", !ma::equal(Braid(3, {1, 2}), Braid(3, {2, 1})).value() && ma::equal(Braid(3, {1, 2, 1}), Braid(3, {2, 1, 2})).value());
  s.require("a letter that is not a generator is refused", [] {
    try {
      Braid(3, {3});
      return false;
    } catch (const std::invalid_argument&) {
      return true;
    }
  }());

  // ---- two methods agree ----
  {
    std::mt19937_64 rng(20261112);
    std::size_t trivial_found = 0, trivial_total = 0, disagreements = 0, compared = 0, nontrivial = 0;
    for (std::size_t n = 3; n <= 6; n++) {
      for (int k = 0; k < 120; k++) {
        const Braid b(n, trivial_word(rng, n, 2 + rng() % 7));
        trivial_total++;
        const auto by_handles = ma::trivial(b);
        if (by_handles && *by_handles) trivial_found++;
        if (const auto images = ma::artin_images(b)) {
          compared++;
          if (!by_handles || *by_handles != identity_images(*images)) disagreements++;
        }
      }
      for (int k = 0; k < 150; k++) {
        const Braid b(n, random_word(rng, n, 1 + rng() % 12));
        const auto by_handles = ma::trivial(b);
        const auto images = ma::artin_images(b);
        if (!images) continue;
        compared++;
        if (!identity_images(*images)) nontrivial++;
        if (!by_handles || *by_handles != identity_images(*images)) disagreements++;
      }
    }
    s.equal("480 words made of relators are found trivial by handle reduction", trivial_found, trivial_total);
    s.require("handle reduction and Artin's action agree on every word compared", disagreements == 0 && compared > 1000 && nontrivial > 400,
              std::to_string(compared) + " words, " + std::to_string(nontrivial) + " of them not trivial, " + std::to_string(disagreements) + " disagreements");
  }

  // ---- Garside's element ----
  for (std::size_t n = 3; n <= 5; n++) {
    const Braid delta = garside(n), twist = delta * delta;
    bool central = true, flips = true;
    for (std::size_t i = 1; i < n; i++) {
      const Braid g(n, {static_cast<int>(i)});
      central = central && ma::equal(twist * g, g * twist).value();
      flips = flips && ma::equal(delta * g * ma::inverse(delta), Braid(n, {static_cast<int>(n - i)})).value();
    }
    Braid cycle(n);
    for (std::size_t i = 1; i < n; i++) cycle.word.push_back(static_cast<int>(i));
    s.require("B_" + std::to_string(n) + ": the full twist is central, the half twist flips the generators, (s1...s" + std::to_string(n - 1) + ")^" + std::to_string(n) + " is the full twist",
              central && flips && ma::equal(ma::power(cycle, static_cast<long>(n)), twist).value() && ma::pure(twist) && !ma::pure(delta));
  }

  // ---- what a word shows ----
  {
    const Braid b(4, {1, 2, -3, 1});
    // Positions after each letter hold the strands (1,0,2,3), (1,2,0,3), (1,2,3,0), (2,1,3,0): strand 0 ends at 3.
    s.equal("the permutation of s1 s2 s3^-1 s1", ma::permutation(b), std::vector<std::size_t>{3, 1, 0, 2});
    s.equal("its exponent sum", ma::exponent_sum(b), 2L);
    const Braid turn(2, {1, 1});
    s.require("s1^2 is pure and its two strands cross twice: they wind once", ma::pure(turn) && ma::crossings(turn)[0][1] == 2 && ma::crossings(turn)[1][0] == 2);
    const Braid three(3, {1, 1, 2, 2});  // the middle strand winds once about each of the others
    const auto c = ma::crossings(three);
    s.require("s1^2 s2^2: strands 0 and 1 wind once, 1 and 2 once, 0 and 2 not at all", c[0][1] == 2 && c[1][2] == 2 && c[0][2] == 0);
    s.equal("forgetting the middle strand of s1^2 s2^2 leaves the identity of B_2", ma::forget(three, 1).word, std::vector<int>{});
    s.equal("forgetting the last strand leaves s1^2", ma::forget(three, 2).word, std::vector<int>{1, 1});
    s.require("forgetting commutes with the product on pure braids", ma::equal(ma::forget(three * three, 0), ma::forget(three, 0) * ma::forget(three, 0)).value());
  }

  // ---- Dehornoy's order ----
  {
    const Braid one(3), s1(3, {1}), s2(3, {2});
    s.require("1 < s1, s2 < s1, s1^-1 < 1, s1 = s1", ma::compare(one, s1).value() == -1 && ma::compare(s2, s1).value() == -1 && ma::compare(ma::inverse(s1), one).value() == -1 && ma::compare(s1, s1).value() == 0);
    std::mt19937_64 rng(5);
    bool antisymmetric = true, invariant = true;
    for (int k = 0; k < 200; k++) {
      const Braid a(4, random_word(rng, 4, 1 + rng() % 8)), b(4, random_word(rng, 4, 1 + rng() % 8)), c(4, random_word(rng, 4, 1 + rng() % 5));
      const int ab = ma::compare(a, b).value();
      antisymmetric = antisymmetric && ab == -ma::compare(b, a).value();
      invariant = invariant && ab == ma::compare(c * a, c * b).value();  // the order is invariant on the left
    }
    s.require("on 200 pairs the order is antisymmetric and invariant under multiplication on the left", antisymmetric && invariant);
  }

  // ---- Dynnikov's coordinates ----
  {
    std::mt19937_64 rng(11);
    bool relations = true, twist_fixes = true;
    for (std::size_t n = 3; n <= 7; n++) {
      for (int k = 0; k < 200; k++) {
        ma::Dynnikov u;
        for (std::size_t i = 0; i + 2 < n; i++) u.a.push_back(Integer(static_cast<long long>(rng() % 41) - 20)), u.b.push_back(Integer(static_cast<long long>(rng() % 41) - 20));
        const int i = static_cast<int>(1 + rng() % (n - 1));
        relations = relations && ma::act(Braid(n, {i, -i}), u) == u && ma::act(Braid(n, {-i, i}), u) == u;
        if (i + 1 <= static_cast<int>(n) - 1) relations = relations && ma::act(Braid(n, {i, i + 1, i}), u) == ma::act(Braid(n, {i + 1, i, i + 1}), u);
        for (int j = 1; j <= static_cast<int>(n) - 1; j++)
          if (std::abs(i - j) >= 2) relations = relations && ma::act(Braid(n, {i, j}), u) == ma::act(Braid(n, {j, i}), u);
        if (k < 20) twist_fixes = twist_fixes && ma::act(garside(n) * garside(n), u) == u;
      }
    }
    s.require("the relations of the braid group hold on 1000 random curve systems, 3 to 7 strands", relations);
    s.require("the full twist fixes every curve system", twist_fixes);
    const double golden = std::log((3 + std::sqrt(5.0)) / 2);
    s.check("the growth rate of s1 s2^-1 by Dynnikov's coordinates against log((3 + sqrt 5) / 2)", growth_by_dynnikov(Braid(3, {1, -2}), 60) - golden, 1e-12);
    s.check("the same by the words of Artin's action", growth_by_free_group(Braid(3, {1, -2}), 14) - golden, 1e-3);
    double lo = 2, hi = 3;  // the largest root of x^4 - 2x^3 - 2x + 1
    for (int k = 0; k < 80; k++) {
      const double x = 0.5 * (lo + hi);
      (x * x * x * x - 2 * x * x * x - 2 * x + 1 < 0 ? lo : hi) = x;
    }
    s.check("the growth rate of s1 s2 s3^-1 against the largest root of x^4 - 2x^3 - 2x + 1", growth_by_dynnikov(Braid(4, {1, 2, -3}), 80) - std::log(lo), 1e-12);
    s.check("the same by the words of Artin's action", growth_by_free_group(Braid(4, {1, 2, -3}), 14) - std::log(lo), 1e-2);
    s.check("a braid of finite order does not grow: s1 s2 on three strands", growth_by_dynnikov(Braid(3, {1, 2}), 60), 1e-12);
  }

  // ---- Burau and Alexander ----
  {
    const auto B = [](std::size_t n, std::vector<int> w) { return ma::burau(Braid(n, std::move(w))); };
    s.require("Burau keeps the relations exactly", B(3, {1, 2, 1}) == B(3, {2, 1, 2}) && B(4, {1, 3}) == B(4, {3, 1}) && B(3, {1, -1}) == ma::laurent_identity(3) && B(3, {-2, 2}) == ma::laurent_identity(3));
    std::mt19937_64 rng(3);
    bool rows = true, determinants = true, faithful3 = true;
    for (int k = 0; k < 120; k++) {
      const std::size_t n = 2 + rng() % 4;
      const Braid b(n, random_word(rng, n, 1 + rng() % 9));
      const auto m = ma::burau(b);
      for (const auto& row : m) {
        Laurent sum;
        for (const auto& x : row) sum = sum + x;
        rows = rows && sum == Laurent(1);
      }
      const long e = ma::exponent_sum(b);
      determinants = determinants && ma::determinant(m) == Laurent::monomial(Integer(e % 2 == 0 ? 1 : -1), e);
      const Braid c(3, k % 3 == 0 ? trivial_word(rng, 3, 4) : random_word(rng, 3, 1 + rng() % 10));
      faithful3 = faithful3 && (ma::burau(c) == ma::laurent_identity(3)) == ma::trivial(c).value();
    }
    s.require("its rows sum to 1 and its determinant is (-t)^e on 120 random braids", rows && determinants);
    s.require("on three strands the matrix is the identity exactly for the trivial braids of a sample", faithful3);
    s.equal("the Alexander polynomial of the unknot, closure of s1", ma::alexander(Braid(2, {1})).str(), std::string("1"));
    s.equal("of the Hopf link, closure of s1^2", ma::alexander(Braid(2, {1, 1})).str(), std::string("1 - t"));
    s.equal("of the trefoil, closure of s1^3", ma::alexander(Braid(2, {1, 1, 1})).str(), std::string("1 - t + t^2"));
    s.equal("of the figure-eight knot, closure of (s1 s2^-1)^2", ma::alexander(Braid(3, {1, -2, 1, -2})).str(), std::string("1 - 3t + t^2"));
    s.equal("of the knot 5_1, closure of s1^5", ma::alexander(Braid(2, {1, 1, 1, 1, 1})).str(), std::string("1 - t + t^2 - t^3 + t^4"));
    s.equal("the trefoil as a closure on three strands, s1 s2 s1 s2", ma::alexander(Braid(3, {1, 2, 1, 2})).str(), std::string("1 - t + t^2"));
  }

  return s.done();
}
