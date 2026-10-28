// Checks of morphe/arithmetic/cyclotomic.hpp: the cyclotomic polynomials (degree φ(m), monic,
// Π_{d | m} Φ_d = x^m − 1, Φ_p and Φ_{2m}(x) = Φ_m(−x), Φ_m(1), and 105 the least order with a
// coefficient outside {−1, 0, 1}); the field axioms of ℚ(ζ_m) on random numbers, inverses
// among them; roots of unity; the Galois group (σ_k ring homomorphisms, σ_k σ_l = σ_{kl},
// conjugation, norms and traces rational); lifting between fields and mixed arithmetic;
// quadratic subfields found exactly (√2, √−3, √5, Gauss sums with g² = (−1)^{(p−1)/2} p); and the
// complex values used only for display.

#include <morphe/arithmetic/cyclotomic.hpp>
#include <morphe/check.hpp>

#include <complex>
#include <numbers>
#include <random>
#include <string>
#include <vector>

namespace ma = morphe::arithmetic;
using ma::Cyclotomic;
using ma::Integer;
using ma::Rational;

namespace {

// Integer polynomials, low degree first.
using IPoly = std::vector<Integer>;
IPoly mul(const IPoly& a, const IPoly& b) {
  IPoly p(a.size() + b.size() - 1, Integer(0));
  for (std::size_t i = 0; i < a.size(); i++)
    for (std::size_t j = 0; j < b.size(); j++) p[i + j] += a[i] * b[j];
  return p;
}

Cyclotomic random_number(std::mt19937_64& g, unsigned m) {
  std::vector<Rational> c(ma::euler_phi(m));
  for (auto& x : c) x = Rational(Integer(static_cast<long long>(g() % 11) - 5), Integer(static_cast<long long>(1 + g() % 4)));
  return Cyclotomic::from_basis(m, c);
}

// The Legendre symbol (a/p) by Euler's criterion.
int legendre(long long a, long long p) {
  long long r = 1, b = ((a % p) + p) % p, e = (p - 1) / 2;
  if (b == 0) return 0;
  while (e) {
    if (e & 1) r = r * b % p;
    b = b * b % p;
    e >>= 1;
  }
  return r == 1 ? 1 : -1;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(23);

  // ---- the cyclotomic polynomials ------------------------------------------------------------------------------------------
  {
    int bad_degree = 0, bad_product = 0, bad_prime = 0, bad_even = 0, bad_at_one = 0;
    unsigned first_large = 0;
    for (unsigned m = 1; m <= 120; m++) {
      const IPoly P = ma::cyclotomic_polynomial(m);
      if (P.size() != ma::euler_phi(m) + 1 || P.back() != Integer(1)) bad_degree++;
      // Π_{d | m} Φ_d = x^m − 1.
      IPoly prod{Integer(1)};
      for (unsigned d = 1; d <= m; d++)
        if (m % d == 0) prod = mul(prod, ma::cyclotomic_polynomial(d));
      IPoly want(m + 1, Integer(0));
      want[0] = -1, want[m] = 1;
      if (prod != want) bad_product++;
      // Φ_m(1): p for m a power of a prime p, 1 for any other m > 1.
      if (m > 1) {
        Integer at_one = 0;
        for (const auto& c : P) at_one += c;
        unsigned q = 2, x = m;
        while (x % q) q++;
        while (x % q == 0) x /= q;
        if (at_one != Integer(x == 1 ? static_cast<long long>(q) : 1)) bad_at_one++;
      }
      for (const auto& c : P)
        if (first_large == 0 && (c > Integer(1) || c < Integer(-1))) first_large = m;
    }
    for (unsigned p : {2u, 3u, 5u, 7u, 11u, 13u, 31u, 97u}) {
      const IPoly P = ma::cyclotomic_polynomial(p);
      if (P != IPoly(p, Integer(1))) bad_prime++;
    }
    for (unsigned m : {3u, 5u, 9u, 15u, 21u, 45u}) {
      IPoly P = ma::cyclotomic_polynomial(m);
      for (std::size_t i = 1; i < P.size(); i += 2) P[i] = -P[i];
      if (P != ma::cyclotomic_polynomial(2 * m)) bad_even++;
    }
    s.equal("cyclotomic: Φ_m is monic of degree φ(m), m ≤ 120 (failures)", bad_degree, 0);
    s.equal("cyclotomic: Π_{d | m} Φ_d = x^m − 1, m ≤ 120 (failures)", bad_product, 0);
    s.equal("cyclotomic: Φ_m(1) = p for m = p^k and 1 for other m > 1 (failures)", bad_at_one, 0);
    s.equal("cyclotomic: Φ_p = 1 + x + … + x^{p−1} for eight primes (failures)", bad_prime, 0);
    s.equal("cyclotomic: Φ_{2m}(x) = Φ_m(−x) for odd m > 1 (failures)", bad_even, 0);
    s.equal("cyclotomic: 105 is the least m whose Φ_m has a coefficient outside {−1, 0, 1}", first_large, 105u);
    const IPoly P105 = ma::cyclotomic_polynomial(105);
    s.require("cyclotomic: … and that coefficient is −2, at x^7 and x^41", P105[7] == Integer(-2) && P105[41] == Integer(-2));
  }
  // ---- the field axioms ---------------------------------------------------------------------------------------------------
  {
    const std::vector<unsigned> orders = {1, 2, 3, 4, 5, 7, 8, 9, 12, 15, 16, 20, 24, 30};
    int ring = 0, inverse = 0, division = 0;
    for (unsigned m : orders)
      for (int t = 0; t < 12; t++) {
        const Cyclotomic a = random_number(g, m), b = random_number(g, m), c = random_number(g, m);
        if (!((a + b) + c == a + (b + c)) || !((a * b) * c == a * (b * c)) || !(a * (b + c) == a * b + a * c) || !(a * b == b * a) ||
            !(a + b == b + a) || !(a - a == Cyclotomic(0)) || !((a - b) + b == a))
          ring++;
        if (!a.is_zero() && !(a * a.inverse() == Cyclotomic(1))) inverse++;
        if (!b.is_zero() && !((a / b) * b == a)) division++;
      }
    s.equal("cyclotomic: associativity, commutativity, distributivity on 168 random triples in 14 fields (failures)", ring, 0);
    s.equal("cyclotomic: a a⁻¹ = 1 by the extended Euclidean algorithm modulo Φ_m (failures)", inverse, 0);
    s.equal("cyclotomic: (a / b) b = a (failures)", division, 0);
  }
  // ---- roots of unity -----------------------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (unsigned m = 1; m <= 40; m++) {
      const Cyclotomic z = Cyclotomic::zeta(m);
      if (!(z.power(m) == Cyclotomic(1))) bad++;
      for (unsigned k = 1; k < m; k++)
        if (z.power(k) == Cyclotomic(1)) bad++;
      Cyclotomic sum = 0;
      for (unsigned k = 0; k < m; k++) sum += Cyclotomic::zeta(m, k);
      if (!(sum == Cyclotomic(m == 1 ? 1 : 0))) bad++;
      if (!(Cyclotomic::zeta(m, -1) == z.inverse()) || !(Cyclotomic::zeta(m, 3) == z * z * z)) bad++;
    }
    s.equal("cyclotomic: ζ_m has order exactly m, Σ_k ζ_m^k = 0, ζ^{−1} = ζ⁻¹, for m ≤ 40 (failures)", bad, 0);
    s.require("cyclotomic: ζ_2 = −1 and ζ_4² = −1, held as rationals", Cyclotomic::zeta(2) == Cyclotomic(-1) && Cyclotomic::zeta(2).is_rational() &&
                                                                       (Cyclotomic::zeta(4) * Cyclotomic::zeta(4)).is_rational());
  }
  // ---- the Galois group ---------------------------------------------------------------------------------------------------
  {
    int hom = 0, compose = 0, conj = 0, norms = 0;
    for (unsigned m : {5u, 7u, 8u, 9u, 12u, 15u, 20u}) {
      std::vector<long long> units;
      for (long long k = 1; k < m; k++)
        if (std::gcd(k, static_cast<long long>(m)) == 1) units.push_back(k);
      for (int t = 0; t < 6; t++) {
        const Cyclotomic a = random_number(g, m), b = random_number(g, m);
        for (long long k : units) {
          if (!(a.galois(k) * b.galois(k) == (a * b).galois(k)) || !(a.galois(k) + b.galois(k) == (a + b).galois(k))) hom++;
          for (long long l : units)
            if (!(a.galois(k).galois(l) == a.galois(k * l))) compose++;
        }
        if (!(a.conj() == a.galois(-1)) || !(a.conj().conj() == a) || !((a * a.conj()).conj() == a * a.conj())) conj++;
        Cyclotomic norm = 1, trace = 0;
        for (long long k : units) norm *= a.galois(k), trace += a.galois(k);
        if (!norm.is_rational() || !trace.is_rational() || (!a.is_zero() && norm.is_zero())) norms++;
      }
    }
    s.equal("cyclotomic: each σ_k is a ring homomorphism (failures)", hom, 0);
    s.equal("cyclotomic: σ_l σ_k = σ_{kl} (failures)", compose, 0);
    s.equal("cyclotomic: conjugation is σ_{−1}, an involution, and a ā is real (failures)", conj, 0);
    s.equal("cyclotomic: norms and traces down to ℚ are rational, and non-zero norms of non-zero numbers (failures)", norms, 0);
    bool threw = false;
    try {
      (void)Cyclotomic::zeta(12).galois(2);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    s.require("cyclotomic: σ_2 of ℚ(ζ_12) is refused, gcd(2, 12) ≠ 1", threw);
  }
  // ---- lifting and mixed fields -------------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (auto [m, M] : std::vector<std::pair<unsigned, unsigned>>{{3, 12}, {4, 12}, {5, 15}, {6, 30}, {8, 24}, {9, 36}}) {
      for (int t = 0; t < 8; t++) {
        const Cyclotomic a = random_number(g, m), b = random_number(g, m);
        if (!(a.lift(M) + b.lift(M) == (a + b).lift(M)) || !(a.lift(M) * b.lift(M) == (a * b).lift(M)) || !(a.lift(M) == a)) bad++;
      }
      if (!(Cyclotomic::zeta(m).lift(M) == Cyclotomic::zeta(M, M / m))) bad++;
      if (Cyclotomic::zeta(m).lift(M).field() != M) bad++;
    }
    s.equal("cyclotomic: lifting to ℚ(ζ_M) is a homomorphism, ζ_m ↦ ζ_M^{M/m} (failures)", bad, 0);
    const Cyclotomic mixed = Cyclotomic::zeta(3) + Cyclotomic::zeta(4);
    s.require("cyclotomic: ζ_3 + ζ_4 is held in ℚ(ζ_12) and equals ζ_12^4 + ζ_12^3",
              mixed.field() == 12 && mixed == Cyclotomic::zeta(12, 4) + Cyclotomic::zeta(12, 3));
    s.require("cyclotomic: ζ_3 held in ℚ(ζ_3) and in ℚ(ζ_6) are equal; ζ_6 = −ζ_3²",
              Cyclotomic::zeta(3) == Cyclotomic::zeta(6, 2) && Cyclotomic::zeta(6) == -Cyclotomic::zeta(3, 2));
    const auto c = Cyclotomic::zeta(4).coefficients_in(12);
    s.require("cyclotomic: the canonical form of i in ℚ(ζ_12) has φ(12) = 4 coordinates", c.size() == 4);
  }
  // ---- quadratic subfields --------------------------------------------------------------------------------------------------
  {
    const Cyclotomic r2 = Cyclotomic::zeta(8) + Cyclotomic::zeta(8, -1);
    s.require("cyclotomic: (ζ_8 + ζ_8⁻¹)² = 2", r2 * r2 == Cyclotomic(2) && !r2.is_rational());
    const Cyclotomic r3 = Cyclotomic::zeta(3) - Cyclotomic::zeta(3, 2);
    s.require("cyclotomic: (ζ_3 − ζ_3²)² = −3", r3 * r3 == Cyclotomic(-3));
    const Cyclotomic phi = Cyclotomic(1) + Cyclotomic(2) * (Cyclotomic::zeta(5) + Cyclotomic::zeta(5, 4));  // 1 + 2 cos(2π/5)·2 = √5
    s.require("cyclotomic: (1 + 2(ζ_5 + ζ_5⁴))² = 5", phi * phi == Cyclotomic(5));
    // Gauss sums g = Σ_a (a/p) ζ_p^a: g² = (−1)^{(p−1)/2} p (Ireland and Rosen 1990, ch. 6).
    int bad = 0;
    for (long long p : {3, 5, 7, 11, 13, 17, 19, 23, 29, 31}) {
      std::vector<Rational> c(static_cast<std::size_t>(p));
      for (long long a = 0; a < p; a++) c[static_cast<std::size_t>(a)] = legendre(a, p);
      const Cyclotomic gs = Cyclotomic::from_powers(static_cast<unsigned>(p), c);
      if (!(gs * gs == Cyclotomic(((p - 1) / 2) % 2 ? -p : p))) bad++;
    }
    s.equal("cyclotomic: Gauss sums square to (−1)^{(p−1)/2} p for ten odd primes (failures)", bad, 0);
  }
  // ---- complex values, for display only -------------------------------------------------------------------------------------
  {
    double worst = 0;
    for (unsigned m : {5u, 12u, 30u})
      for (long long k = 0; k < m; k++) {
        const std::complex<double> want = std::polar(1.0, 2 * std::numbers::pi * static_cast<double>(k) / m);
        worst = std::max(worst, std::abs(Cyclotomic::zeta(m, k).to_complex() - want));
      }
    s.check("cyclotomic: to_complex(ζ_m^k) against exp(2πik/m), for display (largest deviation)", worst, 1e-12);
    s.equal("cyclotomic: str() of 1 − 2ζ_12³", (Cyclotomic(1) - Cyclotomic(2) * Cyclotomic::zeta(12, 3)).str(), std::string("1 - 2*ζ12^3"));
  }
  return s.done();
}
