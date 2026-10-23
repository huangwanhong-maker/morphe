// Checks of morphe/arithmetic: integers against 128-bit machine arithmetic where it does not
// overflow; division, beyond it, by the identity a = q b + r with |r| < |b|, on operands
// chosen to exercise every branch of Knuth's Algorithm D; known values (2^100, 25!, Catalan
// numbers); gcd identities; rationals by telescoping and harmonic sums; residues and primes.

#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/modular.hpp>
#include <morphe/arithmetic/rational.hpp>
#include <morphe/check.hpp>

#include <cstdint>
#include <random>
#include <string>

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;
namespace ma = morphe::arithmetic;

namespace {

std::string str128(__int128 v) {
  if (v == 0) return "0";
  const bool neg = v < 0;
  unsigned __int128 m = neg ? static_cast<unsigned __int128>(-(v + 1)) + 1 : static_cast<unsigned __int128>(v);
  std::string s;
  while (m) { s.insert(s.begin(), static_cast<char>('0' + static_cast<int>(m % 10))); m /= 10; }
  return neg ? "-" + s : s;
}

// A random integer of up to `limbs` limbs, with runs of all-ones and all-zeros limbs, which
// are the hard cases of the trial quotient.
Integer random_integer(std::mt19937_64& g, int limbs) {
  Integer::Limbs m(1 + g() % limbs);
  for (auto& x : m) {
    const int kind = g() % 4;
    x = kind == 0 ? 0xffffffffu : kind == 1 ? 0u : static_cast<std::uint32_t>(g());
  }
  return Integer::make(g() & 1, m);
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(20261023);

  // ---- against __int128 ----------------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (int k = 0; k < 20000; k++) {
      const long long a = static_cast<long long>(g()) >> (g() % 63), b = static_cast<long long>(g()) >> (g() % 63);
      const __int128 A = a, B = b;
      if ((Integer(a) + Integer(b)).str() != str128(A + B)) bad++;
      if ((Integer(a) - Integer(b)).str() != str128(A - B)) bad++;
      if ((Integer(a) * Integer(b)).str() != str128(A * B)) bad++;
      if (b != 0 && (Integer(a) / Integer(b)).str() != str128(A / B)) bad++;
      if (b != 0 && (Integer(a) % Integer(b)).str() != str128(A % B)) bad++;
      if ((Integer(a) < Integer(b)) != (a < b)) bad++;
    }
    s.equal("arithmetic: +, -, *, /, %, < of 20000 random pairs of words agree with 128-bit arithmetic (disagreements)", bad, 0);
    bad = 0;
    for (int k = 0; k < 20000; k++) {
      const __int128 A = (static_cast<__int128>(static_cast<long long>(g())) << (g() % 60)) ^ static_cast<long long>(g());
      const __int128 B = static_cast<long long>(g()) >> (g() % 63);
      const Integer a = Integer::parse(str128(A)), b = Integer::parse(str128(B));
      if (a.str() != str128(A)) bad++;
      if ((a + b).str() != str128(A + B) || (a - b).str() != str128(A - B)) bad++;
      if (B != 0 && ((a / b).str() != str128(A / B) || (a % b).str() != str128(A % B))) bad++;
    }
    s.equal("arithmetic: parsing, +, -, /, % of integers up to 2^123 agree with 128-bit arithmetic (disagreements)", bad, 0);
    s.equal("arithmetic: the most negative word is negated into limbs", (-Integer(std::numeric_limits<long long>::min())).str(), std::string("9223372036854775808"));
    s.equal("arithmetic: and back into a word", (-(-Integer(std::numeric_limits<long long>::min()))).fits_int64(), true);
  }
  // ---- division beyond 128 bits -------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (int k = 0; k < 5000; k++) {
      const Integer a = random_integer(g, 12), b = random_integer(g, 7);
      if (b.is_zero()) continue;
      Integer q, r;
      divmod(a, b, q, r);
      if (q * b + r != a || !(ma::abs(r) < ma::abs(b)) || (r.sign() != 0 && r.sign() != a.sign())) bad++;
      if (ma::mod(a, b).sign() < 0 || ma::floor_div(a, b) * b + ma::mod(a, b) != a) bad++;
      if ((a * b) / b != a || (a * b) % b != 0) bad++;
      if ((a + b) - b != a) bad++;
    }
    s.equal("arithmetic: a = q b + r, |r| < |b|, sign of r, Euclidean pair and (a b) / b = a, on 5000 pairs of up to 12 and 7 limbs (failures)", bad, 0);
  }
  // ---- known values ----------------------------------------------------------------------------------------------------
  s.equal("arithmetic: 2^100", ma::pow(Integer(2), 100).str(), std::string("1267650600228229401496703205376"));
  s.equal("arithmetic: 25!", ma::factorial(25).str(), std::string("15511210043330985984000000"));
  s.equal("arithmetic: 15!! = 2027025", ma::double_factorial(15), Integer(2027025));
  s.equal("arithmetic: C(100, 50)", ma::binomial(100, 50).str(), std::string("100891344545564193334812497256"));
  {
    // The Catalan numbers by C(2n, n) / (n + 1) and by their recurrence, to n = 60.
    std::vector<Integer> c{1};
    for (int n = 1; n <= 60; n++) {
      Integer sum = 0;
      for (int i = 0; i < n; i++) sum += c[i] * c[n - 1 - i];
      c.push_back(sum);
    }
    s.equal("arithmetic: Catalan numbers by the recurrence and by C(2n, n) / (n + 1), n = 60", c[60], ma::exact_div(ma::binomial(120, 60), Integer(61)));
  }
  s.equal("arithmetic: the decimal text of 10^40 - 1 reads back", Integer::parse(std::string(40, '9')) + 1, ma::pow(Integer(10), 40));
  // ---- gcd -------------------------------------------------------------------------------------------------------------
  {
    // gcd(2^a - 1, 2^b - 1) = 2^gcd(a, b) - 1; consecutive Fibonacci numbers are coprime.
    bool ok = true;
    for (int a = 1; a <= 90; a += 7)
      for (int b = 1; b <= 120; b += 11) {
        const Integer x = ma::pow(Integer(2), a) - 1, y = ma::pow(Integer(2), b) - 1;
        if (ma::gcd(x, y) != ma::pow(Integer(2), std::gcd(a, b)) - 1) ok = false;
        const auto bz = ma::extended_gcd(x, -y);
        if (bz.s * x + bz.t * (-y) != bz.g || bz.g != ma::gcd(x, y)) ok = false;
      }
    s.require("arithmetic: gcd(2^a - 1, 2^b - 1) = 2^gcd(a, b) - 1, and Bezout's identity", ok);
    Integer f0 = 0, f1 = 1;
    for (int k = 0; k < 300; k++) { Integer f2 = f0 + f1; f0 = f1; f1 = f2; }
    s.require("arithmetic: consecutive Fibonacci numbers F_300, F_301 are coprime", ma::gcd(f0, f1).is_one(), f1.str().substr(0, 12) + "...");
  }
  // ---- rationals -------------------------------------------------------------------------------------------------------
  {
    Rational t = 0, h = 0;
    for (int k = 1; k <= 200; k++) t += Rational(1, Integer(k) * Integer(k + 1));
    for (int k = 1; k <= 10; k++) h += Rational(1, k);
    s.equal("arithmetic: sum of 1 / (k (k + 1)), k <= 200, telescopes to 200/201", t, Rational(200, 201));
    s.equal("arithmetic: the harmonic number H_10 = 7381/2520", h, Rational::parse("7381/2520"));
    s.equal("arithmetic: rationals in lowest terms with positive denominator", Rational(Integer(6), Integer(-4)).str(), std::string("-3/2"));
    Rational big = Rational(ma::pow(Integer(3), 200), ma::pow(Integer(2), 300));
    s.check("arithmetic: 3^200 / 2^300 to double (relative error)", big.to_double() / std::exp(200 * std::log(3.0) - 300 * std::log(2.0)) - 1, 1e-13);
    s.require("arithmetic: order of rationals", Rational(1, 3) < Rational(1, 2) && Rational(-1, 2) < Rational(-1, 3));
  }
  // ---- residues --------------------------------------------------------------------------------------------------------
  {
    const std::uint64_t p = 2305843009213693951ull;  // the Mersenne prime 2^61 - 1
    s.require("arithmetic: 2^61 - 1 is prime and 2^61 + 1 is not; Fermat's little theorem", ma::is_prime(p) && !ma::is_prime(p + 2) && ma::mod_pow(3, p - 1, p) == 1);
    int primes = 0;
    for (std::uint64_t n = 0; n < 100000; n++) primes += ma::is_prime(n);
    s.equal("arithmetic: pi(10^5) = 9592", primes, 9592);
    bool inv = true;
    for (std::uint64_t a = 1; a < 200; a++) inv = inv && ma::mod_mul(a, ma::mod_inverse(a, 1000003), 1000003) == 1;
    s.require("arithmetic: inverses modulo the prime 1000003", inv);
    s.equal("arithmetic: residue of -7 modulo 5", ma::residue(Integer(-7), 5), std::uint64_t(3));
  }
  return s.done();
}
