// Checks of morphe/algebra/characters.hpp, the Dixon–Schneider character table. On every group
// of a list (cyclic, elementary abelian, dihedral, quaternion, alternating, symmetric, the
// Frobenius group of order 21, SL(2, 3), PSL(2, 7), direct products, from order 1 to 720): both orthogonality
// relations exactly; Σ χ(1)² = |G|; as many characters as classes; χ(g⁻¹) = conj χ(g); linear
// characters are homomorphisms; the table is stable under the Galois group of ℚ(ζ_e); the class
// multiplication coefficients counted in the group equal Burnside's character formula; the
// regular character is Σ χ(1) χ; Frobenius–Schur indicators lie in {−1, 0, 1} and count square
// roots; as many real characters as real classes. The known tables by their invariants, computed
// from the definitions: C_n takes every n-th root of unity at a generator, V_4 is ±1, D_4 and Q_8
// share one table and are told apart by the indicator, A_4's linear characters are the cube
// roots of unity, A_5's take the golden ratio, the Frobenius group's and PSL(2, 7)'s (−1 ± √−7)/2,
// SL(2, 3) has a quaternionic character. The class matrices against the coefficients. S_4, S_5 and
// S_6 against the Murnaghan–Nakayama characters, and direct products against products of
// characters. The characteristic polynomial modulo p against determinants over ℤ.

#include <morphe/algebra/characters.hpp>
#include <morphe/algebra/determinant.hpp>
#include <morphe/algebra/symmetric_representations.hpp>
#include <morphe/check.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace ma = morphe::algebra;
using morphe::arithmetic::Cyclotomic;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

Cyclotomic C(long long v) { return Cyclotomic(v); }

struct Generic {
  int orthogonality = 0, degrees = 0, inverses = 0, linear = 0, galois = 0, burnside = 0, regular = 0, indicators = 0, real = 0;
};

// The properties every character table has, each counted as failures.
void generic(const ma::FiniteGroup& G, const ma::CharacterTable& T, Generic& f, std::string& where) {
  const std::size_t r = T.classes.size();
  const long long n = static_cast<long long>(G.order());
  const auto before = f.orthogonality + f.degrees + f.inverses + f.linear + f.galois + f.burnside + f.regular + f.indicators + f.real;
  if (T.size() != r) f.degrees++;
  // Row orthogonality Σ_j h_j χ(g_j) conj ψ(g_j) = |G| δ, and column orthogonality
  // Σ_χ χ(g_j) conj χ(g_k) = δ_jk |G| / h_j.
  for (std::size_t x = 0; x < T.size(); x++)
    for (std::size_t y = 0; y < T.size(); y++) {
      Cyclotomic sum = 0;
      for (std::size_t j = 0; j < r; j++) sum += C(static_cast<long long>(T.class_sizes[j])) * T(x, j) * T(y, j).conj();
      if (!(sum == C(x == y ? n : 0))) f.orthogonality++;
    }
  for (std::size_t j = 0; j < r; j++)
    for (std::size_t k = 0; k < r; k++) {
      Cyclotomic sum = 0;
      for (std::size_t x = 0; x < T.size(); x++) sum += T(x, j) * T(x, k).conj();
      if (!(sum == C(j == k ? n / static_cast<long long>(T.class_sizes[j]) : 0))) f.orthogonality++;
    }
  long long squares = 0;
  for (std::size_t x = 0; x < T.size(); x++) squares += static_cast<long long>(T.degree(x) * T.degree(x));
  if (squares != n) f.degrees++;
  for (std::size_t x = 0; x < T.size(); x++)
    for (std::size_t j = 0; j < r; j++)
      if (!(T(x, T.inverse_class[j]) == T(x, j).conj())) f.inverses++;
  // Linear characters are homomorphisms G → ℂ^×.
  if (n <= 200)
    for (std::size_t x = 0; x < T.size(); x++) {
      if (T.degree(x) != 1) continue;
      for (std::uint32_t a = 0; a < G.n; a++)
        for (std::uint32_t b = 0; b < G.n; b++)
          if (!(T.value(x, G.mul(a, b)) == T.value(x, a) * T.value(x, b))) f.linear++;
    }
  // Galois stability: σ_k ∘ χ is a row, and σ_k permutes the rows.
  for (long long k = 1; k <= static_cast<long long>(T.exponent); k++) {
    if (std::gcd(k, static_cast<long long>(T.exponent)) != 1) continue;
    std::vector<char> hit(T.size(), 0);
    try {
      for (std::size_t x = 0; x < T.size(); x++) hit[ma::galois_conjugate(T, x, k)] = 1;
    } catch (const std::logic_error&) {
      f.galois++;
    }
    if (std::count(hit.begin(), hit.end(), 1) != static_cast<long>(T.size())) f.galois++;
  }
  // a_{jkl} = (h_j h_k / |G|) Σ_χ χ(g_j) χ(g_k) conj χ(g_l) / χ(1).
  const auto a = ma::class_multiplication_coefficients(G, T);
  for (std::size_t j = 0; j < r; j++) {
    const auto M = ma::class_matrix(G, T, j);
    for (std::size_t x = 0; x < r * r; x++)
      if (M[x] != a[j * r * r + x]) f.burnside++;
  }
  for (std::size_t j = 0; j < r; j++)
    for (std::size_t k = 0; k < r; k++)
      for (std::size_t l = 0; l < r; l++) {
        Cyclotomic sum = 0;
        for (std::size_t x = 0; x < T.size(); x++)
          sum += T(x, j) * T(x, k) * T(x, l).conj() / C(static_cast<long long>(T.degree(x)));
        sum *= Cyclotomic(Rational(Integer(static_cast<long long>(T.class_sizes[j] * T.class_sizes[k])), Integer(n)));
        if (!(sum == C(static_cast<long long>(a[(j * r + k) * r + l])))) f.burnside++;
      }
  // The regular character (|G|, 0, …, 0) = Σ χ(1) χ.
  std::vector<Cyclotomic> reg(r, C(0));
  reg[0] = C(n);
  const auto m = ma::decompose(T, reg);
  for (std::size_t x = 0; x < T.size(); x++)
    if (!(m[x] == C(static_cast<long long>(T.degree(x))))) f.regular++;
  // Frobenius–Schur: ν ∈ {−1, 0, 1}, and Σ_χ ν_χ χ(h) = #{g : g² = h}.
  std::vector<Rational> nu;
  for (std::size_t x = 0; x < T.size(); x++) {
    nu.push_back(ma::frobenius_schur(G, T, x));
    if (!(nu.back() == Rational(1) || nu.back() == Rational(0) || nu.back() == Rational(-1))) f.indicators++;
  }
  for (std::size_t j = 0; j < r; j++) {
    long long roots = 0;
    for (std::uint32_t g = 0; g < G.n; g++)
      if (G.mul(g, g) == T.representatives[j]) roots++;
    Cyclotomic sum = 0;
    for (std::size_t x = 0; x < T.size(); x++) sum += Cyclotomic(nu[x]) * T(x, j);
    if (!(sum == C(roots))) f.indicators++;
  }
  // Real characters (ν ≠ 0) and real classes (K = K⁻¹) are equal in number.
  std::size_t real_chars = 0, real_classes = 0;
  for (std::size_t x = 0; x < T.size(); x++) {
    bool real = true;
    for (std::size_t j = 0; j < r; j++)
      if (!(T(x, j) == T(x, j).conj())) real = false;
    real_chars += real;
    if (real != !nu[x].is_zero()) f.real++;
  }
  for (std::size_t j = 0; j < r; j++) real_classes += T.inverse_class[j] == j;
  if (real_chars != real_classes) f.real++;
  if (f.orthogonality + f.degrees + f.inverses + f.linear + f.galois + f.burnside + f.regular + f.indicators + f.real != before) where += G.name + " ";
}

std::vector<std::size_t> degrees(const ma::CharacterTable& T) {
  std::vector<std::size_t> d;
  for (std::size_t x = 0; x < T.size(); x++) d.push_back(T.degree(x));
  return d;
}
bool all_rational(const ma::CharacterTable& T) {
  for (const auto& row : T.values)
    for (const auto& v : row)
      if (!v.is_rational()) return false;
  return true;
}

// Two tables equal up to an order of the classes (keeping their sizes) and of the characters.
bool same_table(const ma::CharacterTable& A, const ma::CharacterTable& B) {
  const std::size_t r = A.classes.size();
  if (B.classes.size() != r || A.order != B.order) return false;
  std::vector<std::size_t> pi(r);
  std::iota(pi.begin(), pi.end(), 0);
  do {
    if (pi[0] != 0) continue;
    bool sizes = true;
    for (std::size_t j = 0; j < r; j++) sizes &= A.class_sizes[j] == B.class_sizes[pi[j]];
    if (!sizes) continue;
    std::vector<char> used(r, 0);
    bool all = true;
    for (std::size_t x = 0; x < r && all; x++) {
      bool found = false;
      for (std::size_t y = 0; y < r && !found; y++) {
        if (used[y]) continue;
        bool eq = true;
        for (std::size_t j = 0; j < r && eq; j++) eq = A(x, j) == B(y, pi[j]);
        if (eq) used[y] = 1, found = true;
      }
      all = found;
    }
    if (all) return true;
  } while (std::next_permutation(pi.begin(), pi.end()));
  return false;
}

// The rows of a table as functions on the elements.
std::vector<std::vector<Cyclotomic>> by_element(const ma::CharacterTable& T) {
  std::vector<std::vector<Cyclotomic>> rows;
  for (std::size_t x = 0; x < T.size(); x++) {
    std::vector<Cyclotomic> row;
    for (std::uint32_t g = 0; g < T.order; g++) row.push_back(T.value(x, g));
    rows.push_back(std::move(row));
  }
  return rows;
}

// SL(2, 3), acting on the eight non-zero vectors of 𝔽_3², generated by [[1, 1], [0, 1]] and [[0, −1], [1, 0]].
ma::FiniteGroup sl_2_3() {
  auto act = [](int a, int b, int c, int d) {
    ma::Permutation p(8);
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++) {
        if (x == 0 && y == 0) continue;
        p[static_cast<std::size_t>(3 * x + y - 1)] = static_cast<std::uint32_t>(3 * ((a * x + b * y) % 3) + (c * x + d * y) % 3 - 1);
      }
    return p;
  };
  return ma::permutation_group({act(1, 1, 0, 1), act(0, 2, 1, 0)}, 8, "SL(2,3)").group;
}
// PSL(2, 7) on the projective line over 𝔽_7 (7 = ∞), generated by x ↦ x + 1 and x ↦ −1/x.
ma::FiniteGroup psl_2_7() {
  return ma::permutation_group({{1, 2, 3, 4, 5, 6, 0, 7}, {7, 6, 3, 2, 5, 4, 1, 0}}, 8, "PSL(2,7)").group;
}

// Each of `want` matches exactly one of `got`, and conversely.
bool same_rows(const std::vector<std::vector<Cyclotomic>>& got, const std::vector<std::vector<Cyclotomic>>& want) {
  if (got.size() != want.size()) return false;
  std::vector<char> used(got.size(), 0);
  for (const auto& w : want) {
    bool found = false;
    for (std::size_t y = 0; y < got.size() && !found; y++)
      if (!used[y] && got[y] == w) used[y] = 1, found = true;
    if (!found) return false;
  }
  return true;
}

}  // namespace

int main() {
  morphe::Suite s;

  // ---- the characteristic polynomial modulo p ---------------------------------------------------------------------------------
  {
    std::mt19937_64 g(5);
    int bad = 0;
    for (int t = 0; t < 60; t++) {
      const std::uint64_t p = std::vector<std::uint64_t>{61, 97, 101, 7}[t % 4];
      const std::size_t d = 1 + g() % 7;
      ma::detail::FpRows A(d, std::vector<std::uint64_t>(d));
      for (auto& row : A)
        for (auto& x : row) x = t % 3 == 0 ? (g() % 3 == 0 ? g() % p : 0) : g() % p;
      const auto poly = ma::detail::fp_charpoly(A, p);
      if (poly.size() != d + 1 || poly.back() != 1) bad++;
      for (std::uint64_t lambda : {0ull, 1ull, 5ull, 17ull}) {
        ma::Matrix B(d, d);
        for (std::size_t i = 0; i < d; i++)
          for (std::size_t j = 0; j < d; j++) B(i, j) = Integer(static_cast<long long>(i == j ? lambda : 0) - static_cast<long long>(A[i][j]));
        if (morphe::arithmetic::residue(ma::determinant(B), p) != ma::detail::fp_evaluate(poly, lambda % p, p)) bad++;
      }
    }
    s.equal("characters: the Hessenberg characteristic polynomial modulo p equals det(λI − A) over ℤ, 60 matrices (failures)", bad, 0);
  }

  // ---- the generic properties -------------------------------------------------------------------------------------------------
  const ma::FiniteGroup C2 = ma::cyclic(2), C3 = ma::cyclic(3), C4 = ma::cyclic(4);
  const ma::FiniteGroup V4 = ma::direct_product(C2, C2);
  const ma::FiniteGroup C2_3 = ma::direct_product(V4, C2);
  const ma::FiniteGroup C2_5 = ma::direct_product(ma::direct_product(C2_3, C2), C2);
  const ma::FiniteGroup F21 = ma::permutation_group({{1, 2, 3, 4, 5, 6, 0}, {0, 2, 4, 6, 1, 3, 5}}, 7, "F21").group;
  std::vector<ma::FiniteGroup> groups = {ma::FiniteGroup{}, C2, C3, C4, ma::cyclic(6), ma::cyclic(12), V4, C2_3, C2_5,
                                         ma::symmetric(3), ma::dihedral(4), ma::quaternion(), ma::alternating(4), ma::dihedral(5),
                                         ma::dihedral(6), ma::direct_product(C3, C3), F21, ma::symmetric(4), ma::alternating(5),
                                         ma::direct_product(ma::quaternion(), C3), ma::direct_product(ma::symmetric(3), C3),
                                         ma::symmetric(5), ma::alternating(6), ma::symmetric(6), sl_2_3(), psl_2_7()};
  std::map<std::string, ma::CharacterTable> tables;
  {
    Generic f;
    std::string where;
    for (const auto& G : groups) {
      const auto T = ma::character_table(G);
      generic(G, T, f, where);
      tables.emplace(G.name, T);
    }
    s.require("characters: as many characters as classes, and Σ χ(1)² = |G|, on 26 groups of orders 1 to 720", f.degrees == 0, where);
    s.equal("characters: both orthogonality relations, exactly (failures)", f.orthogonality, 0);
    s.equal("characters: χ(g⁻¹) = conj χ(g) (failures)", f.inverses, 0);
    s.equal("characters: linear characters are homomorphisms, on the groups of order ≤ 200 (failures)", f.linear, 0);
    s.equal("characters: the Galois group of ℚ(ζ_e) permutes the characters (failures)", f.galois, 0);
    s.equal("characters: the class multiplication coefficients equal (h_j h_k/|G|) Σ_χ χ(g_j)χ(g_k)conj χ(g_l)/χ(1) (failures)", f.burnside, 0);
    s.equal("characters: the regular character is Σ_χ χ(1) χ (failures)", f.regular, 0);
    s.equal("characters: Frobenius–Schur indicators lie in {−1, 0, 1} and count the square roots of every element (failures)", f.indicators, 0);
    s.equal("characters: as many real characters as real classes, the real ones those with ν ≠ 0 (failures)", f.real, 0);
    s.require("characters: the primes are ≡ 1 mod e and exceed 2√|G|", [&] {
      for (const auto& [name, T] : tables)
        if ((T.prime - 1) % T.exponent != 0 || T.prime * T.prime <= 4 * T.order) return false;
      return true;
    }());
  }

  // ---- known tables by their invariants ------------------------------------------------------------------------------------------
  {
    // C_n: every character is linear, and χ ↦ χ(generator) is a bijection onto the n-th roots of unity.
    int bad = 0;
    for (std::size_t n = 1; n <= 12; n++) {
      const auto T = ma::character_table(ma::cyclic(n));
      if (T.size() != n) bad++;
      std::vector<char> hit(n, 0);
      for (std::size_t x = 0; x < T.size(); x++) {
        if (T.degree(x) != 1) bad++;
        for (std::size_t k = 0; k < n; k++)
          if (T.value(x, n > 1 ? 1u : 0u) == Cyclotomic::zeta(static_cast<unsigned>(n), static_cast<long long>(k))) hit[k]++;
      }
      if (std::count(hit.begin(), hit.end(), 1) != static_cast<long>(n)) bad++;
    }
    s.equal("characters: C_n for n ≤ 12 takes each n-th root of unity at a generator exactly once (failures)", bad, 0);

    const auto& Tv = tables.at(V4.name);
    bool pm = Tv.size() == 4 && all_rational(Tv);
    for (const auto& row : Tv.values)
      for (const auto& v : row) pm &= v == C(1) || v == C(-1);
    s.require("characters: V_4 = C_2 × C_2 has four linear characters with values ±1", pm && degrees(Tv) == std::vector<std::size_t>{1, 1, 1, 1});

    const auto& T3 = tables.at("S3");
    s.require("characters: S_3 has degrees (1, 1, 2) and a rational table", degrees(T3) == std::vector<std::size_t>{1, 1, 2} && all_rational(T3));

    const ma::FiniteGroup D4 = ma::dihedral(4), Q8 = ma::quaternion();
    const auto &Td = tables.at(D4.name), &Tq = tables.at(Q8.name);
    s.require("characters: D_4 and Q_8 have degrees (1, 1, 1, 1, 2) and rational tables",
              degrees(Td) == std::vector<std::size_t>{1, 1, 1, 1, 2} && degrees(Tq) == degrees(Td) && all_rational(Td) && all_rational(Tq));
    s.require("characters: D_4 and Q_8 have the same table, up to the order of classes and characters", same_table(Td, Tq));
    s.require("characters: … and are told apart by the indicator of the degree-2 character, +1 for D_4 and −1 for Q_8",
              ma::frobenius_schur(D4, Td, 4) == Rational(1) && ma::frobenius_schur(Q8, Tq, 4) == Rational(-1));
    s.require("characters: D_4 and C_2^3 have different tables", !same_table(Td, tables.at(C2_3.name)));

    const ma::FiniteGroup A4 = ma::alternating(4);
    const auto& Ta = tables.at(A4.name);
    std::vector<Cyclotomic> at3;
    std::uint32_t three = 0;
    while (A4.element_order(three) != 3) three++;
    for (std::size_t x = 0; x < Ta.size(); x++)
      if (Ta.degree(x) == 1) at3.push_back(Ta.value(x, three));
    const Cyclotomic w = Cyclotomic::zeta(3);
    const bool cube = at3.size() == 3 && std::count(at3.begin(), at3.end(), C(1)) == 1 && std::count(at3.begin(), at3.end(), w) == 1 &&
                      std::count(at3.begin(), at3.end(), w * w) == 1;
    s.require("characters: A_4 has degrees (1, 1, 1, 3), its linear characters take 1, ζ_3, ζ_3² at a 3-cycle",
              degrees(Ta) == std::vector<std::size_t>{1, 1, 1, 3} && cube && !all_rational(Ta));

    const ma::FiniteGroup A5 = ma::alternating(5);
    const auto& T5 = tables.at(A5.name);
    std::uint32_t five = 0;
    while (A5.element_order(five) != 5) five++;
    int golden = 0;
    for (std::size_t x = 0; x < T5.size(); x++)
      if (T5.degree(x) == 3) {
        const Cyclotomic v = T5.value(x, five);
        golden += !v.is_rational() && v * v - v - C(1) == C(0);
      }
    s.require("characters: A_5 has degrees (1, 3, 3, 4, 5), the degree-3 characters take roots of x² − x − 1 at a 5-cycle",
              degrees(T5) == std::vector<std::size_t>{1, 3, 3, 4, 5} && golden == 2);

    const auto& Tf = tables.at(F21.name);
    std::uint32_t seven = 0;
    while (F21.element_order(seven) != 7) seven++;
    int b7 = 0;
    for (std::size_t x = 0; x < Tf.size(); x++)
      if (Tf.degree(x) == 3) {
        const Cyclotomic v = Tf.value(x, seven);
        b7 += v * v + v + C(2) == C(0);
      }
    s.require("characters: the Frobenius group of order 21 has degrees (1, 1, 1, 3, 3), the degree-3 characters taking roots of x² + x + 2",
              degrees(Tf) == std::vector<std::size_t>{1, 1, 1, 3, 3} && b7 == 2);
  }

  {
    const ma::FiniteGroup S = sl_2_3();
    const auto& T = tables.at(S.name);
    int quaternionic = 0, nonreal = 0;
    for (std::size_t x = 0; x < T.size(); x++)
      if (T.degree(x) == 2) {
        const Rational nu = ma::frobenius_schur(S, T, x);
        quaternionic += nu == Rational(-1), nonreal += nu == Rational(0);
      }
    s.require("characters: SL(2, 3) has order 24, degrees (1, 1, 1, 2, 2, 2, 3), and one quaternionic and two non-real characters of degree 2",
              S.order() == 24 && degrees(T) == std::vector<std::size_t>{1, 1, 1, 2, 2, 2, 3} && quaternionic == 1 && nonreal == 2);
    const ma::FiniteGroup P = psl_2_7();
    const auto& Tp = tables.at(P.name);
    std::uint32_t seven = 0;
    while (P.element_order(seven) != 7) seven++;
    int roots = 0;
    for (std::size_t x = 0; x < Tp.size(); x++)
      if (Tp.degree(x) == 3) {
        const Cyclotomic v = Tp.value(x, seven);
        roots += v * v + v + C(2) == C(0);
      }
    s.require("characters: PSL(2, 7) has order 168, degrees (1, 3, 3, 6, 7, 8), its degree-3 characters taking roots of x² + x + 2",
              P.order() == 168 && degrees(Tp) == std::vector<std::size_t>{1, 3, 3, 6, 7, 8} && roots == 2);
  }

  // ---- symmetric groups against Murnaghan–Nakayama ---------------------------------------------------------------------------------
  for (std::size_t n : {4u, 5u, 6u}) {
    const ma::PermutationGroup P = ma::symmetric_permutations(n);
    const auto T = ma::character_table(P.group);
    std::vector<std::vector<Cyclotomic>> want;
    for (const auto& lambda : ma::partitions(static_cast<int>(n))) {
      std::vector<Cyclotomic> row;
      for (std::size_t j = 0; j < T.classes.size(); j++)
        row.push_back(Cyclotomic(ma::murnaghan_nakayama(lambda, ma::cycle_type(P.elements[T.representatives[j]]))));
      want.push_back(std::move(row));
    }
    s.require("characters: the table of S_" + std::to_string(n) + " is the Murnaghan–Nakayama table, row for row", same_rows(T.values, want) && all_rational(T));
  }

  // ---- direct products -------------------------------------------------------------------------------------------------------------
  {
    int bad = 0;
    std::string where;
    const std::vector<std::pair<ma::FiniteGroup, ma::FiniteGroup>> pairs = {
        {C2, C2}, {ma::symmetric(3), C3}, {ma::quaternion(), C2}, {C4, C3}, {ma::dihedral(4), ma::symmetric(3)}, {ma::alternating(4), C4}};
    for (const auto& [G, H] : pairs) {
      const auto TG = ma::character_table(G), TH = ma::character_table(H), TP = ma::character_table(ma::direct_product(G, H));
      std::vector<std::vector<Cyclotomic>> want;
      for (std::size_t x = 0; x < TG.size(); x++)
        for (std::size_t y = 0; y < TH.size(); y++) {
          std::vector<Cyclotomic> row;
          for (std::uint32_t g = 0; g < G.n; g++)
            for (std::uint32_t h = 0; h < H.n; h++) row.push_back(TG.value(x, g) * TH.value(y, h));
          want.push_back(std::move(row));
        }
      if (!same_rows(by_element(TP), want)) bad++, where += G.name + "×" + H.name + " ";
    }
    s.require("characters: the table of G × H is {χ ⊗ ψ}, for six products up to order 48", bad == 0, where);
  }
  return s.done();
}
