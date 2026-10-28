// Checks of morphe/algebra/symmetric_representations.hpp. Partitions counted against Euler's
// pentagonal recurrence; standard tableaux counted by the hook length formula, Σ_λ d_λ² = n! and
// Σ_λ d_λ the number of involutions; Young's seminormal matrices against the Coxeter relations
// for every partition of n ≤ 8, against ρ(σ τ) = ρ(σ) ρ(τ) on all of S_4 and on samples of S_6
// and S_8, and their traces against the Murnaghan–Nakayama rule (n ≤ 7); the rule against both
// orthogonality relations for n ≤ 9. The Fourier transform on S_4, S_5 and S_6: inversion,
// convolution to products, Plancherel; matrix units E_{ST} E_{UV} = δ_{TU} E_{SV}, the primitive
// idempotents E_{TT} of minimal left ideals summing to the central idempotent of λ computed
// from the Dixon–Schneider table.

#include <morphe/algebra/characters.hpp>
#include <morphe/algebra/group_algebra.hpp>
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
using Mat = ma::FieldMatrix<Rational>;

namespace {

ma::Permutation compose(const ma::Permutation& p, const ma::Permutation& q) {
  ma::Permutation r(q.size());
  for (std::size_t i = 0; i < q.size(); i++) r[i] = p[q[i]];
  return r;
}
ma::Permutation random_permutation(std::mt19937_64& g, std::size_t n) {
  ma::Permutation p(n);
  std::iota(p.begin(), p.end(), 0u);
  std::shuffle(p.begin(), p.end(), g);
  return p;
}
// A permutation of cycle type μ: consecutive blocks, each cycled.
ma::Permutation of_type(const ma::Partition& mu) {
  ma::Permutation p;
  std::uint32_t start = 0;
  for (int len : mu) {
    for (int k = 0; k < len; k++) p.push_back(start + static_cast<std::uint32_t>((k + 1) % len));
    start += static_cast<std::uint32_t>(len);
  }
  return p;
}
// z_μ = Π_i i^{m_i} m_i!, the order of the centraliser of a permutation of type μ.
Integer z_mu(const ma::Partition& mu) {
  std::map<int, unsigned long long> m;
  for (int part : mu) m[part]++;
  Integer z = 1;
  for (auto [part, k] : m) z *= morphe::arithmetic::pow(Integer(part), k) * morphe::arithmetic::factorial(k);
  return z;
}
std::vector<Rational> random_function(std::mt19937_64& g, std::size_t N) {
  std::vector<Rational> f(N);
  for (auto& x : f) x = Rational(Integer(static_cast<long long>(g() % 9) - 4), Integer(static_cast<long long>(1 + g() % 3)));
  return f;
}
bool is_zero(const std::vector<Rational>& x) {
  for (const auto& v : x)
    if (!v.is_zero()) return false;
  return true;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(29);

  // ---- partitions, tableaux and hooks -------------------------------------------------------------------------------------------
  {
    // p(n) by Euler's pentagonal recurrence p(n) = Σ_{k ≥ 1} (−1)^{k+1} (p(n − k(3k−1)/2) + p(n − k(3k+1)/2)).
    std::vector<long long> p(31, 0);
    p[0] = 1;
    for (int n = 1; n <= 30; n++)
      for (int k = 1;; k++) {
        const int a = n - k * (3 * k - 1) / 2, b = n - k * (3 * k + 1) / 2;
        if (a < 0) break;
        const long long sign = k % 2 ? 1 : -1;
        p[static_cast<std::size_t>(n)] += sign * (p[static_cast<std::size_t>(a)] + (b >= 0 ? p[static_cast<std::size_t>(b)] : 0));
      }
    int bad = 0;
    for (int n = 0; n <= 30; n++)
      if (static_cast<long long>(ma::partitions(n).size()) != p[static_cast<std::size_t>(n)]) bad++;
    s.equal("symmetric: the partitions of n ≤ 30 number p(n), by Euler's pentagonal recurrence (failures)", bad, 0);
    int hooks = 0, squares = 0, involutions = 0;
    Integer t0 = 1, t1 = 1;  // involutions: t(n) = t(n−1) + (n−1) t(n−2)
    for (int n = 1; n <= 10; n++) {
      if (n >= 2) {
        const Integer t2 = t1 + Integer(n - 1) * t0;
        t0 = t1, t1 = t2;
      }
      Integer sum_sq = 0, sum = 0;
      for (const auto& lambda : ma::partitions(n)) {
        const Integer d = ma::hook_length_dimension(lambda);
        if (Integer(static_cast<long long>(ma::standard_tableaux(lambda).size())) != d) hooks++;
        if (ma::murnaghan_nakayama(lambda, ma::Partition(static_cast<std::size_t>(n), 1)) != d) hooks++;
        sum_sq += d * d, sum += d;
      }
      if (sum_sq != morphe::arithmetic::factorial(static_cast<unsigned long long>(n))) squares++;
      if (sum != t1) involutions++;
    }
    s.equal("symmetric: standard tableaux, the hook length formula and χ^λ(1) agree for every λ ⊢ n ≤ 10 (failures)", hooks, 0);
    s.equal("symmetric: Σ_λ d_λ² = n! for n ≤ 10 (failures)", squares, 0);
    s.equal("symmetric: Σ_λ d_λ is the number of involutions of S_n, n ≤ 10 (failures)", involutions, 0);
  }

  // ---- the Murnaghan–Nakayama rule ------------------------------------------------------------------------------------------------
  {
    int rows = 0, cols = 0;
    for (int n = 1; n <= 9; n++) {
      const auto P = ma::partitions(n);
      std::vector<std::vector<Integer>> X(P.size(), std::vector<Integer>(P.size()));
      for (std::size_t a = 0; a < P.size(); a++)
        for (std::size_t b = 0; b < P.size(); b++) X[a][b] = ma::murnaghan_nakayama(P[a], P[b]);
      const Integer N = morphe::arithmetic::factorial(static_cast<unsigned long long>(n));
      for (std::size_t m = 0; m < P.size(); m++)
        for (std::size_t k = 0; k < P.size(); k++) {
          Integer c = 0;
          for (std::size_t a = 0; a < P.size(); a++) c += X[a][m] * X[a][k];
          if (c != (m == k ? z_mu(P[m]) : Integer(0))) cols++;
        }
      for (std::size_t a = 0; a < P.size(); a++)
        for (std::size_t b = 0; b < P.size(); b++) {
          Rational r = 0;
          for (std::size_t m = 0; m < P.size(); m++) r += Rational(X[a][m] * X[b][m], z_mu(P[m]));
          if (r != Rational(a == b ? 1 : 0)) rows++;
        }
    }
    s.equal("symmetric: Σ_λ χ^λ(μ) χ^λ(ν) = δ_{μν} z_μ for n ≤ 9 (failures)", cols, 0);
    s.equal("symmetric: Σ_μ χ^λ(μ) χ^κ(μ) / z_μ = δ_{λκ} for n ≤ 9 (failures)", rows, 0);
    s.require("symmetric: χ^λ is independent of the order of the parts of μ (λ = (4, 2, 1), μ = (1, 3, 2, 1) and (3, 2, 1, 1))",
              ma::murnaghan_nakayama({4, 2, 1}, {1, 3, 2, 1}) == ma::murnaghan_nakayama({4, 2, 1}, {3, 2, 1, 1}));
  }

  // ---- Young's seminormal form -----------------------------------------------------------------------------------------------------
  {
    int coxeter = 0;
    for (int n = 2; n <= 8; n++)
      for (const auto& lambda : ma::partitions(n)) {
        const ma::SeminormalRep R(lambda);
        const Mat I = Mat::identity(R.dimension());
        std::vector<Mat> S;
        for (int i = 0; i + 1 < n; i++) S.push_back(R.generator(static_cast<std::size_t>(i)));
        for (std::size_t i = 0; i < S.size(); i++) {
          if (S[i] * S[i] != I) coxeter++;
          if (i + 1 < S.size()) {
            const Mat B = S[i] * S[i + 1];
            if (B * B * B != I) coxeter++;
          }
          for (std::size_t j = i + 2; j < S.size(); j++) {
            const Mat B = S[i] * S[j];
            if (B * B != I) coxeter++;
          }
        }
      }
    s.equal("symmetric: ρ(s_i)² = (ρ(s_i)ρ(s_{i+1}))³ = (ρ(s_i)ρ(s_j))² = 1 for every λ ⊢ n ≤ 8 (failures)", coxeter, 0);

    int hom4 = 0;
    const ma::PermutationGroup P4 = ma::symmetric_permutations(4);
    for (const auto& lambda : ma::partitions(4)) {
      const ma::SeminormalRep R(lambda);
      std::vector<Mat> rho;
      for (const auto& p : P4.elements) rho.push_back(R(p));
      for (std::uint32_t a = 0; a < 24; a++)
        for (std::uint32_t b = 0; b < 24; b++)
          if (rho[P4.group.mul(a, b)] != rho[a] * rho[b] || R(compose(P4.elements[a], P4.elements[b])) != rho[a] * rho[b]) hom4++;
    }
    s.equal("symmetric: ρ_λ(σ τ) = ρ_λ(σ) ρ_λ(τ) on all of S_4, every λ (failures)", hom4, 0);
    int hom6 = 0, hom8 = 0;
    for (const auto& lambda : ma::partitions(6)) {
      const ma::SeminormalRep R(lambda);
      for (int t = 0; t < 30; t++) {
        const auto a = random_permutation(g, 6), b = random_permutation(g, 6);
        if (R(compose(a, b)) != R(a) * R(b)) hom6++;
      }
    }
    for (const ma::Partition& lambda : {ma::Partition{4, 2, 1, 1}, ma::Partition{3, 3, 2}, ma::Partition{5, 3}}) {
      const ma::SeminormalRep R(lambda);
      for (int t = 0; t < 4; t++) {
        const auto a = random_permutation(g, 8), b = random_permutation(g, 8);
        if (R(compose(a, b)) != R(a) * R(b)) hom8++;
      }
    }
    s.equal("symmetric: ρ_λ(σ τ) = ρ_λ(σ) ρ_λ(τ) on 30 random pairs for each λ ⊢ 6 (failures)", hom6, 0);
    s.equal("symmetric: … and on 4 pairs for (4, 2, 1, 1), (3, 3, 2), (5, 3) of 8, of dimensions 90, 42, 28 (failures)", hom8, 0);

    int traces = 0;
    for (int n = 1; n <= 7; n++)
      for (const auto& lambda : ma::partitions(n)) {
        const ma::SeminormalRep R(lambda);
        for (const auto& mu : ma::partitions(n))
          if (R(of_type(mu)).trace() != Rational(ma::murnaghan_nakayama(lambda, mu))) traces++;
      }
    const ma::PermutationGroup P5 = ma::symmetric_permutations(5);
    for (const auto& lambda : ma::partitions(5)) {
      const ma::SeminormalRep R(lambda);
      for (const auto& p : P5.elements)
        if (R(p).trace() != Rational(ma::murnaghan_nakayama(lambda, ma::cycle_type(p)))) traces++;
    }
    s.equal("symmetric: tr ρ_λ = χ^λ by Murnaghan–Nakayama at every cycle type, n ≤ 7, and at every element of S_5 (failures)", traces, 0);
    s.require("symmetric: S_n with its permutations is numbered as symmetric(n), n ≤ 6", [] {
      for (std::size_t n = 1; n <= 6; n++)
        if (ma::symmetric_permutations(n).group.table != ma::symmetric(n).table) return false;
      return true;
    }());
  }

  // ---- the Fourier transform on S_n ------------------------------------------------------------------------------------------------
  for (std::size_t n : {4u, 5u}) {
    const ma::SnFourier F(n);
    const ma::GroupAlgebra<Rational> A(F.group().group);
    const std::size_t N = F.order();
    int inversion = 0, convolution = 0, plancherel = 0, deltas = 0;
    for (int t = 0; t < 3; t++) {
      const auto f = random_function(g, N), h = random_function(g, N);
      const auto fh = ma::fourier_sn(F, f), hh = ma::fourier_sn(F, h);
      if (ma::inverse_fourier_sn(F, fh) != f) inversion++;
      const auto ch = ma::fourier_sn(F, A.multiply(f, h));
      for (std::size_t l = 0; l < ch.size(); l++)
        if (ch[l] != fh[l] * hh[l]) convolution++;
      Rational left = 0, right = 0;
      for (std::size_t x = 0; x < N; x++) left += f[x] * h[F.group().group.inv(static_cast<std::uint32_t>(x))];
      for (std::size_t l = 0; l < ch.size(); l++)
        right += Rational(static_cast<long long>(F.representations()[l].dimension())) * (fh[l] * hh[l]).trace();
      if (left * Rational(static_cast<long long>(N)) != right) plancherel++;
    }
    for (std::uint32_t x = 0; x < N; x += 7) {
      std::vector<Rational> delta(N, Rational(0));
      delta[x] = 1;
      const auto dh = ma::fourier_sn(F, delta);
      for (std::size_t l = 0; l < dh.size(); l++)
        if (dh[l] != F.representations()[l](F.group().elements[x])) deltas++;
    }
    const std::string Sn = "S_" + std::to_string(n);
    s.equal("symmetric: f = (1/n!) Σ_λ d_λ tr(ρ_λ(σ⁻¹) f̂(λ)) recovers f on " + Sn + " (failures)", inversion, 0);
    s.equal("symmetric: (f ∗ g)^ = f̂ ĝ on " + Sn + ", the product of the group algebra (failures)", convolution, 0);
    s.equal("symmetric: Σ_σ f(σ) g(σ⁻¹) = (1/n!) Σ_λ d_λ tr(f̂ ĝ) on " + Sn + " (failures)", plancherel, 0);
    s.equal("symmetric: the transform of a point mass δ_σ is ρ_λ(σ) on " + Sn + " (failures)", deltas, 0);
  }
  {
    const ma::SnFourier F(6);
    const auto f = random_function(g, F.order());
    s.require("symmetric: the Fourier transform on S_6 (11 blocks, 720 elements) is inverted exactly", ma::inverse_fourier_sn(F, ma::fourier_sn(F, f)) == f);
  }

  // ---- matrix units and minimal idempotents ------------------------------------------------------------------------------------------
  for (std::size_t n : {3u, 4u, 5u}) {
    const ma::SnFourier F(n);
    const ma::GroupAlgebra<Rational> A(F.group().group);
    const auto T = ma::character_table(F.group().group);
    int units = 0, central = 0, minimal = 0, unit_fourier = 0;
    auto one = A.zero();
    for (std::size_t l = 0; l < F.representations().size(); l++) {
      const auto& R = F.representations()[l];
      const std::size_t d = R.dimension();
      // E_{ST} E_{UV} = δ_{TU} E_{SV}, on S_3 and S_4 for all units, on S_5 for the diagonal ones.
      std::vector<std::vector<std::vector<Rational>>> E(d, std::vector<std::vector<Rational>>(d));
      for (std::size_t a = 0; a < d; a++)
        for (std::size_t b = 0; b < d; b++)
          if (n < 5 || a == b) E[a][b] = F.matrix_unit(l, a, b);
      for (std::size_t a = 0; a < d; a++)
        for (std::size_t b = 0; b < d; b++)
          for (std::size_t c = 0; c < d; c++)
            for (std::size_t e = 0; e < d; e++) {
              if (n == 5 && !(a == b && c == e)) continue;
              const auto P = A.multiply(E[a][b], E[c][e]);
              if (b == c ? P != E[a][e] : !is_zero(P)) units++;
            }
      // Σ_T E_{TT} is the central idempotent of the character of λ in the Dixon–Schneider table.
      auto sum = A.zero();
      for (std::size_t a = 0; a < d; a++) sum = A.add(sum, E[a][a]);
      one = A.add(one, sum);
      std::size_t chi = T.size();
      for (std::size_t x = 0; x < T.size(); x++) {
        bool match = true;
        for (std::size_t j = 0; j < T.classes.size(); j++)
          match &= T(x, j) == Cyclotomic(ma::murnaghan_nakayama(R.shape(), ma::cycle_type(F.group().elements[T.representatives[j]])));
        if (match) chi = x;
      }
      if (chi == T.size() || ma::to_cyclotomic(sum) != ma::central_idempotent(T, chi)) central++;
      // ℚ[S_n] E_{TT} is a minimal left ideal, of dimension d_λ.
      for (std::size_t a = 0; a < d; a++)
        if (A.ideal_dimension(E[a][a]) != d) minimal++;
      if (n < 5) {
        const auto Fh = ma::fourier_sn(F, E[0][d - 1]);
        for (std::size_t k = 0; k < Fh.size(); k++) {
          Mat want(F.representations()[k].dimension(), F.representations()[k].dimension());
          if (k == l) want(0, d - 1) = 1;
          if (Fh[k] != want) unit_fourier++;
        }
      }
    }
    const std::string Sn = "S_" + std::to_string(n);
    s.equal("symmetric: the matrix units of " + Sn + " multiply as E_{ST} E_{UV} = δ_{TU} E_{SV} (failures)", units, 0);
    s.equal("symmetric: Σ_T E^λ_{TT} = e_λ, the central idempotent from the character table, on " + Sn + " (failures)", central, 0);
    s.equal("symmetric: each E_{TT} generates a minimal left ideal of dimension d_λ in ℚ[" + Sn + "] (failures)", minimal, 0);
    s.require("symmetric: Σ_λ Σ_T E^λ_{TT} = 1 in ℚ[" + Sn + "]", one == A.one());
    if (n < 5) s.equal("symmetric: the transform of E^λ_{ST} is the matrix unit e_{ST} of block λ on " + Sn + " (failures)", unit_fourier, 0);
  }
  return s.done();
}
