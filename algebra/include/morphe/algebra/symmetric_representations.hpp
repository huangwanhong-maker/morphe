// The irreducible representations of the symmetric group S_n, exactly: Young's seminormal form,
// the Murnaghan–Nakayama rule, the hook length formula, and the Fourier transform on S_n.
//
// The irreducible complex representations of S_n are indexed by the partitions λ of n, and each
// is realised over ℚ. Young's seminormal form realises the representation ρ_λ on a space with
// one basis vector v_T for each standard Young tableau T of shape λ. Write c_T(k) = column −
// row of the cell holding k (its content), and s_i = (i i+1) for the Coxeter generators. Then
//
//   ρ(s_i) v_T = v_T                         if i and i+1 lie in one row of T,
//   ρ(s_i) v_T = −v_T                        if they lie in one column,
//   ρ(s_i) v_T = α v_T + v_{s_i T}           if i+1 lies in a lower row than i,
//   ρ(s_i) v_T = α v_T + (1 − α²) v_{s_i T}  if i+1 lies in a higher row than i,
//
// where α = 1 / (c_T(i+1) − c_T(i)) and s_i T exchanges i and i+1 (a standard tableau in the
// last two cases). The diagonal is forced by the Jucys–Murphy elements X_k = Σ_{j<k} (j k),
// which act on v_T by c_T(k) and satisfy X_{i+1} s_i = s_i X_i + 1; the off-diagonal pair is
// Young's orthogonal form (coefficient √(1 − α²) both ways) after rescaling v_T by
// λ_T = Π √(1 − 1/(c_T(k) − c_T(j))²) over the pairs j < k with k in a higher row than j, which
// changes by exactly the factor √(1 − α²) under s_i (Okounkov and Vershik 1996, "A new approach
// to representation theory of symmetric groups", Selecta Math. (N.S.) 2, 581–605). Every entry
// is rational. Remark: the opposite choice (1 when i+1 lies higher, 1 − α² when lower) is the
// other consistent rescaling, by the pairs with k in a lower row and |c_T(k) − c_T(j)| ≥ 2, and
// gives an equivalent representation; the convention here is fixed so that results are
// reproducible. The matrices satisfy the Coxeter relations s_i² = (s_i s_{i+1})³ = (s_i s_j)² = 1
// (|i − j| ≥ 2), so ρ(σ) is the product along any word for σ: here the word of the bubble sort,
// at most n(n−1)/2 generators, each applied in O(d) per column as it has at most two non-zero
// entries in each row. The checks hold ρ to the relations, to ρ(σ τ) = ρ(σ) ρ(τ), and its trace
// to the Murnaghan–Nakayama rule.
//
// The Murnaghan–Nakayama rule: χ^λ(μ) = Σ_ξ (−1)^{ht ξ} χ^{λ∖ξ}(μ ∖ μ_1), the sum over the rim
// hooks ξ of length μ_1 of λ, ht ξ the number of its rows less one (Sagan 2001, The Symmetric
// Group, 2nd ed., Springer, §4.10). It is computed on β-numbers: with ℓ = len λ, the beads
// β_i = λ_i + ℓ − 1 − i; removing a rim hook of length k moves a bead from b to an empty b − k ≥ 0,
// and its height is the number of beads strictly between. Memoised on (beads, parts used).
//
// The hook length formula: dim ρ_λ = n! / Π_{cells} h(i, j), h the hook length (Frame, Robinson and
// Thrall 1954, "The hook graphs of the symmetric group", Canad. J. Math. 6, 316–324; Sagan 2001,
// §3.10). Also the number of standard tableaux of shape λ, which the checks compare.
//
// The Fourier transform on S_n (Diaconis 1988, Group Representations in Probability and
// Statistics, IMS, ch. 2), for f: S_n → ℚ,
//
//   f̂(λ) = Σ_σ f(σ) ρ_λ(σ),   f(σ) = (1/n!) Σ_λ d_λ tr(ρ_λ(σ⁻¹) f̂(λ)),
//
// and with (f ∗ g)(σ) = Σ_τ f(σ τ⁻¹) g(τ), the product of the group algebra, (f ∗ g)^ = f̂ ĝ. The
// inverse transform of the matrix unit e_{ST} of block λ is E^λ_{ST} = (d_λ/n!) Σ_σ ρ_λ(σ⁻¹)_{TS} σ;
// by the convolution theorem these are matrix units of ℚ[S_n], E_{ST} E_{UV} = δ_{TU} E_{SV}, the
// explicit isomorphism ℚ[S_n] ≅ ⊕_λ M_{d_λ}(ℚ), and the E_{TT} are primitive idempotents whose
// sum over T is the central idempotent of λ. Each transform costs O(n! Σ_λ d_λ²) = O((n!)²)
// operations on rationals, the images ρ_λ(σ) being formed by one sparse product each along a
// breadth-first walk of the Cayley graph of the s_i. Permutations are numbered as symmetric(n)
// numbers them, and compose applicatively, (σ τ)(i) = σ(τ(i)).

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/rational.hpp>

#include "field_matrix.hpp"
#include "finite_group.hpp"

namespace morphe::algebra {

using arithmetic::Rational;
using Partition = std::vector<int>;  // parts in decreasing order, all positive

inline bool is_partition(const Partition& lambda) {
  for (std::size_t i = 0; i < lambda.size(); i++)
    if (lambda[i] <= 0 || (i > 0 && lambda[i] > lambda[i - 1])) return false;
  return true;
}
inline int partition_size(const Partition& lambda) { return std::accumulate(lambda.begin(), lambda.end(), 0); }

// The conjugate partition λ': λ'_j = #{i : λ_i > j}.
inline Partition conjugate(const Partition& lambda) {
  Partition c(lambda.empty() ? 0 : static_cast<std::size_t>(lambda[0]), 0);
  for (int part : lambda)
    for (int j = 0; j < part; j++) c[static_cast<std::size_t>(j)]++;
  return c;
}

// The partitions of n in reverse lexicographic order: (n) first, (1, …, 1) last.
inline std::vector<Partition> partitions(int n) {
  if (n < 0) throw std::invalid_argument("morphe::algebra::partitions: n < 0");
  std::vector<Partition> out;
  Partition cur;
  std::function<void(int, int)> go = [&](int left, int largest) {
    if (left == 0) {
      out.push_back(cur);
      return;
    }
    for (int k = std::min(left, largest); k >= 1; k--) {
      cur.push_back(k);
      go(left - k, k);
      cur.pop_back();
    }
  };
  go(n, n);
  return out;
}

// The cycle type of a permutation, as a partition.
inline Partition cycle_type(const Permutation& p) {
  std::vector<char> seen(p.size(), 0);
  Partition mu;
  for (std::size_t i = 0; i < p.size(); i++) {
    if (seen[i]) continue;
    int len = 0;
    for (std::size_t j = i; !seen[j]; j = p[j]) seen[j] = 1, len++;
    mu.push_back(len);
  }
  std::sort(mu.rbegin(), mu.rend());
  return mu;
}

// dim ρ_λ = n! / Π h(i, j), h(i, j) = λ_i − j + λ'_j − i − 1 (0-based cells).
inline arithmetic::Integer hook_length_dimension(const Partition& lambda) {
  if (!is_partition(lambda)) throw std::invalid_argument("morphe::algebra::hook_length_dimension: not a partition");
  const Partition c = conjugate(lambda);
  arithmetic::Integer hooks = 1;
  for (std::size_t i = 0; i < lambda.size(); i++)
    for (int j = 0; j < lambda[i]; j++) hooks *= arithmetic::Integer((lambda[i] - j) + (c[static_cast<std::size_t>(j)] - static_cast<int>(i)) - 1);
  return arithmetic::exact_div(arithmetic::factorial(static_cast<unsigned long long>(partition_size(lambda))), hooks);
}

// χ^λ at a permutation of cycle type μ (any order of the parts), by the Murnaghan–Nakayama rule.
inline arithmetic::Integer murnaghan_nakayama(const Partition& lambda, const Partition& mu) {
  if (!is_partition(lambda)) throw std::invalid_argument("morphe::algebra::murnaghan_nakayama: λ is not a partition");
  for (int part : mu)
    if (part <= 0) throw std::invalid_argument("morphe::algebra::murnaghan_nakayama: a part of μ is not positive");
  if (partition_size(lambda) != std::accumulate(mu.begin(), mu.end(), 0))
    throw std::invalid_argument("morphe::algebra::murnaghan_nakayama: λ and μ are partitions of different numbers");
  const int l = static_cast<int>(lambda.size());
  std::vector<int> beads;  // increasing
  for (int i = l - 1; i >= 0; i--) beads.push_back(lambda[static_cast<std::size_t>(i)] + (l - 1 - i));
  std::map<std::pair<std::vector<int>, std::size_t>, arithmetic::Integer> memo;
  std::function<arithmetic::Integer(const std::vector<int>&, std::size_t)> go = [&](const std::vector<int>& b, std::size_t k) {
    if (k == mu.size()) return arithmetic::Integer(1);
    const auto key = std::make_pair(b, k);
    const auto it = memo.find(key);
    if (it != memo.end()) return it->second;
    const int len = mu[k];
    arithmetic::Integer total = 0;
    for (std::size_t i = 0; i < b.size(); i++) {
      const int to = b[i] - len;
      if (to < 0 || std::binary_search(b.begin(), b.end(), to)) continue;
      // Beads strictly between to and b[i]: those at indices in (index of to, i).
      const std::size_t below = static_cast<std::size_t>(std::lower_bound(b.begin(), b.end(), to) - b.begin());
      const std::size_t between = i - below;
      std::vector<int> next = b;
      next.erase(next.begin() + static_cast<std::ptrdiff_t>(i));
      next.insert(std::lower_bound(next.begin(), next.end(), to), to);
      const arithmetic::Integer v = go(next, k + 1);
      total += between % 2 ? -v : v;
    }
    memo.emplace(key, total);
    return total;
  };
  return go(beads, 0);
}

// The standard Young tableaux of shape λ, each as its row word: w[k] is the row (0-based) of
// the cell holding k, for the entries 0, …, n−1; in lexicographic order of the words.
inline std::vector<std::vector<int>> standard_tableaux(const Partition& lambda) {
  if (!is_partition(lambda)) throw std::invalid_argument("morphe::algebra::standard_tableaux: not a partition");
  const int n = partition_size(lambda);
  std::vector<std::vector<int>> out;
  std::vector<int> filled(lambda.size(), 0), w(static_cast<std::size_t>(n));
  std::function<void(int)> go = [&](int k) {
    if (k == n) {
      out.push_back(w);
      return;
    }
    for (std::size_t r = 0; r < lambda.size(); r++)
      if (filled[r] < lambda[r] && (r == 0 || filled[r - 1] > filled[r])) {
        w[static_cast<std::size_t>(k)] = static_cast<int>(r);
        filled[r]++;
        go(k + 1);
        filled[r]--;
      }
  };
  go(0);
  return out;
}

class SeminormalRep {
 public:
  explicit SeminormalRep(Partition lambda) : lambda_(std::move(lambda)) {
    tableaux_ = standard_tableaux(lambda_);
    n_ = partition_size(lambda_);
    std::map<std::vector<int>, std::size_t> index;
    for (std::size_t t = 0; t < tableaux_.size(); t++) index[tableaux_[t]] = t;
    const std::size_t d = tableaux_.size();
    contents_.assign(d, std::vector<int>(static_cast<std::size_t>(n_)));
    for (std::size_t t = 0; t < d; t++) {
      std::vector<int> length(lambda_.size(), 0);
      for (int k = 0; k < n_; k++) {
        const int row = tableaux_[t][static_cast<std::size_t>(k)];
        contents_[t][static_cast<std::size_t>(k)] = length[static_cast<std::size_t>(row)]++ - row;
      }
    }
    for (int i = 0; i + 1 < n_; i++) {
      Generator g;
      g.diagonal.resize(d);
      g.partner.assign(d, d);
      g.off.resize(d);
      for (std::size_t t = 0; t < d; t++) {
        const auto& w = tableaux_[t];
        const int a = contents_[t][static_cast<std::size_t>(i)], b = contents_[t][static_cast<std::size_t>(i) + 1];
        const Rational alpha(arithmetic::Integer(1), arithmetic::Integer(b - a));
        g.diagonal[t] = alpha;
        const int ri = w[static_cast<std::size_t>(i)], rj = w[static_cast<std::size_t>(i) + 1];
        if (ri == rj || a + ri == b + rj) continue;  // one row (then α = 1) or one column (α = −1)
        std::vector<int> u = w;
        std::swap(u[static_cast<std::size_t>(i)], u[static_cast<std::size_t>(i) + 1]);
        g.partner[t] = index.at(u);
        // Row t of ρ(s_i) holds, at column s_i T, the coefficient of v_T in ρ(s_i) v_{s_i T}:
        // 1 when i+1 lies lower than i in s_i T, that is higher in T, and 1 − α² otherwise.
        g.off[t] = rj < ri ? Rational(1) : Rational(1) - alpha * alpha;
      }
      generators_.push_back(std::move(g));
    }
  }

  const Partition& shape() const { return lambda_; }
  int n() const { return n_; }
  std::size_t dimension() const { return tableaux_.size(); }
  // The standard tableaux indexing the basis, as row words.
  const std::vector<std::vector<int>>& tableaux() const { return tableaux_; }
  // c_T(k), the content of the cell holding k in the tableau of index t.
  int content(std::size_t t, int k) const { return contents_.at(t).at(static_cast<std::size_t>(k)); }

  // M ← ρ(s_i) M, s_i = (i i+1), 0 ≤ i < n − 1.
  void apply_generator(std::size_t i, FieldMatrix<Rational>& M) const {
    const Generator& g = generators_.at(i);
    const std::size_t d = dimension();
    if (M.rows != d) throw std::invalid_argument("morphe::algebra::SeminormalRep: a matrix of the wrong height");
    // Rows are changed alone (one row or one column of the tableau) or in partnered pairs.
    for (std::size_t r = 0; r < d; r++) {
      const std::size_t q = g.partner[r];
      if (q == d) {
        if (g.diagonal[r] != Rational(1))
          for (std::size_t c = 0; c < M.cols; c++) M(r, c) = -M(r, c);
        continue;
      }
      if (q < r) continue;
      for (std::size_t c = 0; c < M.cols; c++) {
        if (M(r, c).is_zero() && M(q, c).is_zero()) continue;
        Rational x = g.diagonal[r] * M(r, c) + g.off[r] * M(q, c);
        Rational y = g.diagonal[q] * M(q, c) + g.off[q] * M(r, c);
        M(r, c) = std::move(x);
        M(q, c) = std::move(y);
      }
    }
  }
  // ρ(s_i) as a matrix.
  FieldMatrix<Rational> generator(std::size_t i) const {
    FieldMatrix<Rational> M = FieldMatrix<Rational>::identity(dimension());
    apply_generator(i, M);
    return M;
  }
  // ρ(σ) for a permutation of {0, …, n−1}: σ = s_{a_k} ∘ ⋯ ∘ s_{a_1}, the a's the exchanges of the
  // bubble sort of σ (σ ∘ s_{a_1} ∘ ⋯ ∘ s_{a_k} = 1), so ρ(σ) = ρ(s_{a_k}) ⋯ ρ(s_{a_1}).
  FieldMatrix<Rational> operator()(const Permutation& sigma) const {
    if (sigma.size() != static_cast<std::size_t>(n_)) throw std::invalid_argument("morphe::algebra::SeminormalRep: a permutation of another degree");
    Permutation s = sigma;
    FieldMatrix<Rational> M = FieldMatrix<Rational>::identity(dimension());
    for (bool moved = true; moved;) {
      moved = false;
      for (std::size_t i = 0; i + 1 < s.size(); i++)
        if (s[i] > s[i + 1]) {
          std::swap(s[i], s[i + 1]);
          apply_generator(i, M);
          moved = true;
        }
    }
    return M;
  }

 private:
  struct Generator {
    std::vector<Rational> diagonal, off;
    std::vector<std::size_t> partner;  // dimension() when there is none
  };
  Partition lambda_;
  int n_ = 0;
  std::vector<std::vector<int>> tableaux_;
  std::vector<std::vector<int>> contents_;
  std::vector<Generator> generators_;
};

// S_n with its permutations, numbered as symmetric(n) numbers them (n ≤ 7): the same generators,
// the transposition (0 1) and the cycle i ↦ i + 1, closed in the same order.
inline PermutationGroup symmetric_permutations(std::size_t n) {
  if (n > 7) throw std::invalid_argument("morphe::algebra::symmetric_permutations: n > 7");
  if (n <= 1) {
    PermutationGroup P;
    P.group = symmetric(n);
    P.elements = {Permutation(n, 0)};
    return P;
  }
  Permutation t(n), c(n);
  std::iota(t.begin(), t.end(), 0u);
  std::swap(t[0], t[1]);
  for (std::size_t i = 0; i < n; i++) c[i] = static_cast<std::uint32_t>((i + 1) % n);
  return permutation_group({t, c}, n, "S" + std::to_string(n));
}

// The Fourier transform on S_n, n ≤ 7: functions are vectors indexed as the elements of group().
class SnFourier {
 public:
  explicit SnFourier(std::size_t n) : n_(n), P_(symmetric_permutations(n)) {
    for (const auto& lambda : partitions(static_cast<int>(n))) reps_.emplace_back(lambda);
    for (std::size_t k = 0; k < P_.elements.size(); k++) index_[P_.elements[k]] = k;
  }
  std::size_t n() const { return n_; }
  const PermutationGroup& group() const { return P_; }
  std::size_t order() const { return P_.elements.size(); }
  // One seminormal representation for each partition of n, in the order of partitions(n).
  const std::vector<SeminormalRep>& representations() const { return reps_; }

  // ρ_λ(σ) for every σ, along a breadth-first walk of the Cayley graph: ρ(s_i σ) = ρ(s_i) ρ(σ).
  std::vector<FieldMatrix<Rational>> images(std::size_t lambda) const {
    const SeminormalRep& R = reps_.at(lambda);
    const std::size_t N = order();
    std::vector<FieldMatrix<Rational>> out(N);
    std::vector<char> done(N, 0);
    out[0] = FieldMatrix<Rational>::identity(R.dimension());
    done[0] = 1;
    std::deque<std::size_t> queue{0};
    while (!queue.empty()) {
      const std::size_t k = queue.front();
      queue.pop_front();
      for (std::size_t i = 0; i + 1 < n_; i++) {
        Permutation s = P_.elements[k];
        for (auto& x : s) x = x == i ? static_cast<std::uint32_t>(i + 1) : (x == i + 1 ? static_cast<std::uint32_t>(i) : x);
        const std::size_t j = index_.at(s);
        if (done[j]) continue;
        out[j] = out[k];
        R.apply_generator(i, out[j]);
        done[j] = 1;
        queue.push_back(j);
      }
    }
    return out;
  }
  // f̂(λ) = Σ_σ f(σ) ρ_λ(σ), for every λ.
  std::vector<FieldMatrix<Rational>> transform(const std::vector<Rational>& f) const {
    check(f);
    std::vector<FieldMatrix<Rational>> out;
    for (std::size_t lambda = 0; lambda < reps_.size(); lambda++) {
      const auto img = images(lambda);
      const std::size_t d = reps_[lambda].dimension();
      FieldMatrix<Rational> F(d, d);
      for (std::size_t s = 0; s < order(); s++) {
        if (f[s].is_zero()) continue;
        for (std::size_t x = 0; x < d * d; x++)
          if (!img[s].a[x].is_zero()) F.a[x] += f[s] * img[s].a[x];
      }
      out.push_back(std::move(F));
    }
    return out;
  }
  // f(σ) = (1/n!) Σ_λ d_λ tr(ρ_λ(σ⁻¹) f̂(λ)).
  std::vector<Rational> inverse(const std::vector<FieldMatrix<Rational>>& fhat) const {
    if (fhat.size() != reps_.size()) throw std::invalid_argument("morphe::algebra::SnFourier::inverse: one block for each partition");
    const std::size_t N = order();
    std::vector<Rational> f(N, Rational(0));
    for (std::size_t lambda = 0; lambda < reps_.size(); lambda++) {
      const std::size_t d = reps_[lambda].dimension();
      if (fhat[lambda].rows != d || fhat[lambda].cols != d) throw std::invalid_argument("morphe::algebra::SnFourier::inverse: a block of the wrong size");
      const auto img = images(lambda);
      const Rational dl(static_cast<long long>(d));
      for (std::size_t s = 0; s < N; s++) {
        const auto& A = img[P_.group.inv(static_cast<std::uint32_t>(s))];
        Rational t = 0;
        for (std::size_t i = 0; i < d; i++)
          for (std::size_t j = 0; j < d; j++)
            if (!A(i, j).is_zero() && !fhat[lambda](j, i).is_zero()) t += A(i, j) * fhat[lambda](j, i);
        f[s] += dl * t;
      }
    }
    const Rational scale = Rational(1) / Rational(static_cast<long long>(N));
    for (auto& x : f) x *= scale;
    return f;
  }
  // E^λ_{ST} = (d_λ / n!) Σ_σ ρ_λ(σ⁻¹)_{TS} σ, the inverse transform of the matrix unit e_{ST}.
  std::vector<Rational> matrix_unit(std::size_t lambda, std::size_t S, std::size_t T) const {
    const auto img = images(lambda);
    const std::size_t d = reps_.at(lambda).dimension();
    if (S >= d || T >= d) throw std::out_of_range("morphe::algebra::SnFourier::matrix_unit");
    const Rational scale(arithmetic::Integer(static_cast<long long>(d)), arithmetic::Integer(static_cast<long long>(order())));
    std::vector<Rational> e(order());
    for (std::size_t s = 0; s < order(); s++) e[s] = scale * img[P_.group.inv(static_cast<std::uint32_t>(s))](T, S);
    return e;
  }
  // E^λ_{TT}, a primitive idempotent of ℚ[S_n].
  std::vector<Rational> minimal_idempotent(std::size_t lambda, std::size_t T) const { return matrix_unit(lambda, T, T); }

 private:
  void check(const std::vector<Rational>& f) const {
    if (f.size() != order()) throw std::invalid_argument("morphe::algebra::SnFourier: a function of the wrong length");
  }
  std::size_t n_;
  PermutationGroup P_;
  std::vector<SeminormalRep> reps_;
  std::map<Permutation, std::size_t> index_;
};

inline std::vector<FieldMatrix<Rational>> fourier_sn(const SnFourier& F, const std::vector<Rational>& f) { return F.transform(f); }
inline std::vector<Rational> inverse_fourier_sn(const SnFourier& F, const std::vector<FieldMatrix<Rational>>& fhat) { return F.inverse(fhat); }

}  // namespace morphe::algebra
