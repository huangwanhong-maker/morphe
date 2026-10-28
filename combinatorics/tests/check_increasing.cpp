// Checks of morphe/combinatorics/increasing.hpp: the heaviest strictly increasing subsequence
// against every subsequence enumerated, on short random sequences with repeated values and
// weights of both signs; against the quadratic recurrence on long ones; the longest increasing
// subsequence the same with unit weights; chains of points strictly increasing in both
// coordinates, with ties in each, against the quadratic recurrence over all pairs; the theorem of
// Erdős and Szekeres (1935, "A combinatorial problem in geometry", Compositio Math. 2, 463–470),
// that every sequence of rs + 1 distinct numbers has an increasing subsequence of r + 1 terms or
// a decreasing one of s + 1, on random permutations, and its sharpness on the permutation of r
// decreasing blocks of s; and the answer the same when computed twice.

#include <morphe/check.hpp>
#include <morphe/combinatorics/increasing.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace mc = morphe::combinatorics;

namespace {

// Whether idx is a strictly increasing subsequence of v, and its weight.
template <class T>
bool increasing(const std::vector<T>& v, const std::vector<std::size_t>& idx) {
  for (std::size_t t = 1; t < idx.size(); t++)
    if (!(idx[t - 1] < idx[t]) || !(v[idx[t - 1]] < v[idx[t]])) return false;
  return true;
}
std::int64_t weight(const std::vector<std::int64_t>& w, const std::vector<std::size_t>& idx) {
  std::int64_t s = 0;
  for (auto i : idx) s += w[i];
  return s;
}

// The heaviest non-empty increasing subsequence by enumeration of all 2^n − 1 subsequences.
std::int64_t brute(const std::vector<int>& v, const std::vector<std::int64_t>& w) {
  const std::size_t n = v.size();
  std::int64_t best = std::numeric_limits<std::int64_t>::min();
  for (std::uint32_t mask = 1; mask < (1u << n); mask++) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < n; i++)
      if (mask >> i & 1) idx.push_back(i);
    if (increasing(v, idx)) best = std::max(best, weight(w, idx));
  }
  return best;
}

// The quadratic recurrence f(i) = w_i + max(0, max_{j < i, v_j < v_i} f(j)).
template <class T>
std::int64_t quadratic(const std::vector<T>& v, const std::vector<std::int64_t>& w) {
  std::vector<std::int64_t> f(v.size());
  std::int64_t best = std::numeric_limits<std::int64_t>::min();
  for (std::size_t i = 0; i < v.size(); i++) {
    std::int64_t m = 0;
    for (std::size_t j = 0; j < i; j++)
      if (v[j] < v[i]) m = std::max(m, f[j]);
    f[i] = w[i] + m;
    best = std::max(best, f[i]);
  }
  return best;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(1935);

  // Against enumeration: 3000 sequences of 1 to 13 values in [0, 6), weights in [−3, 9].
  {
    std::size_t bad = 0, invalid = 0;
    for (int trial = 0; trial < 3000; trial++) {
      const std::size_t n = 1 + g() % 13;
      std::vector<int> v(n);
      std::vector<std::int64_t> w(n);
      for (auto& x : v) x = static_cast<int>(g() % 6);
      for (auto& x : w) x = static_cast<std::int64_t>(g() % 13) - 3;
      const auto idx = mc::heaviest_increasing_subsequence(v, w);
      if (idx.empty() || !increasing(v, idx)) invalid++;
      else if (weight(w, idx) != brute(v, w)) bad++;
    }
    s.equal("on 3000 short sequences with repeated values and weights of both signs, the subsequences returned are non-empty and strictly increasing (failures)", invalid, std::size_t{0});
    s.equal("and their weight is the greatest of every subsequence enumerated (failures)", bad, std::size_t{0});
  }

  // Against the quadratic recurrence on long sequences; unit weights for the longest.
  {
    std::size_t bad = 0, bad_lis = 0;
    for (int trial = 0; trial < 40; trial++) {
      const std::size_t n = 500 + g() % 2500;
      std::vector<std::int64_t> v(n), w(n);
      const std::uint64_t range = trial % 2 ? 50 : 1'000'000;
      for (auto& x : v) x = static_cast<std::int64_t>(g() % range);
      for (auto& x : w) x = static_cast<std::int64_t>(g() % 1000) - 100;
      const auto idx = mc::heaviest_increasing_subsequence(v, w);
      if (!increasing(v, idx) || weight(w, idx) != quadratic(v, w)) bad++;
      const auto lis = mc::longest_increasing_subsequence(v);
      if (!increasing(v, lis) || static_cast<std::int64_t>(lis.size()) != quadratic(v, std::vector<std::int64_t>(n, 1))) bad_lis++;
    }
    s.equal("on 40 sequences of 500 to 3000 values the heaviest increasing subsequence weighs what the quadratic recurrence gives (failures)", bad, std::size_t{0});
    s.equal("and the longest has the length the recurrence gives with unit weights (failures)", bad_lis, std::size_t{0});
  }

  // Chains of points, strictly increasing in both coordinates, with ties in each.
  {
    std::size_t bad = 0;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 1 + g() % 300;
      std::vector<int> x(n), y(n);
      std::vector<std::int64_t> w(n);
      for (std::size_t i = 0; i < n; i++) {
        x[i] = static_cast<int>(g() % 40);
        y[i] = static_cast<int>(g() % 40);
        w[i] = static_cast<std::int64_t>(g() % 20) + 1;
      }
      const auto chain = mc::heaviest_chain(x, y, w);
      bool ok = !chain.empty();
      for (std::size_t t = 1; t < chain.size(); t++) ok = ok && x[chain[t - 1]] < x[chain[t]] && y[chain[t - 1]] < y[chain[t]];
      // The quadratic recurrence over all pairs, in order of x.
      std::vector<std::size_t> order(n);
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
      std::vector<std::int64_t> f(n);
      std::int64_t best = 0;
      for (std::size_t a = 0; a < n; a++) {
        std::int64_t m = 0;
        for (std::size_t b = 0; b < a; b++)
          if (x[order[b]] < x[order[a]] && y[order[b]] < y[order[a]]) m = std::max(m, f[b]);
        f[a] = w[order[a]] + m;
        best = std::max(best, f[a]);
      }
      if (!ok || weight(w, chain) != best) bad++;
    }
    s.equal("on 200 sets of points with ties in x and in y, the heaviest chain strictly increasing in both is a chain of the greatest weight (failures)", bad, std::size_t{0});
  }

  // Erdős–Szekeres on random permutations, and its sharpness.
  {
    std::size_t bad = 0, sharp = 0;
    for (std::size_t r = 1; r <= 6; r++)
      for (std::size_t q = 1; q <= 6; q++) {
        for (int trial = 0; trial < 20; trial++) {
          std::vector<std::int64_t> p(r * q + 1);
          std::iota(p.begin(), p.end(), std::int64_t{0});
          std::shuffle(p.begin(), p.end(), g);
          std::vector<std::int64_t> neg(p.size());
          for (std::size_t i = 0; i < p.size(); i++) neg[i] = -p[i];
          const std::size_t up = mc::longest_increasing_subsequence(p).size(), down = mc::longest_increasing_subsequence(neg).size();
          if (up < r + 1 && down < q + 1) bad++;
        }
        // r decreasing blocks of q: block b holds b q + q − 1, …, b q.
        std::vector<std::int64_t> e;
        for (std::size_t b = 0; b < r; b++)
          for (std::size_t j = 0; j < q; j++) e.push_back(static_cast<std::int64_t>(b * q + q - 1 - j));
        std::vector<std::int64_t> neg(e.size());
        for (std::size_t i = 0; i < e.size(); i++) neg[i] = -e[i];
        if (mc::longest_increasing_subsequence(e).size() != r || mc::longest_increasing_subsequence(neg).size() != q) sharp++;
      }
    s.equal("Erdős–Szekeres: 720 random permutations of rs + 1 terms (1 ≤ r, s ≤ 6) each have an increasing subsequence of r + 1 or a decreasing one of s + 1 (failures)", bad, std::size_t{0});
    s.equal("and r decreasing blocks of s terms, rs in all, have neither: their longest increasing subsequence is r and their longest decreasing s (failures)", sharp, std::size_t{0});
  }

  // A function of the input.
  {
    std::vector<std::int64_t> v(5000), w(5000);
    for (auto& x : v) x = static_cast<std::int64_t>(g() % 100);
    for (auto& x : w) x = static_cast<std::int64_t>(g() % 7);
    s.require("computed twice, the heaviest increasing subsequence is the same, ties broken the same way",
              mc::heaviest_increasing_subsequence(v, w) == mc::heaviest_increasing_subsequence(v, w));
    s.require("the empty sequence has the empty subsequence", mc::heaviest_increasing_subsequence(std::vector<int>{}, {}).empty());
  }
  return s.done();
}
