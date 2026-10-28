// Increasing subsequences: the longest and the heaviest strictly increasing subsequence of a
// sequence, and the heaviest chain of points strictly increasing in both coordinates.
//
// - A subsequence of v_0, …, v_{n−1} is a choice of indices i_1 < … < i_m; it is (strictly)
//   increasing when v_{i_1} < … < v_{i_m}. With weights w_i its weight is Σ_j w_{i_j}; the
//   longest increasing subsequence is the heaviest one when every weight is 1.
// - The recurrence. The heaviest increasing subsequence ending at i weighs
//
//     f(i) = w_i + max(0, max{ f(j) : j < i, v_j < v_i }),
//
//   since such a subsequence is i alone or extends one ending at some j < i with v_j < v_i; the
//   heaviest of all weighs max_i f(i). Taken in index order, with a Fenwick tree (Fenwick 1994,
//   "A new data structure for cumulative frequency tables", Software: Practice and Experience
//   24, 327–336) over the ranks of the values holding prefix maxima of f, each step is one query
//   and one update of O(log n): O(n log n) time with the sort of the values, O(n) space. For unit
//   weights the bound is Fredman's (1975, "On computing the length of longest increasing
//   subsequences", Discrete Math. 11, 29–35), who showed it optimal among comparison
//   algorithms; for weights, the heaviest increasing subsequence and its O(n log n) solution are
//   Jacobson and Vo's (1992, "Heaviest increasing/common subsequence problems", CPM, LNCS 644).
// - Chains of points. A chain of points (x_i, y_i) strictly increasing in both coordinates is an
//   increasing subsequence of the y's once the points are sorted by increasing x and, within one
//   x, by decreasing y, so that two points of one x are never both taken.
// - Ties. The predecessor of i is, among the j attaining the maximum, the one of least index; the
//   answer ends at the i of greatest f, the least such index. The answer is therefore a function
//   of the input, the same on every machine.
//
// Weights are 64-bit integers of any sign, and the sums must fit in 64 bits. The values need only
// a strict weak order (operator<), equal values being incomparable for the strict increase. The
// subsequence returned is non-empty whenever the sequence is: of a sequence of negative weights
// it is the single heaviest element. Every function terminates on every input.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace morphe::combinatorics {

namespace detail {

// The best (f, index) pair: greater f, then smaller index.
struct Best {
  std::int64_t f = std::numeric_limits<std::int64_t>::min();
  std::size_t index = std::numeric_limits<std::size_t>::max();
  bool better_than(const Best& o) const { return f != o.f ? f > o.f : index < o.index; }
};

// Prefix maxima over [0, r) by a Fenwick tree of n positions.
class PrefixMax {
 public:
  explicit PrefixMax(std::size_t n) : tree_(n + 1) {}
  void update(std::size_t position, const Best& b) {
    for (std::size_t i = position + 1; i < tree_.size(); i += i & (~i + 1))
      if (b.better_than(tree_[i])) tree_[i] = b;
  }
  Best query(std::size_t r) const {  // the best over positions [0, r)
    Best b;
    for (std::size_t i = r; i > 0; i -= i & (~i + 1))
      if (tree_[i].better_than(b)) b = tree_[i];
    return b;
  }

 private:
  std::vector<Best> tree_;
};

}  // namespace detail

// The indices, increasing, of the heaviest strictly increasing subsequence of v under weights w.
template <class T>
std::vector<std::size_t> heaviest_increasing_subsequence(const std::vector<T>& v, const std::vector<std::int64_t>& w) {
  const std::size_t n = v.size();
  if (w.size() != n) throw std::invalid_argument("morphe::combinatorics::heaviest_increasing_subsequence: one weight for each value");
  if (n == 0) return {};
  // The rank of each value among the distinct values, equal values sharing one.
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
  std::vector<std::size_t> rank(n);
  std::size_t r = 0;
  for (std::size_t t = 0; t < n; t++) {
    if (t > 0 && v[order[t - 1]] < v[order[t]]) r++;
    rank[order[t]] = r;
  }
  detail::PrefixMax tree(r + 1);
  std::vector<std::int64_t> f(n);
  std::vector<std::size_t> predecessor(n, n);
  detail::Best end;
  for (std::size_t i = 0; i < n; i++) {
    const detail::Best b = tree.query(rank[i]);  // over strictly smaller values
    if (b.index < n && b.f > 0) {
      f[i] = w[i] + b.f;
      predecessor[i] = b.index;
    } else {
      f[i] = w[i];
    }
    const detail::Best here{f[i], i};
    tree.update(rank[i], here);
    if (here.better_than(end)) end = here;
  }
  std::vector<std::size_t> out;
  for (std::size_t i = end.index; i < n; i = predecessor[i]) out.push_back(i);
  std::reverse(out.begin(), out.end());
  return out;
}

// The indices of a longest strictly increasing subsequence of v.
template <class T>
std::vector<std::size_t> longest_increasing_subsequence(const std::vector<T>& v) {
  return heaviest_increasing_subsequence(v, std::vector<std::int64_t>(v.size(), 1));
}

// The indices, in chain order, of the heaviest chain of points (x_i, y_i) strictly increasing in
// both coordinates.
template <class X, class Y>
std::vector<std::size_t> heaviest_chain(const std::vector<X>& x, const std::vector<Y>& y, const std::vector<std::int64_t>& w) {
  const std::size_t n = x.size();
  if (y.size() != n || w.size() != n) throw std::invalid_argument("morphe::combinatorics::heaviest_chain: one y and one weight for each x");
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (x[a] < x[b]) return true;
    if (x[b] < x[a]) return false;
    return y[b] < y[a];
  });
  std::vector<Y> ys(n);
  std::vector<std::int64_t> ws(n);
  for (std::size_t t = 0; t < n; t++) {
    ys[t] = y[order[t]];
    ws[t] = w[order[t]];
  }
  std::vector<std::size_t> chain = heaviest_increasing_subsequence(ys, ws);
  for (auto& t : chain) t = order[t];
  return chain;
}

}  // namespace morphe::combinatorics
