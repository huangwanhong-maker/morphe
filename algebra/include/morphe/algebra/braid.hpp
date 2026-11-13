// Braid groups: words in Artin's generators, the word problem, the action on the free group, Dynnikov's
// coordinates, and the Burau representation with the Alexander polynomial of a closure.
//
// The group. B_n is presented by σ_1, …, σ_{n−1} with σ_i σ_j = σ_j σ_i for |i − j| ≥ 2 and
// σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1} (Artin 1925; Kassel and Turaev 2008, Braid Groups, Springer, §1.1). It is
// the fundamental group of the space of n unordered distinct points of the plane (Fox and Neuwirth 1962).
//
// Conventions, fixed here once.
//   - A word is read from the left in time: its first letter is the first crossing. The letter +i is σ_i, −i is
//     σ_i⁻¹, 1 ≤ i ≤ n − 1.
//   - σ_i is the counterclockwise half-exchange of the strands at positions i and i + 1 of the line of
//     projection: the left point passes below the right one (González-Meneses 2011, "Basic results on braid
//     groups", Ann. Math. Blaise Pascal 18, §3.2; Guillemot and Lairez 2026, arXiv:2601.23073). A full
//     counterclockwise turn of two strands is σ_1², so the exponent sum of a pure braid on two strands is twice
//     its winding number. Thiffeault (2010), braidlab and Bigelow take σ_i clockwise: their words are read here
//     with every letter inverted.
//   - Positions are counted from 0 in the code: the permutation of a braid sends the position a strand starts at
//     to the position it ends at.
//
// The word problem, twice.
//   - Dehornoy's handle reduction (Dehornoy 1997, "A fast method for comparing braids", Adv. Math. 125, 200–235;
//     Dehornoy, Dynnikov, Rolfsen and Wiest 2008, Ordering Braids, AMS, ch. 3). A σ_i-handle is a subword
//     σ_i^e u σ_i^{−e} in which u holds only letters of index above i. The handle that ends first holds no handle
//     inside it, so the letters σ_{i+1} in its u all have one sign d, and it is reduced by deleting its two ends
//     and replacing each σ_{i+1}^d by σ_{i+1}^{−e} σ_i^d σ_{i+1}^e, which is the same braid by the braid relation.
//     Every sequence of such reductions ends (Dehornoy's theorem), in a word with no handle; in such a word the
//     letters of the lowest index all have one sign, and a word that is not empty is then not the identity (the
//     braid is σ-positive or σ-negative in Dehornoy's order). So a braid is trivial exactly when its reduced word
//     is empty. The number of steps is not bounded by a proved polynomial; a bound on steps is taken, and
//     passing it is reported, never answered.
//   - Artin's action on the free group F_n = ⟨x_1, …, x_n⟩ (Artin 1947; Kassel and Turaev §1.5):
//     σ_i: x_i ↦ x_i x_{i+1} x_i⁻¹, x_{i+1} ↦ x_i, the others fixed. The action is faithful, so two braids are
//     equal exactly when the images of the generators are equal as reduced words. The images can grow
//     exponentially in the length of the braid, so this decides short words only, within a bound on length; it
//     is here as the independent check of the first method, and for the growth rate below.
//
// Dynnikov's coordinates (Dynnikov 2002, Russ. Math. Surv. 57, 592–594; the rules as Thiffeault 2010,
// "Braids of entangled particle trajectories", Chaos 20, 017516, arXiv:0906.3647, Eqs. 17–19, and Hall and
// Yurttaş 2009). A closed curve system on the disc with n punctures is 2n − 4 integers (a_1 … a_{n−2},
// b_1 … b_{n−2}), and a generator acts on them by rules that use only +, − , max and min, so that a braid of
// length ℓ acts in O(ℓ) operations on integers without bound. Thiffeault's rules are for his clockwise σ_i; this
// library's σ_i is his σ_i⁻¹. The rules are checked here by the braid relations on random coordinates and by the
// growth rate against the free group's.
//
// Burau (Burau 1936; Kassel and Turaev §3.1). The unreduced representation sends σ_i to the identity with the
// block [[1 − t, t], [1, 0]] at rows i, i + 1, over ℤ[t, t⁻¹]. Its matrices have rows that sum to 1 and
// determinant (−t)^e, e the exponent sum. It is faithful for n ≤ 3 and not for n ≥ 5 (Moody 1991; Long and Paton
// 1993; Bigelow 1999); so it does not decide equality, and is not used to. A first minor of I − ψ(β) is the
// Alexander polynomial of the closure of β, up to a unit (Burau 1936; Kassel and Turaev §3.4).
#pragma once

#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/laurent.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::algebra {

using arithmetic::Integer;
using arithmetic::Laurent;

struct Braid {
  std::size_t strands = 1;
  std::vector<int> word;  // +i for σ_i, −i for σ_i⁻¹, 1 ≤ i ≤ strands − 1

  Braid() = default;
  explicit Braid(std::size_t n, std::vector<int> w = {}) : strands(n), word(std::move(w)) { validate(); }
  void validate() const {
    if (strands == 0) throw std::invalid_argument("morphe::algebra::Braid: a braid has at least one strand");
    for (int l : word)
      if (l == 0 || static_cast<std::size_t>(std::abs(l)) >= strands)
        throw std::invalid_argument("morphe::algebra::Braid: the letter " + std::to_string(l) + " is not a generator of B_" + std::to_string(strands));
  }
  std::size_t length() const { return word.size(); }

  // As it is written: "s1 s2^-1 s1", or "1" for the empty word.
  std::string str() const {
    if (word.empty()) return "1";
    std::string out;
    for (int l : word) out += (out.empty() ? "" : " ") + std::string("s") + std::to_string(std::abs(l)) + (l < 0 ? "^-1" : "");
    return out;
  }
};

inline Braid operator*(const Braid& a, const Braid& b) {
  if (a.strands != b.strands) throw std::invalid_argument("morphe::algebra::Braid: a product of braids on " + std::to_string(a.strands) + " and " + std::to_string(b.strands) + " strands");
  Braid c = a;
  c.word.insert(c.word.end(), b.word.begin(), b.word.end());
  return c;
}
inline Braid inverse(const Braid& a) {
  Braid c(a.strands);
  for (auto it = a.word.rbegin(); it != a.word.rend(); ++it) c.word.push_back(-*it);
  return c;
}
inline Braid power(const Braid& a, long k) {
  const Braid base = k < 0 ? inverse(a) : a;
  Braid c(a.strands);
  for (long i = 0; i < std::labs(k); i++) c.word.insert(c.word.end(), base.word.begin(), base.word.end());
  return c;
}

// The word with adjacent inverse letters cancelled: the same braid, and the reduced word of the free group.
inline std::vector<int> free_reduced(const std::vector<int>& w) {
  std::vector<int> out;
  for (int l : w) {
    if (!out.empty() && out.back() == -l) out.pop_back();
    else out.push_back(l);
  }
  return out;
}

// ---- what a word shows at once ---------------------------------------------------------------------------

// Where each strand ends: p[k] is the final position of the strand that starts at position k.
inline std::vector<std::size_t> permutation(const Braid& b) {
  std::vector<std::size_t> at(b.strands);  // at[position] = the strand there
  std::iota(at.begin(), at.end(), std::size_t{0});
  for (int l : b.word) std::swap(at[static_cast<std::size_t>(std::abs(l)) - 1], at[static_cast<std::size_t>(std::abs(l))]);
  std::vector<std::size_t> p(b.strands);
  for (std::size_t position = 0; position < b.strands; position++) p[at[position]] = position;
  return p;
}
inline bool pure(const Braid& b) {
  const auto p = permutation(b);
  for (std::size_t k = 0; k < p.size(); k++)
    if (p[k] != k) return false;
  return true;
}
// The exponent sum: the homomorphism B_n → ℤ that sends every σ_i to 1.
inline long exponent_sum(const Braid& b) {
  long e = 0;
  for (int l : b.word) e += l > 0 ? 1 : -1;
  return e;
}
// c[i][j]: the signed number of crossings between the strands that start at positions i and j. For a pure
// braid c[i][j] / 2 is the winding number of the two strands about each other.
inline std::vector<std::vector<long>> crossings(const Braid& b) {
  std::vector<std::size_t> at(b.strands);
  std::iota(at.begin(), at.end(), std::size_t{0});
  std::vector<std::vector<long>> c(b.strands, std::vector<long>(b.strands, 0));
  for (int l : b.word) {
    const std::size_t i = static_cast<std::size_t>(std::abs(l)) - 1;
    const long e = l > 0 ? 1 : -1;
    c[at[i]][at[i + 1]] += e, c[at[i + 1]][at[i]] += e;
    std::swap(at[i], at[i + 1]);
  }
  return c;
}
// The braid with the strand that starts at position `strand` deleted: a word on one strand fewer. On pure braids
// this is the homomorphism induced by forgetting a point (Fadell and Neuwirth 1962).
inline Braid forget(const Braid& b, std::size_t strand) {
  if (strand >= b.strands || b.strands < 2) throw std::invalid_argument("morphe::algebra::forget: no such strand");
  Braid out(b.strands - 1);
  std::size_t position = strand;  // of the deleted strand
  for (int l : b.word) {
    const std::size_t i = static_cast<std::size_t>(std::abs(l)) - 1;
    if (i == position) position = i + 1;
    else if (i + 1 == position) position = i;
    else out.word.push_back((l > 0 ? 1 : -1) * static_cast<int>(i < position ? i + 1 : i));
  }
  return out;
}

// ---- the word problem by handle reduction ----------------------------------------------------------------

struct HandleReduction {
  std::vector<int> word;    // with no handle, when `finished`
  std::size_t steps = 0;
  bool finished = true;     // false: the bound on steps or on length was passed
};

inline HandleReduction handle_reduce(std::vector<int> w, std::size_t max_steps = 1000000, std::size_t max_length = 10000000) {
  HandleReduction r;
  for (;;) {
    // The handle that ends first: for each position q, the nearest earlier letter of index ≤ |w[q]| is its other
    // end exactly when it is the inverse of w[q].
    std::size_t p = 0, q = 0;
    bool found = false;
    for (q = 1; q < w.size() && !found; q++) {
      const int index = std::abs(w[q]);
      for (std::size_t k = q; k-- > 0;) {
        if (std::abs(w[k]) > index) continue;
        if (w[k] == -w[q]) p = k, found = true;
        break;
      }
      if (found) break;
    }
    if (!found) break;
    if (++r.steps > max_steps || w.size() > max_length) {
      r.finished = false;
      break;
    }
    const int i = std::abs(w[p]), e = w[p] > 0 ? 1 : -1;
    std::vector<int> next(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(p));
    for (std::size_t k = p + 1; k < q; k++) {
      if (std::abs(w[k]) == i + 1) {
        const int d = w[k] > 0 ? 1 : -1;
        next.push_back(-e * (i + 1)), next.push_back(d * i), next.push_back(e * (i + 1));
      } else {
        next.push_back(w[k]);
      }
    }
    next.insert(next.end(), w.begin() + static_cast<std::ptrdiff_t>(q) + 1, w.end());
    w.swap(next);
  }
  r.word = std::move(w);
  return r;
}

// Whether the braid is the identity; nothing when the bound on the reduction was passed.
inline std::optional<bool> trivial(const Braid& b, std::size_t max_steps = 1000000) {
  const HandleReduction r = handle_reduce(b.word, max_steps);
  if (!r.finished) return std::nullopt;
  return r.word.empty();
}
inline std::optional<bool> equal(const Braid& a, const Braid& b, std::size_t max_steps = 1000000) { return trivial(inverse(a) * b, max_steps); }

// Dehornoy's order: −1, 0 or +1 as a < b, a = b or a > b; a < b when a⁻¹ b is σ-positive (the letters of the
// lowest index in a reduced word of it are positive).
inline std::optional<int> compare(const Braid& a, const Braid& b, std::size_t max_steps = 1000000) {
  const HandleReduction r = handle_reduce((inverse(a) * b).word, max_steps);
  if (!r.finished) return std::nullopt;
  if (r.word.empty()) return 0;
  int lowest = 0;
  for (int l : r.word)
    if (lowest == 0 || std::abs(l) < std::abs(lowest)) lowest = l;
  return lowest > 0 ? -1 : 1;
}

// ---- Artin's action on the free group --------------------------------------------------------------------

// The images of x_1 … x_n under the automorphism of the braid, as reduced words (+j for x_j, −j for x_j⁻¹);
// nothing when an image grows beyond `max_length` letters. Braids are equal exactly when these are.
inline std::optional<std::vector<std::vector<int>>> artin_images(const Braid& b, std::size_t max_length = 2000000) {
  const std::size_t n = b.strands;
  std::vector<std::vector<int>> image(n);
  for (std::size_t j = 0; j < n; j++) image[j] = {static_cast<int>(j) + 1};
  const auto inverse_of = [](const std::vector<int>& w) {
    std::vector<int> out(w.rbegin(), w.rend());
    for (int& l : out) l = -l;
    return out;
  };
  const auto product = [](std::initializer_list<const std::vector<int>*> parts) {
    std::vector<int> out;
    for (const auto* part : parts)
      for (int l : *part) {
        if (!out.empty() && out.back() == -l) out.pop_back();
        else out.push_back(l);
      }
    return out;
  };
  // The composite φ_{l_1} ∘ … ∘ φ_{l_k} is built from the last letter to the first: the images under the letters
  // already taken are substituted into the next generator's.
  for (auto it = b.word.rbegin(); it != b.word.rend(); ++it) {
    const std::size_t i = static_cast<std::size_t>(std::abs(*it)) - 1;
    const std::vector<int> xi = image[i], xj = image[i + 1];
    if (*it > 0) {
      const std::vector<int> inv = inverse_of(xi);
      image[i] = product({&xi, &xj, &inv});  // x_i ↦ x_i x_{i+1} x_i⁻¹
      image[i + 1] = xi;                     // x_{i+1} ↦ x_i
    } else {
      const std::vector<int> inv = inverse_of(xj);
      image[i] = xj;                             // x_i ↦ x_{i+1}
      image[i + 1] = product({&inv, &xi, &xj});  // x_{i+1} ↦ x_{i+1}⁻¹ x_i x_{i+1}
    }
    if (image[i].size() > max_length || image[i + 1].size() > max_length) return std::nullopt;
  }
  return image;
}

// ---- Dynnikov's coordinates ------------------------------------------------------------------------------

// (a_1 … a_{n−2}, b_1 … b_{n−2}) of a curve system on the disc with n punctures.
struct Dynnikov {
  std::vector<Integer> a, b;
  friend bool operator==(const Dynnikov& x, const Dynnikov& y) { return x.a == y.a && x.b == y.b; }
  // Σ |a_i| + |b_i|: a norm whose growth under a braid has the rate of the braid's entropy.
  Integer norm() const {
    Integer s(0);
    for (const auto& v : a) s = s + arithmetic::abs(v);
    for (const auto& v : b) s = s + arithmetic::abs(v);
    return s;
  }
};

// The curve system from which entropy is measured: a_i = 0, b_i = −1 (Thiffeault 2010, §III.B: the loops that
// generate the fundamental group).
inline Dynnikov generating_loops(std::size_t strands) {
  Dynnikov u;
  if (strands >= 3) u.a.assign(strands - 2, Integer(0)), u.b.assign(strands - 2, Integer(-1));
  return u;
}

namespace braid_detail {
inline Integer pos(const Integer& x) { return x.sign() > 0 ? x : Integer(0); }
inline Integer neg(const Integer& x) { return x.sign() < 0 ? x : Integer(0); }
// Thiffeault's rules for his clockwise generator σ_i (clockwise = true) and for its inverse.
inline void act(Dynnikov& u, std::size_t i, bool clockwise, std::size_t n) {
  if (n < 3) return;
  auto& a = u.a;
  auto& b = u.b;
  const std::size_t last = n - 3;  // the index of a_{n−2}, b_{n−2}
  if (i == 1) {
    const Integer a1 = a[0], b1 = b[0];
    if (clockwise) {
      a[0] = -b1 + pos(a1 + pos(b1));
      b[0] = a1 + pos(b1);
    } else {
      a[0] = b1 - pos(pos(b1) - a1);
      b[0] = pos(b1) - a1;
    }
  } else if (i == n - 1) {
    const Integer an = a[last], bn = b[last];
    if (clockwise) {
      a[last] = -bn + neg(an + neg(bn));
      b[last] = an + neg(bn);
    } else {
      a[last] = bn - neg(neg(bn) - an);
      b[last] = neg(bn) - an;
    }
  } else {
    const std::size_t k = i - 2;  // a_{i−1}, b_{i−1} are at k; a_i, b_i at k + 1
    const Integer a0 = a[k], a1 = a[k + 1], b0 = b[k], b1 = b[k + 1];
    if (clockwise) {
      const Integer c = a0 - a1 - pos(b1) + neg(b0);
      a[k] = a0 - pos(b0) - pos(pos(b1) + c);
      b[k] = b1 + neg(c);
      a[k + 1] = a1 - neg(b1) - neg(neg(b0) - c);
      b[k + 1] = b0 - neg(c);
    } else {
      const Integer d = a0 - a1 + pos(b1) - neg(b0);
      a[k] = a0 + pos(b0) + pos(pos(b1) - d);
      b[k] = b1 - pos(d);
      a[k + 1] = a1 + neg(b1) + neg(neg(b0) + d);
      b[k + 1] = b0 + pos(d);
    }
  }
}
}  // namespace braid_detail

// The curve system after the braid has acted on it, letter by letter in time.
inline Dynnikov act(const Braid& braid, Dynnikov u) {
  if (u.a.size() + 2 != braid.strands && braid.strands >= 3) throw std::invalid_argument("morphe::algebra::act: coordinates of another number of strands");
  // This library's σ_i is counterclockwise: Thiffeault's σ_i⁻¹.
  for (int l : braid.word) braid_detail::act(u, static_cast<std::size_t>(std::abs(l)), l < 0, braid.strands);
  return u;
}

// ---- the Burau representation ----------------------------------------------------------------------------

using LaurentMatrix = std::vector<std::vector<Laurent>>;

inline LaurentMatrix laurent_identity(std::size_t n) {
  LaurentMatrix m(n, std::vector<Laurent>(n));
  for (std::size_t i = 0; i < n; i++) m[i][i] = Laurent(1);
  return m;
}
inline LaurentMatrix operator*(const LaurentMatrix& x, const LaurentMatrix& y) {
  const std::size_t n = x.size();
  LaurentMatrix z(n, std::vector<Laurent>(n));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < n; k++) {
      if (x[i][k].is_zero()) continue;
      for (std::size_t j = 0; j < n; j++)
        if (!y[k][j].is_zero()) z[i][j] = z[i][j] + x[i][k] * y[k][j];
    }
  return z;
}

// The unreduced Burau matrix of a braid, n × n over ℤ[t, t⁻¹].
inline LaurentMatrix burau(const Braid& b) {
  const std::size_t n = b.strands;
  LaurentMatrix m = laurent_identity(n);
  for (int l : b.word) {
    const std::size_t i = static_cast<std::size_t>(std::abs(l)) - 1;
    LaurentMatrix g = laurent_identity(n);
    if (l > 0) {
      g[i][i] = Laurent(1) - Laurent::t(), g[i][i + 1] = Laurent::t();
      g[i + 1][i] = Laurent(1), g[i + 1][i + 1] = Laurent(0);
    } else {
      g[i][i] = Laurent(0), g[i][i + 1] = Laurent(1);
      g[i + 1][i] = Laurent::t(-1), g[i + 1][i + 1] = Laurent(1) - Laurent::t(-1);
    }
    m = m * g;
  }
  return m;
}

// The determinant of a square matrix over ℤ[t, t⁻¹], by expansion along the first row: for the small matrices
// of braids on a few strands (the cost grows as n!).
inline Laurent determinant(const LaurentMatrix& m) {
  const std::size_t n = m.size();
  if (n == 0) return Laurent(1);
  if (n == 1) return m[0][0];
  if (n > 9) throw std::invalid_argument("morphe::algebra::determinant: an expansion by minors is for matrices of order at most 9");
  Laurent d;
  for (std::size_t j = 0; j < n; j++) {
    if (m[0][j].is_zero()) continue;
    LaurentMatrix minor(n - 1, std::vector<Laurent>(n - 1));
    for (std::size_t r = 1; r < n; r++)
      for (std::size_t c = 0, k = 0; c < n; c++)
        if (c != j) minor[r - 1][k++] = m[r][c];
    const Laurent term = m[0][j] * determinant(minor);
    d = j % 2 == 0 ? d + term : d - term;
  }
  return d;
}

// The Alexander polynomial of the closure of a braid, normalised to begin with a positive constant: the first
// minor of I − ψ(β) with its first row and column deleted. For a link of several components it is the
// one-variable polynomial.
inline Laurent alexander(const Braid& b) {
  const std::size_t n = b.strands;
  if (n < 2) return Laurent(1);
  const LaurentMatrix psi = burau(b);
  LaurentMatrix m(n - 1, std::vector<Laurent>(n - 1));
  for (std::size_t i = 1; i < n; i++)
    for (std::size_t j = 1; j < n; j++) m[i - 1][j - 1] = (i == j ? Laurent(1) : Laurent(0)) - psi[i][j];
  return determinant(m).normalised();
}

}  // namespace morphe::algebra
