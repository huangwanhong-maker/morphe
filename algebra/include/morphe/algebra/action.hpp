// Actions of finite groups on finite sets and on the words over them.
//
// A left action of a finite group G on X = {0, …, m − 1} is a homomorphism G → Sym(X), stored as
// the table x ↦ g·x for every g. The orbits partition X, the stabiliser G_x of a point is a
// subgroup, and |G·x| |G_x| = |G| (the orbit–stabiliser theorem); the number of orbits is the
// mean number of fixed points,
//
//   |X/G| = (1/|G|) Σ_g |Fix(g)|                         (Burnside's lemma; Rotman 1995,
//                                                         An Introduction to the Theory of
//                                                         Groups, 4th ed., Theorem 3.22),
//
// which is at least |X|/|G|, since Fix(1) = X, with equality exactly when every stabiliser is
// trivial (the action is free).
//
// On words. An action on an alphabet Σ extends to Σ* letter by letter, g·(x_1…x_n) =
// (g·x_1)…(g·x_n), an action by automorphisms of the free monoid. A sign ε: G → {±1}, a
// homomorphism, lets the elements with ε(g) = −1 act by anti-automorphisms instead, reversing
// the word as well: g·(x_1…x_n) = (g·x_n)…(g·x_1). Since reversal commutes with every letter
// permutation, g ↦ (letters, ε(g)) is a homomorphism into Sym(Σ) × ℤ/2, and this is again an
// action on Σ*; it is the action of the strand symmetries of DNA (complement by letters,
// reversal by sign). A word's canonical form is the least of its images in the lexicographic
// order, its orbit's representative; on words of length n the fixed points of g number
//
//   ε(g) = +1:  |Fix_Σ(g)|^n,
//   ε(g) = −1:  |Fix_Σ(g²)|^⌊n/2⌋ · |Fix_Σ(g)|^(n mod 2),
//
// the second because w = g·w pairs position i with n − 1 − i (w_i = g·w_{n−1−i}, so
// w_{n−1−i} = g²·w_{n−1−i}) and fixes the middle letter of an odd word; with Burnside's lemma
// these count the orbits of Σ^n exactly, without enumerating them.
//
// Every operation is exact: tables of integers. Building an action from the images of all the
// elements checks the homomorphism property in O(|G|² |X|).

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/arithmetic/integer.hpp>

#include "finite_group.hpp"

namespace morphe::algebra {

struct GroupAction {
  const FiniteGroup* group = nullptr;
  std::size_t degree = 0;                 // |X|
  std::vector<std::uint32_t> table;       // table[g * degree + x] = g·x
  std::vector<std::int8_t> sign;          // ε(g) ∈ {+1, −1}; on words, −1 reverses

  std::uint32_t act(std::uint32_t g, std::uint32_t x) const { return table[static_cast<std::size_t>(g) * degree + x]; }
  bool reverses(std::uint32_t g) const { return sign[g] < 0; }
  std::size_t order() const { return group->order(); }

  // From the permutation of X of every element (images[g][x] = g·x) and, optionally, the sign of
  // every element; throws std::invalid_argument unless it is an action and ε a homomorphism.
  static GroupAction from_images(const FiniteGroup& G, std::size_t degree, const std::vector<std::vector<std::uint32_t>>& images,
                                 std::vector<std::int8_t> sign = {}) {
    GroupAction a;
    a.group = &G;
    a.degree = degree;
    if (images.size() != G.order()) throw std::invalid_argument("morphe::algebra::GroupAction: one image per element is needed");
    a.table.resize(G.order() * degree);
    for (std::uint32_t g = 0; g < G.order(); g++) {
      if (images[g].size() != degree) throw std::invalid_argument("morphe::algebra::GroupAction: an image is not a map of X");
      std::vector<char> seen(degree, 0);
      for (std::uint32_t x = 0; x < degree; x++) {
        const std::uint32_t y = images[g][x];
        if (y >= degree || seen[y]++) throw std::invalid_argument("morphe::algebra::GroupAction: an image is not a permutation");
        a.table[static_cast<std::size_t>(g) * degree + x] = y;
      }
    }
    a.sign = sign.empty() ? std::vector<std::int8_t>(G.order(), 1) : std::move(sign);
    if (a.sign.size() != G.order()) throw std::invalid_argument("morphe::algebra::GroupAction: one sign per element is needed");
    for (std::uint32_t x = 0; x < degree; x++)
      if (a.act(FiniteGroup::identity(), x) != x) throw std::invalid_argument("morphe::algebra::GroupAction: the identity moves a point");
    if (a.sign[FiniteGroup::identity()] != 1) throw std::invalid_argument("morphe::algebra::GroupAction: the identity reverses");
    for (std::uint32_t g = 0; g < G.order(); g++)
      for (std::uint32_t h = 0; h < G.order(); h++) {
        if (a.sign[G.mul(g, h)] != a.sign[g] * a.sign[h]) throw std::invalid_argument("morphe::algebra::GroupAction: the sign is not a homomorphism");
        for (std::uint32_t x = 0; x < degree; x++)
          if (a.act(G.mul(g, h), x) != a.act(g, a.act(h, x))) throw std::invalid_argument("morphe::algebra::GroupAction: (gh)·x ≠ g·(h·x)");
      }
    return a;
  }

  // The orbit of a point, in increasing order.
  std::vector<std::uint32_t> orbit(std::uint32_t x) const {
    std::vector<std::uint32_t> out;
    for (std::uint32_t g = 0; g < order(); g++) out.push_back(act(g, x));
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }
  // The orbit of each point as the index of its orbit, the orbits numbered by their least point.
  std::vector<std::uint32_t> orbit_index() const {
    std::vector<std::uint32_t> id(degree, UINT32_MAX);
    std::uint32_t next = 0;
    for (std::uint32_t x = 0; x < degree; x++) {
      if (id[x] != UINT32_MAX) continue;
      for (std::uint32_t y : orbit(x)) id[y] = next;
      next++;
    }
    return id;
  }
  std::size_t orbit_count() const {
    const auto id = orbit_index();
    return id.empty() ? 0 : *std::max_element(id.begin(), id.end()) + 1;
  }
  std::vector<std::uint32_t> stabiliser(std::uint32_t x) const {
    std::vector<std::uint32_t> out;
    for (std::uint32_t g = 0; g < order(); g++)
      if (act(g, x) == x) out.push_back(g);
    return out;
  }
  std::size_t fixed_points(std::uint32_t g) const {
    std::size_t n = 0;
    for (std::uint32_t x = 0; x < degree; x++) n += act(g, x) == x;
    return n;
  }
  // Burnside's count of the orbits of X: (1/|G|) Σ_g |Fix(g)|, exact.
  std::size_t burnside() const {
    std::size_t total = 0;
    for (std::uint32_t g = 0; g < order(); g++) total += fixed_points(g);
    if (total % order() != 0) throw std::logic_error("morphe::algebra::GroupAction: Burnside's sum is not divisible by |G|");
    return total / order();
  }

  // ---- words ------------------------------------------------------------------------------------
  template <class Word>
  Word act_on_word(std::uint32_t g, const Word& w) const {
    Word out(w.size(), typename Word::value_type{});
    const std::size_t n = w.size();
    for (std::size_t i = 0; i < n; i++) {
      const auto x = static_cast<std::uint32_t>(w[reverses(g) ? n - 1 - i : i]);
      out[i] = static_cast<typename Word::value_type>(act(g, x));
    }
    return out;
  }
  // The least image of the word, lexicographically, and an element that carries the word there.
  template <class Word>
  Word canonical(const Word& w, std::uint32_t* witness = nullptr) const {
    Word best = w;
    std::uint32_t at = FiniteGroup::identity();
    for (std::uint32_t g = 1; g < order(); g++) {
      Word image = act_on_word(g, w);
      if (image < best) { best = std::move(image); at = g; }
    }
    if (witness) *witness = at;
    return best;
  }
  // The fixed words of length n of one element, and the number of orbits of Σ^n by Burnside's
  // lemma, exactly (without bound).
  arithmetic::Integer fixed_words(std::uint32_t g, std::size_t n) const {
    if (!reverses(g)) return pow(arithmetic::Integer(static_cast<long long>(fixed_points(g))), n);
    const std::uint32_t g2 = group->mul(g, g);
    return pow(arithmetic::Integer(static_cast<long long>(fixed_points(g2))), n / 2) * pow(arithmetic::Integer(static_cast<long long>(fixed_points(g))), n % 2);
  }
  arithmetic::Integer word_orbits(std::size_t n) const {
    arithmetic::Integer total(0);
    for (std::uint32_t g = 0; g < order(); g++) total = total + fixed_words(g, n);
    return total / arithmetic::Integer(static_cast<long long>(order()));
  }

 private:
  static arithmetic::Integer pow(arithmetic::Integer base, std::size_t e) {
    arithmetic::Integer out(1);
    for (std::size_t k = 0; k < e; k++) out = out * base;
    return out;
  }
};

// The action of a group on itself by left multiplication, the regular action: free and
// transitive, so X is a torsor and the choice of a base point identifies it with G.
inline GroupAction regular_action(const FiniteGroup& G) {
  std::vector<std::vector<std::uint32_t>> images(G.order(), std::vector<std::uint32_t>(G.order()));
  for (std::uint32_t g = 0; g < G.order(); g++)
    for (std::uint32_t h = 0; h < G.order(); h++) images[g][h] = G.mul(g, h);
  return GroupAction::from_images(G, G.order(), images);
}

}  // namespace morphe::algebra
