// Finitely presented groups: words, presentations, Tietze simplification, abelianisation.
//
// A word in the free group F(x_1, ..., x_n) is a sequence of letters, the letter +i standing
// for x_i and -i for its inverse (1-based, so that 0 is never a letter). A presentation
// ⟨x_1, ..., x_n | r_1, ..., r_m⟩ is the quotient of F by the normal closure of the
// relators. Words are multiplied left to right, so the word (1, 2) is x_1 x_2.
//
// simplify() applies Tietze transformations (Lyndon and Schupp 1977, Combinatorial Group
// Theory, §II.2; Magnus, Karrass and Solitar 1966, §1.5) that preserve the group:
// - free and cyclic reduction of relators; removal of trivial and repeated relators (up to
//   cyclic permutation and inversion);
// - elimination of a generator x that occurs exactly once in some relator r = u x^ε v:
//   then x = (v u)^{-ε} in the group, which is substituted everywhere and r is dropped.
// The shortest relator with a generator occurring once is used first, and an elimination
// that would make the relators more than four times as long in all is not made.
// It always terminates (each elimination removes a generator, the other steps shorten
// relators) and returns a presentation of the same group, but it is a heuristic: no normal
// form is claimed, since the isomorphism problem for finitely presented groups is
// undecidable (Adian 1957; Rabin 1958). The kept generators are a subset of the original
// ones, so the new generators are old generators, and kept[] says which.
//
// abelianisation() gives G / [G, G] ≅ ℤ^rank ⊕ ⊕ ℤ/t_i from the matrix of exponent sums of
// the relators, reduced sparsely and then put in Smith normal form.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "chain_complex.hpp"
#include "smith.hpp"

namespace morphe::algebra {

using Word = std::vector<int>;

inline Word inverse(const Word& w) {
  Word r(w.rbegin(), w.rend());
  for (auto& x : r) x = -x;
  return r;
}

// The freely reduced form: no letter next to its inverse.
inline Word reduce(const Word& w) {
  Word r;
  for (int x : w) {
    if (!r.empty() && r.back() == -x) r.pop_back();
    else r.push_back(x);
  }
  return r;
}

// The cyclically reduced form of a freely reduced word: first and last letters not inverse.
inline Word cyclic_reduce(Word w) {
  w = reduce(w);
  std::size_t i = 0, j = w.size();
  while (j - i >= 2 && w[i] == -w[j - 1]) i++, j--;
  return Word(w.begin() + static_cast<long>(i), w.begin() + static_cast<long>(j));
}

// A canonical representative of a cyclic word up to rotation and inversion: the
// lexicographically least rotation of w or of its inverse.
inline Word canonical_cyclic(const Word& w) {
  Word best;
  for (const Word& v : {w, inverse(w)})
    for (std::size_t s = 0; s < v.size(); s++) {
      Word rot(v.begin() + static_cast<long>(s), v.end());
      rot.insert(rot.end(), v.begin(), v.begin() + static_cast<long>(s));
      if (best.empty() || rot < best) best = rot;
    }
  return best;
}

struct Presentation {
  std::size_t generators = 0;
  std::vector<Word> relators;
  std::vector<std::string> names;  // one per generator; may be empty

  std::string name(std::size_t i) const { return i < names.size() ? names[i] : "x" + std::to_string(i + 1); }
  std::string word(const Word& w) const {
    if (w.empty()) return "1";
    std::string s;
    for (std::size_t k = 0; k < w.size(); k++) s += (k ? " " : "") + name(static_cast<std::size_t>(std::abs(w[k])) - 1) + (w[k] < 0 ? "⁻¹" : "");
    return s;
  }
  std::string str() const {
    std::string s = "⟨";
    for (std::size_t i = 0; i < generators; i++) s += (i ? ", " : "") + name(i);
    s += " | ";
    for (std::size_t k = 0; k < relators.size(); k++) s += (k ? ", " : "") + word(relators[k]);
    return s + "⟩";
  }
};

struct Simplified {
  Presentation presentation;
  std::vector<std::size_t> kept;  // kept[i]: the original index of new generator i
};

inline Simplified simplify(const Presentation& P) {
  // Work with original generator numbers; eliminated ones are substituted away. Relators
  // are indexed by the generators they contain, and those in which some generator occurs
  // once are kept ordered by length, so that an elimination touches only the relators that
  // contain the eliminated generator: O(occurrences × length) each.
  std::vector<Word> rel;
  std::vector<char> alive;
  std::vector<std::set<std::size_t>> holders(P.generators + 1);  // holders[g]: relators containing g
  std::set<std::pair<std::size_t, std::size_t>> candidates;       // (length, relator) with a generator occurring once
  std::map<Word, std::size_t> canonical;                          // canonical form -> a live relator
  std::vector<Word> form;                                          // each relator's canonical form
  std::size_t total = 0;
  auto single = [](const Word& w) {
    std::map<int, int> count;
    for (int x : w) count[std::abs(x)]++;
    for (std::size_t p = 0; p < w.size(); p++)
      if (count[std::abs(w[p])] == 1) return p;
    return SIZE_MAX;
  };
  auto remove = [&](std::size_t k) {
    alive[k] = 0;
    total -= rel[k].size();
    candidates.erase({rel[k].size(), k});
    for (int x : rel[k]) holders[static_cast<std::size_t>(std::abs(x))].erase(k);
    const auto it = canonical.find(form[k]);
    if (it != canonical.end() && it->second == k) canonical.erase(it);
  };
  // Installs w as relator k (new or replacing), dropping it when trivial or a repeat.
  auto install = [&](std::size_t k, Word w) {
    w = cyclic_reduce(w);
    if (w.empty()) return;
    Word c = canonical_cyclic(w);
    if (canonical.count(c)) return;
    if (k == rel.size()) rel.emplace_back(), alive.push_back(0), form.emplace_back();
    rel[k] = std::move(w);
    alive[k] = 1;
    total += rel[k].size();
    form[k] = c;
    canonical[std::move(c)] = k;
    for (int x : rel[k]) holders[static_cast<std::size_t>(std::abs(x))].insert(k);
    if (single(rel[k]) != SIZE_MAX) candidates.insert({rel[k].size(), k});
  };
  for (const auto& r : P.relators) install(rel.size(), r);
  const std::size_t budget = 4 * total + 256;
  std::vector<char> live(P.generators, 1);
  while (!candidates.empty()) {
    // The shortest relator in which some generator occurs exactly once: r = u x^e v, so
    // x = (v u)^{-e}.
    const std::size_t best = candidates.begin()->second;
    const Word r = rel[best];
    const std::size_t pos = single(r);
    const int x = r[pos];
    const std::size_t g = static_cast<std::size_t>(std::abs(x));
    Word vu(r.begin() + static_cast<long>(pos) + 1, r.end());
    vu.insert(vu.end(), r.begin(), r.begin() + static_cast<long>(pos));
    const Word image = x > 0 ? inverse(vu) : vu;  // the word equal to x_g
    const Word image_inverse = inverse(image);
    // The substituted relators, and the growth they would bring.
    std::vector<std::pair<std::size_t, Word>> changed;
    std::size_t growth = 0;
    for (std::size_t k : holders[g]) {
      if (k == best) continue;
      Word w;
      for (int y : rel[k]) {
        if (static_cast<std::size_t>(std::abs(y)) != g) w.push_back(y);
        else w.insert(w.end(), (y > 0 ? image : image_inverse).begin(), (y > 0 ? image : image_inverse).end());
      }
      growth += w.size() > rel[k].size() ? w.size() - rel[k].size() : 0;
      changed.emplace_back(k, std::move(w));
    }
    if (total + growth > budget) break;  // substitution would grow the presentation too much
    remove(best);
    for (auto& [k, w] : changed) {
      remove(k);
      install(k, std::move(w));
    }
    live[g - 1] = 0;
  }
  // Renumber the live generators.
  std::vector<int> renumber(P.generators, 0);
  Simplified S;
  for (std::size_t i = 0; i < P.generators; i++)
    if (live[i]) {
      S.kept.push_back(i);
      renumber[i] = static_cast<int>(S.kept.size());
      S.presentation.names.push_back(P.name(i));
    }
  S.presentation.generators = S.kept.size();
  for (std::size_t k = 0; k < rel.size(); k++) {
    if (!alive[k]) continue;
    Word w;
    for (int y : rel[k]) w.push_back(y > 0 ? renumber[static_cast<std::size_t>(y) - 1] : -renumber[static_cast<std::size_t>(-y) - 1]);
    S.presentation.relators.push_back(w);
  }
  return S;
}

struct AbelianGroup {
  std::size_t rank = 0;
  std::vector<Integer> torsion;  // invariant factors > 1
  std::string str() const {
    std::string s;
    if (rank) s = rank == 1 ? "ℤ" : "ℤ^" + std::to_string(rank);
    for (const auto& t : torsion) s += (s.empty() ? "" : " ⊕ ") + std::string("ℤ/") + t.str();
    return s.empty() ? "0" : s;
  }
  friend bool operator==(const AbelianGroup&, const AbelianGroup&) = default;
};

// The matrix of exponent sums (rows: generators; columns: relators) presents G^ab as
// ℤ^n / columns, the cokernel of a map ℤ^m → ℤ^n: the H_0 of the two-term chain complex it
// forms, computed by the sparse reduction of chain_complex.hpp before any dense Smith form,
// so that presentations with thousands of generators (those of large complexes) cost little.
inline AbelianGroup abelianisation(const Presentation& P) {
  ChainComplex C;
  C.n = {P.generators, P.relators.size()};
  C.boundary.resize(2);
  C.boundary[1] = SparseMatrix(P.generators, P.relators.size());
  for (std::size_t j = 0; j < P.relators.size(); j++)
    for (int x : P.relators[j]) C.boundary[1].columns[j] = axpy(C.boundary[1].columns[j], x > 0 ? 1 : -1, Vector{{static_cast<std::size_t>(std::abs(x)) - 1, 1}});
  const auto H = homology(C, false);
  return {H[0].rank, H[0].torsion};
}

}  // namespace morphe::algebra
