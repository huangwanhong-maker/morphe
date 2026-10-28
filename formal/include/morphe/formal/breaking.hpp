// Breaking maps: letter-to-letter maps between alphabets, the images of grammars under them, the
// breakings of a torsor, and parsing from a coarse alphabet to a fine one.
//
// A map β: Γ → Σ of alphabets extends to the monoid homomorphism β: Γ* → Σ* letter by letter.
// The image β(Γ) of a context-free grammar, its terminals renamed, generates exactly β(L(Γ))
// (the context-free languages are closed under homomorphism; Hopcroft and Ullman 1979,
// Theorem 6.2, whose proof is this construction), and its derivations are the images of Γ's.
// This makes precise the two phases of a derivation in the GE-CSG draft (Huang 2026, §3.3):
// a grammar over Γ (phase I) and a breaking β (phase II), with L = β(L(Γ)).
//
// Breaking a torsor. When G acts on Σ freely and transitively (a G-torsor; the nucleotides under
// the Klein group, domains/genomics/grammar), each base point σ₀ ∈ Σ gives the bijection
// β_σ₀: G → Σ, g ↦ g·σ₀, which commutes with the left actions: β(hg) = h·β(g). The breakings
// are exactly these |Σ| = |G| maps, and two of them differ by a right translation of G,
// β_{g₀σ₀}(g) = β_σ₀(g g₀): the set of breakings is itself a torsor (the "vacuum manifold" of
// the draft, §6.2, which is G/H for a transitive action with stabiliser H), and a grammar's
// languages under two breakings correspond by the right translation, an automorphism of the
// G-set G, not by an outer automorphism of G.
//
// Coarse to fine. A chain of surjections Σ → Σ₁ → Σ₂ → … (for instance the orbit maps of a
// descending chain of subgroups, the draft's graduated breaking, §6.3) gives coarser grammars
// β_k(Γ). If S ⇒* w in Γ uses nonterminal X over the span (i, j), then the image derivation uses
// X over (i, j) in β_k(Γ) on β_k(w). So the spans on which some nonterminal has both an inside
// and an outside derivation in the coarse grammar contain every span of every complete parse in
// the fine one, and the fine inside algorithm restricted to them gives every weight of the root
// exactly, in any semiring: pruning by the coarse pass loses nothing (Charniak et al. 2006,
// Multilevel coarse-to-fine PCFG parsing, NAACL, prune approximately; here the prune is exact
// because the coarse grammar is an image, not an approximation).

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <set>
#include <stdexcept>
#include <vector>

#include <morphe/algebra/action.hpp>

#include "grammar.hpp"
#include "inside_outside.hpp"
#include "semiring.hpp"

namespace morphe::formal {

struct AlphabetMap {
  std::size_t from = 0, to = 0;
  std::vector<Symbol> image;  // image[a] ∈ [0, to)

  Symbol operator()(Symbol a) const { return image.at(a); }
  std::vector<Symbol> operator()(const std::vector<Symbol>& w) const {
    std::vector<Symbol> out;
    out.reserve(w.size());
    for (Symbol a : w) out.push_back(image.at(a));
    return out;
  }
  AlphabetMap then(const AlphabetMap& next) const {
    if (next.from != to) throw std::invalid_argument("morphe::formal::AlphabetMap::then: the alphabets do not meet");
    AlphabetMap c{from, next.to, {}};
    for (Symbol a : image) c.image.push_back(next.image.at(a));
    return c;
  }
  void validate() const {
    if (image.size() != from) throw std::invalid_argument("morphe::formal::AlphabetMap: one image per letter is needed");
    for (Symbol a : image)
      if (a >= to) throw std::invalid_argument("morphe::formal::AlphabetMap: an image outside the target alphabet");
  }
};

// The orbit map Σ → Σ/G of an action, orbits numbered by their least letter.
inline AlphabetMap orbit_map(const algebra::GroupAction& a) {
  AlphabetMap m{a.degree, a.orbit_count(), {}};
  for (auto id : a.orbit_index()) m.image.push_back(id);
  return m;
}

// The breaking of a torsor at a base point: g ↦ g·σ₀, Γ = G (element numbers) to Σ.
inline AlphabetMap torsor_breaking(const algebra::GroupAction& a, std::uint32_t base) {
  if (a.degree != a.order()) throw std::invalid_argument("morphe::formal::torsor_breaking: the action is not on a set of |G| points");
  AlphabetMap m{a.order(), a.degree, {}};
  for (std::uint32_t g = 0; g < a.order(); g++) m.image.push_back(a.act(g, base));
  std::set<Symbol> distinct(m.image.begin(), m.image.end());
  if (distinct.size() != a.degree) throw std::invalid_argument("morphe::formal::torsor_breaking: the action is not free and transitive");
  return m;
}

// The image of a context-free grammar: its terminals renamed. Rules that become equal are kept
// once (a weighted grammar should add their weights instead: image_weights does).
inline Grammar image(const Grammar& gr, const AlphabetMap& m) {
  m.validate();
  if (m.from != gr.terminals) throw std::invalid_argument("morphe::formal::image: the map is not on the grammar's terminals");
  if (!gr.context_free()) throw std::invalid_argument("morphe::formal::image: the grammar is not context-free");
  Grammar out;
  out.terminals = m.to;
  out.nonterminals = gr.nonterminals;
  const auto map = [&](const Form& f) {
    Form r;
    for (Symbol s : f) r.push_back(gr.terminal(s) ? m(s) : static_cast<Symbol>(s - gr.terminals + m.to));
    return r;
  };
  for (const auto& r : gr.rules) out.rules.push_back({map(r.lhs), map(r.rhs)});
  return out;
}

struct CoarseToFine {
  std::size_t spans = 0, kept = 0;  // spans of the word, and those the coarse pass keeps open
};

// The mask of spans a coarse pass keeps: those over which some nonterminal of the coarse grammar's
// binary form has both inside and outside derivations, on β(w).
inline std::vector<char> coarse_mask(const Grammar& fine, const AlphabetMap& beta, const std::vector<Symbol>& w, CoarseToFine* stats = nullptr) {
  const Grammar coarse = image(fine, beta);
  const BinaryGrammar b = BinaryGrammar::from(coarse);
  const std::vector<bool> ones(coarse.rules.size(), true);
  const auto cw = beta(w);
  const auto in = inside<Boolean>(b, ones, cw);
  const auto out = outside<Boolean>(b, ones, cw, in);
  const std::size_t n = w.size();
  std::vector<char> mask((n + 1) * (n + 1), 0);
  CoarseToFine s;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j <= n; j++) {
      s.spans++;
      for (std::size_t x = 0; x < b.nonterminals; x++)
        if (in.at(x, i, j) && out.at(x, i, j)) { mask[i * (n + 1) + j] = 1; s.kept++; break; }
    }
  if (stats) *stats = s;
  return mask;
}

// A derivation driven by the caller: the leftmost nonterminal is rewritten by the rule choose
// picks among those of its left side (it is given their indices), until the form is a word or
// max_steps rewrites have been made (then std::runtime_error). The generation mode of the draft
// (§5.3.3): the caller's choices may be random, structured, or both.
inline std::vector<Symbol> derive(const Grammar& gr, const std::function<std::size_t(Symbol, const std::vector<std::uint32_t>&)>& choose,
                                  std::size_t max_steps = 100000) {
  if (!gr.context_free()) throw std::invalid_argument("morphe::formal::derive: the grammar is not context-free");
  std::vector<std::vector<std::uint32_t>> by_lhs(gr.nonterminals.size());
  for (std::uint32_t r = 0; r < gr.rules.size(); r++) by_lhs[gr.index(gr.rules[r].lhs[0])].push_back(r);
  Form f{gr.start()};
  for (std::size_t step = 0; step < max_steps; step++) {
    std::size_t at = f.size();
    for (std::size_t i = 0; i < f.size(); i++)
      if (!gr.terminal(f[i])) { at = i; break; }
    if (at == f.size()) return f;
    const auto& options = by_lhs[gr.index(f[at])];
    if (options.empty()) throw std::runtime_error("morphe::formal::derive: a nonterminal without rules");
    const std::uint32_t r = options.at(choose(f[at], options));
    Form next(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(at));
    next.insert(next.end(), gr.rules[r].rhs.begin(), gr.rules[r].rhs.end());
    next.insert(next.end(), f.begin() + static_cast<std::ptrdiff_t>(at + 1), f.end());
    f = std::move(next);
  }
  throw std::runtime_error("morphe::formal::derive: the derivation did not end");
}

}  // namespace morphe::formal
