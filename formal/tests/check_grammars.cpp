// Checks of morphe/formal/{grammar,earley,inside_outside,equivariant,breaking}.hpp.
//
// - Earley's parser decides known languages, with ε-rules among them (balanced brackets,
//   aⁿbⁿ, Aycock and Horspool's nullable grammar), and agrees with recognition by the inside
//   algorithm on random ε-free grammars over every word up to length 6; every parse tree spells
//   its word.
// - The number of parse trees (inside over the counting semiring) equals the number of leftmost
//   derivations counted by brute force; the inside probability and the expected rule counts
//   equal the sums over the enumerated derivations; the Viterbi parse attains the best
//   enumerated weight.
// - Group-invariant grammars (closures of random grammars under S₃ on three letters and under the
//   Klein group on four): Burnside's count is the number of rule orbits, at least |R|/|G|; the
//   equivariant Earley parser accepts exactly what Earley's does on the full rule set, its
//   items standing for the plain items exactly, and never more of them; the draft's
//   orbit-by-position scan accepts exactly π⁻¹(L(Γ/G)), a superset, strictly for S → aa | bb;
//   recognition with spans shared across their orbits agrees with recognition, sharing more
//   than by content alone.
// - Monotone grammars: aⁿbⁿcⁿ decided by the search of sentential forms; on invariant monotone
//   grammars the orbit-canonical search agrees with the plain one, visiting fewer forms.
// - Breaking: an image grammar generates the image of the language; the breakings of a torsor
//   differ by right translations; the coarse pass's mask keeps every complete parse, so the
//   masked inside values at the root equal the unmasked ones; derived words are accepted.

#include <morphe/algebra/action.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/check.hpp>
#include <morphe/formal/breaking.hpp>
#include <morphe/formal/earley.hpp>
#include <morphe/formal/equivariant.hpp>
#include <morphe/formal/grammar.hpp>
#include <morphe/formal/inside_outside.hpp>
#include <morphe/formal/semiring.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace mf = morphe::formal;
namespace ma = morphe::algebra;
using morphe::arithmetic::Integer;
using mf::Form;
using mf::Symbol;

namespace {

struct Random {
  std::uint64_t state;
  std::uint64_t next() {
    std::uint64_t z = (state += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(next() % n); }
  double unit() { return static_cast<double>(next() >> 11) / 9007199254740992.0; }
};

// Every word over t letters of length exactly n.
std::vector<std::vector<Symbol>> words(std::size_t t, std::size_t n) {
  std::vector<std::vector<Symbol>> out;
  std::vector<Symbol> w(n, 0);
  for (;;) {
    out.push_back(w);
    std::size_t i = 0;
    while (i < n && ++w[i] == t) w[i++] = 0;
    if (i == n) break;
  }
  if (n == 0) out.resize(1);
  return out;
}

// A random context-free, ε-free grammar with acyclic unary rules: nonterminal i rewrites only to
// right sides whose single nonterminal (for unary rules) is of higher index.
mf::Grammar random_cfg(Random& r, std::size_t t, std::size_t v, std::size_t rules) {
  mf::Grammar g;
  g.terminals = t;
  for (std::size_t i = 0; i < v; i++) g.add_nonterminal("N" + std::to_string(i));
  for (std::size_t k = 0; k < rules; k++) {
    const std::size_t a = k < v ? k : r.below(v);
    const std::size_t len = 1 + r.below(3);
    Form rhs;
    for (std::size_t i = 0; i < len; i++) {
      if (r.below(2) == 0) rhs.push_back(static_cast<Symbol>(r.below(t)));
      else rhs.push_back(g.nonterminal(r.below(v)));
    }
    if (rhs.size() == 1 && !g.terminal(rhs[0]) && g.index(rhs[0]) <= a) rhs[0] = static_cast<Symbol>(r.below(t));
    g.add(g.nonterminal(a), rhs);
  }
  // Every nonterminal can end in a terminal.
  for (std::size_t i = 0; i < v; i++) g.add(g.nonterminal(i), Form{static_cast<Symbol>(r.below(t))});
  return g;
}

// The number of leftmost derivations of w from a form, for ε-free grammars without unary cycles.
Integer count_leftmost(const mf::Grammar& g, Form form, const std::vector<Symbol>& w, std::size_t at) {
  // Consume the leading terminals.
  std::size_t i = 0;
  while (i < form.size() && g.terminal(form[i])) {
    if (at + i >= w.size() || w[at + i] != form[i]) return Integer(0);
    i++;
  }
  if (i == form.size()) return Integer(at + i == w.size() ? 1 : 0);
  // Each remaining symbol derives at least one letter.
  if (form.size() - i > w.size() - (at + i)) return Integer(0);
  Integer total(0);
  for (const auto& rule : g.rules) {
    if (rule.lhs[0] != form[i]) continue;
    Form next(form.begin() + static_cast<std::ptrdiff_t>(i), form.begin() + static_cast<std::ptrdiff_t>(i));
    next.insert(next.end(), rule.rhs.begin(), rule.rhs.end());
    next.insert(next.end(), form.begin() + static_cast<std::ptrdiff_t>(i + 1), form.end());
    total = total + count_leftmost(g, next, w, at + i);
  }
  return total;
}

// Every leftmost derivation of w, as its sequence of rules.
void derivations(const mf::Grammar& g, Form form, const std::vector<Symbol>& w, std::size_t at, std::vector<std::uint32_t>& used,
                 std::vector<std::vector<std::uint32_t>>& out) {
  std::size_t i = 0;
  while (i < form.size() && g.terminal(form[i])) {
    if (at + i >= w.size() || w[at + i] != form[i]) return;
    i++;
  }
  if (i == form.size()) { if (at + i == w.size()) out.push_back(used); return; }
  if (form.size() - i > w.size() - (at + i)) return;
  for (std::uint32_t r = 0; r < g.rules.size(); r++) {
    if (g.rules[r].lhs[0] != form[i]) continue;
    Form next(g.rules[r].rhs);
    next.insert(next.end(), form.begin() + static_cast<std::ptrdiff_t>(i + 1), form.end());
    used.push_back(r);
    derivations(g, next, w, at + i, used, out);
    used.pop_back();
  }
}

// S₃ acting on {0, 1, 2} by its permutations, and the Klein group acting regularly on four.
ma::GroupAction s3_on_three(const ma::FiniteGroup& s3) {
  // symmetric(3)'s elements are permutations; recover each by its action on a reference word:
  // build the images from the multiplication table restricted to the transpositions is not
  // available, so take the permutation group closure directly.
  std::vector<std::vector<std::uint32_t>> images(s3.order());
  // The table of S₃ is that of its permutation group in some order; act on cosets of a point
  // stabiliser: the elements fixing nothing in common give the natural action. Here we use the
  // left cosets of a subgroup of order 2 (a stabiliser), a transitive action of degree 3.
  std::vector<std::uint32_t> h;
  for (std::uint32_t g = 0; g < s3.order(); g++)
    if (s3.element_order(g) <= 2 && h.size() < 2 && (g == 0 || h.size() == 1)) h.push_back(g);
  std::vector<std::vector<std::uint32_t>> cosets;
  std::vector<std::uint32_t> coset_of(s3.order(), UINT32_MAX);
  for (std::uint32_t g = 0; g < s3.order(); g++) {
    if (coset_of[g] != UINT32_MAX) continue;
    std::vector<std::uint32_t> c;
    for (std::uint32_t x : h) c.push_back(s3.mul(g, x));
    for (std::uint32_t y : c) coset_of[y] = static_cast<std::uint32_t>(cosets.size());
    cosets.push_back(c);
  }
  for (std::uint32_t g = 0; g < s3.order(); g++)
    for (std::uint32_t c = 0; c < cosets.size(); c++) images[g].push_back(coset_of[s3.mul(g, cosets[c][0])]);
  return ma::GroupAction::from_images(s3, cosets.size(), images);
}

}  // namespace

int main() {
  morphe::Suite suite;
  Random rng{20261024};

  // ---- Earley on known languages ---------------------------------------------------------------
  {
    // Balanced brackets: S → ( S ) S | ε, over a = '(' and b = ')'.
    mf::Grammar g;
    g.terminals = 2;
    const Symbol S = g.add_nonterminal("S");
    g.add(S, Form{0, S, 1, S});
    g.add(S, Form{});
    const mf::Earley e(g);
    bool all = true;
    std::size_t tested = 0;
    for (std::size_t n = 0; n <= 10; n++)
      for (const auto& w : words(2, n)) {
        int depth = 0;
        bool ok = true;
        for (Symbol x : w) { depth += x == 0 ? 1 : -1; ok = ok && depth >= 0; }
        ok = ok && depth == 0;
        const auto r = e.parse(w, true);
        all = all && r.accepts == ok;
        if (r.accepts) all = all && mf::yield(g, r.tree, w) == w;
        tested++;
      }
    suite.require("Earley: balanced brackets, with S → ε, on every word up to length 10", all, std::to_string(tested) + " words");
  }
  {
    // aⁿbⁿ: S → a S b | ε.
    mf::Grammar g;
    g.terminals = 2;
    const Symbol S = g.add_nonterminal("S");
    g.add(S, Form{0, S, 1});
    g.add(S, Form{});
    const mf::Earley e(g);
    bool all = true;
    for (std::size_t n = 0; n <= 10; n++)
      for (const auto& w : words(2, n)) {
        bool ok = n % 2 == 0;
        for (std::size_t i = 0; i < n; i++) ok = ok && w[i] == (i < n / 2 ? 0u : 1u);
        all = all && e.accepts(w) == ok;
      }
    suite.require("Earley: aⁿbⁿ on every word up to length 10", all);
  }
  {
    // Aycock and Horspool's example: S → A A x, A → ε; only "x" is in the language, and a
    // completion with origin j finds it only by advancing over the nullable A.
    mf::Grammar g;
    g.terminals = 1;
    const Symbol S = g.add_nonterminal("S");
    const Symbol A = g.add_nonterminal("A");
    g.add(S, Form{A, A, 0});
    g.add(A, Form{});
    const mf::Earley e(g);
    const auto r = e.parse({0}, true);
    suite.require("Earley: S → A A x, A → ε accepts x (the nullable completion)", r.accepts && !e.accepts({}) && !e.accepts({0, 0}) && mf::yield(g, r.tree, {0}) == std::vector<Symbol>{0});
  }

  // ---- random grammars: Earley, inside, counts, probabilities --------------------------------------
  {
    bool agree = true, counts = true, trees = true;
    std::size_t words_tested = 0, accepted = 0;
    for (int trial = 0; trial < 40; trial++) {
      const mf::Grammar g = random_cfg(rng, 2, 3, 7);
      const mf::Earley e(g);
      const auto b = mf::BinaryGrammar::from(g);
      const std::vector<bool> ones(g.rules.size(), true);
      const std::vector<Integer> unit(g.rules.size(), Integer(1));
      for (std::size_t n = 1; n <= 6; n++)
        for (const auto& w : words(2, n)) {
          const auto r = e.parse(w, true);
          const auto in = mf::inside<mf::Boolean>(b, ones, w);
          agree = agree && r.accepts == static_cast<bool>(in.at(0, 0, n));
          if (r.accepts) { accepted++; trees = trees && mf::yield(g, r.tree, w) == w; }
          const auto ic = mf::inside<mf::Counting>(b, unit, w);
          counts = counts && ic.at(0, 0, n) == count_leftmost(g, Form{g.start()}, w, 0);
          words_tested++;
        }
    }
    suite.require("Earley agrees with inside recognition on 40 random grammars", agree, std::to_string(words_tested) + " words, " + std::to_string(accepted) + " accepted");
    suite.require("Earley's parse trees spell their words", trees);
    suite.require("the counting semiring's inside value is the number of leftmost derivations", counts);
  }
  {
    bool inside_ok = true, counts_ok = true, viterbi_ok = true;
    double worst = 0;
    for (int trial = 0; trial < 25; trial++) {
      const mf::Grammar g = random_cfg(rng, 2, 3, 6);
      // Probabilities: each nonterminal's rules normalised.
      std::vector<double> p(g.rules.size());
      std::vector<double> sum(g.nonterminals.size(), 0);
      for (std::size_t r = 0; r < p.size(); r++) { p[r] = 0.1 + rng.unit(); sum[g.index(g.rules[r].lhs[0])] += p[r]; }
      for (std::size_t r = 0; r < p.size(); r++) p[r] /= sum[g.index(g.rules[r].lhs[0])];
      std::vector<double> logp(p.size());
      for (std::size_t r = 0; r < p.size(); r++) logp[r] = std::log(p[r]);
      const auto b = mf::BinaryGrammar::from(g);
      for (std::size_t n = 2; n <= 5; n++)
        for (const auto& w : words(2, n)) {
          std::vector<std::vector<std::uint32_t>> ds;
          std::vector<std::uint32_t> used;
          derivations(g, Form{g.start()}, w, 0, used, ds);
          if (ds.empty()) continue;
          double z = 0, best = -1e300;
          std::vector<double> expect(g.rules.size(), 0);
          for (const auto& d : ds) {
            double pr = 1, lp = 0;
            for (auto r : d) { pr *= p[r]; lp += logp[r]; }
            z += pr;
            best = std::max(best, lp);
            for (auto r : d) expect[r] += pr;
          }
          for (auto& x : expect) x /= z;
          double total = 0;
          const auto c = mf::expected_counts<mf::Arithmetic<double>>(g, b, p, w, &total);
          inside_ok = inside_ok && std::fabs(total - z) <= 1e-12 * std::max(1.0, z);
          worst = std::max(worst, std::fabs(total - z) / z);
          for (std::size_t r = 0; r < c.size(); r++) counts_ok = counts_ok && std::fabs(c[r] - expect[r]) <= 1e-9;
          double vw = 0;
          const auto tree = mf::best_parse<mf::MaxPlus<double>>(g, b, logp, w, &vw);
          viterbi_ok = viterbi_ok && std::fabs(vw - best) <= 1e-9 && mf::yield(g, tree, w) == w;
        }
    }
    suite.require("the inside probability is the sum over the enumerated derivations", inside_ok);
    suite.check("its largest relative deviation", worst, 1e-12);
    suite.require("the expected rule counts (inside–outside) are the enumerated ones", counts_ok);
    suite.require("the Viterbi parse attains the best derivation and spells the word", viterbi_ok);
  }

  // ---- group-invariant grammars ---------------------------------------------------------------------
  {
    // The draft's scan rule, against S → aa | bb under a ↔ b.
    const auto z2 = ma::cyclic(2);
    const auto swap = ma::GroupAction::from_images(z2, 2, {{0, 1}, {1, 0}});
    mf::Grammar g;
    g.terminals = 2;
    const Symbol S = g.add_nonterminal("S");
    g.add(S, Form{0, 0});
    g.add(S, Form{1, 1});
    const mf::GrammarAction act{&swap, {}};
    mf::EquivariantEarley eq(g, act), loose(g, act, false);
    suite.require("S → aa | bb is invariant under a ↔ b, one rule orbit", mf::invariant(g, act) && eq.orbits().representatives.size() == 1);
    suite.require("the equivariant parser rejects ab; the orbit-by-position scan accepts it",
                  !eq.parse({0, 1}).accepts && eq.parse({1, 1}).accepts && loose.parse({0, 1}).accepts);
  }
  {
    const auto s3 = ma::symmetric(3);
    const auto on3 = s3_on_three(s3);
    const auto v4 = ma::direct_product(ma::cyclic(2), ma::cyclic(2));
    const auto regular = ma::regular_action(v4);
    struct Case { const ma::GroupAction* action; std::size_t t; const char* name; };
    for (const Case c : {Case{&on3, 3, "S₃ on three letters"}, Case{&regular, 4, "the Klein group on four letters"}}) {
      bool burnside = true, earley = true, plain_items = true, loose_is_quotient = true, spans = true, sharing = true;
      std::size_t over = 0, tested = 0, rules = 0, orbits = 0, eq_items = 0, plain_total = 0;
      std::size_t shared_group = 0, shared_content = 0;
      for (int trial = 0; trial < 12; trial++) {
        const mf::GrammarAction act{c.action, {}};
        const mf::Grammar g = mf::closure(random_cfg(rng, c.t, 3, 5), act);
        const auto orb = mf::rule_orbits(g, act);
        rules += g.rules.size();
        orbits += orb.representatives.size();
        burnside = burnside && orb.burnside == orb.representatives.size() && orb.representatives.size() * c.action->order() >= g.rules.size();
        mf::EquivariantEarley eq(g, act), loose(g, act, false);
        const mf::Earley plain(g);
        const auto q = mf::orbit_grammar(g, *c.action);
        const mf::Earley quotient(q);
        const auto pi = mf::orbit_map(*c.action);
        const auto trivial_group = ma::cyclic(1);
        const auto trivial = ma::GroupAction::from_images(trivial_group, c.t, {[&] { std::vector<std::uint32_t> id(c.t); for (std::uint32_t x = 0; x < c.t; x++) id[x] = x; return id; }()});
        const mf::GrammarAction none{&trivial, {}};
        for (std::size_t n = 1; n <= (c.t == 3 ? 5u : 4u); n++)
          for (const auto& w : words(c.t, n)) {
            const auto r = eq.parse(w);
            const auto p = plain.parse(w);
            earley = earley && r.accepts == p.accepts;
            plain_items = plain_items && r.plain_equivalent == p.items && r.items <= p.items;
            eq_items += r.items;
            plain_total += p.items;
            const bool l = loose.parse(w).accepts;
            loose_is_quotient = loose_is_quotient && l == quotient.accepts(pi(w));
            over += l && !p.accepts;
            const auto sg = mf::span_recognition(g, act, w);
            const auto sc = mf::span_recognition(g, none, w);
            spans = spans && sg.accepts == p.accepts && sc.accepts == p.accepts;
            sharing = sharing && sg.classes <= sc.classes;
            shared_group += sg.reused;
            shared_content += sc.reused;
            tested++;
          }
      }
      const std::string name = c.name;
      suite.require("Burnside counts the rule orbits, at least |R|/|G| of them (" + name + ")", burnside,
                    std::to_string(rules) + " rules in " + std::to_string(orbits) + " orbits");
      suite.require("the equivariant parser accepts what Earley's does (" + name + ")", earley, std::to_string(tested) + " words");
      suite.require("its items stand for the plain items exactly, and are no more (" + name + ")", plain_items,
                    std::to_string(eq_items) + " items against " + std::to_string(plain_total));
      suite.require("the orbit-by-position scan decides π⁻¹(L(Γ/G)) (" + name + ")", loose_is_quotient,
                    std::to_string(over) + " words accepted outside the language");
      suite.require("spans shared across orbits decide the language (" + name + ")", spans);
      suite.require("sharing by orbit reuses at least as much as sharing by content (" + name + ")", sharing && shared_group >= shared_content,
                    std::to_string(shared_group) + " spans reused against " + std::to_string(shared_content));
    }
  }

  // ---- monotone grammars ------------------------------------------------------------------------------
  {
    // aⁿbⁿcⁿ: S → aSBC | aBC, CB → BC, aB → ab, bB → bb, bC → bc, cC → cc (Hopcroft and Ullman
    // 1979, Example 9.5, in this standard form).
    mf::Grammar g;
    g.terminals = 3;
    const Symbol S = g.add_nonterminal("S"), B = g.add_nonterminal("B"), C = g.add_nonterminal("C");
    g.add(S, Form{0, S, B, C});
    g.add(S, Form{0, B, C});
    g.add(Form{C, B}, Form{B, C});
    g.add(Form{0, B}, Form{0, 1});
    g.add(Form{1, B}, Form{1, 1});
    g.add(Form{1, C}, Form{1, 2});
    g.add(Form{2, C}, Form{2, 2});
    bool all = g.monotone() && !g.context_free();
    for (std::size_t n = 1; n <= 7; n++)
      for (const auto& w : words(3, n)) {
        bool ok = n % 3 == 0;
        for (std::size_t i = 0; i < n && ok; i++) ok = w[i] == (3 * i) / n;
        all = all && mf::derives_monotone(g, w).derives == ok;
      }
    suite.require("the monotone search decides aⁿbⁿcⁿ on every word up to length 7", all);
  }
  {
    const auto v4 = ma::direct_product(ma::cyclic(2), ma::cyclic(2));
    const auto regular = ma::regular_action(v4);
    const mf::GrammarAction act{&regular, {}};
    bool agree = true;
    std::size_t plain_visits = 0, class_visits = 0, tested = 0;
    for (int trial = 0; trial < 8; trial++) {
      mf::Grammar g;
      g.terminals = 4;
      const Symbol S = g.add_nonterminal("S"), A = g.add_nonterminal("A");
      for (int k = 0; k < 3; k++) {
        const Symbol x = static_cast<Symbol>(rng.below(4)), y = static_cast<Symbol>(rng.below(4));
        g.add(S, Form{x, A});
        g.add(Form{x, A}, Form{x, y, A});
        g.add(A, Form{y});
        g.add(Form{A, x}, Form{y, A});
      }
      g = mf::closure(g, act);
      const auto canon = mf::form_canonicaliser(g, act);
      for (std::size_t n = 1; n <= 4; n++)
        for (const auto& w : words(4, n)) {
          const auto plain = mf::derives_monotone(g, w);
          const auto classes = mf::derives_monotone(g, w, canon);
          agree = agree && plain.derives == classes.derives;
          // Both searches explore every derivable form only when the word is not derived.
          if (!plain.derives) { plain_visits += plain.visited; class_visits += classes.visited; }
          tested++;
        }
    }
    suite.require("the orbit-canonical search of forms decides as the plain one (Klein group)", agree, std::to_string(tested) + " words");
    suite.require("it visits fewer forms, at most |G| times fewer", class_visits <= plain_visits && class_visits * 4 >= plain_visits,
                  std::to_string(class_visits) + " classes against " + std::to_string(plain_visits) + " forms");
  }

  // ---- breaking ------------------------------------------------------------------------------------------
  {
    // The image of a grammar generates the image of its language, length by length.
    bool all = true;
    for (int trial = 0; trial < 15; trial++) {
      const mf::Grammar g = random_cfg(rng, 4, 3, 6);
      mf::AlphabetMap beta{4, 2, {0, 1, 1, 0}};
      const mf::Grammar h = mf::image(g, beta);
      const mf::Earley eg(g), eh(h);
      for (std::size_t n = 1; n <= 4; n++) {
        std::set<std::vector<Symbol>> images;
        for (const auto& w : words(4, n))
          if (eg.accepts(w)) images.insert(beta(w));
        for (const auto& u : words(2, n)) all = all && eh.accepts(u) == static_cast<bool>(images.count(u));
      }
    }
    suite.require("an image grammar generates the image of the language", all);
  }
  {
    const auto v4 = ma::direct_product(ma::cyclic(2), ma::cyclic(2));
    const auto regular = ma::regular_action(v4);
    bool translations = true;
    for (std::uint32_t s0 = 0; s0 < 4; s0++)
      for (std::uint32_t g0 = 0; g0 < 4; g0++) {
        const auto b0 = mf::torsor_breaking(regular, s0);
        const auto b1 = mf::torsor_breaking(regular, regular.act(g0, s0));
        for (std::uint32_t g = 0; g < 4; g++) translations = translations && b1(g) == b0(v4.mul(g, g0));
      }
    suite.require("the breakings of a torsor at two base points differ by a right translation", translations);
  }
  {
    bool exact = true;
    std::size_t spans = 0, kept = 0;
    for (int trial = 0; trial < 20; trial++) {
      const mf::Grammar g = random_cfg(rng, 4, 3, 7);
      const mf::AlphabetMap beta{4, 2, {0, 0, 1, 1}};
      const auto b = mf::BinaryGrammar::from(g);
      const std::vector<Integer> unit(g.rules.size(), Integer(1));
      for (std::size_t n = 2; n <= 6; n++)
        for (int k = 0; k < 6; k++) {
          std::vector<Symbol> w(n);
          for (auto& x : w) x = static_cast<Symbol>(rng.below(4));
          mf::CoarseToFine s;
          const auto mask = mf::coarse_mask(g, beta, w, &s);
          const auto full = mf::inside<mf::Counting>(b, unit, w);
          const auto pruned = mf::inside<mf::Counting>(b, unit, w, &mask);
          exact = exact && full.at(0, 0, n) == pruned.at(0, 0, n);
          spans += s.spans;
          kept += s.kept;
        }
    }
    suite.require("coarse to fine: the masked inside count at the root is the full one", exact,
                  std::to_string(kept) + " of " + std::to_string(spans) + " spans kept open");
  }
  {
    // Derivations driven by a choice are words of the language.
    const mf::Grammar g = random_cfg(rng, 2, 3, 7);
    const mf::Earley e(g);
    bool ok = true;
    std::size_t made = 0;
    for (int k = 0; k < 50; k++) {
      std::size_t steps = 0;
      try {
        const auto w = mf::derive(g, [&](Symbol, const std::vector<std::uint32_t>& options) {
          steps++;
          // Prefer the terminal-only rules as the derivation grows, so that it ends.
          if (steps > 20)
            for (std::size_t i = 0; i < options.size(); i++) {
              bool terminal_only = true;
              for (Symbol s : g.rules[options[i]].rhs) terminal_only = terminal_only && g.terminal(s);
              if (terminal_only) return i;
            }
          return static_cast<std::size_t>(rng.below(options.size()));
        });
        ok = ok && e.accepts(w);
        made++;
      } catch (const std::runtime_error&) {
      }
    }
    suite.require("derived words are accepted", ok && made > 0, std::to_string(made) + " words");
  }
  return suite.done();
}
