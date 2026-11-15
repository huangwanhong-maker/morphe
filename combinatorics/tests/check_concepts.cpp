// Checks of formal concept analysis (concepts.hpp): the concepts of random contexts equal the closed sets of attributes
// found among every subset, listed once each in lectic order; the contranominal scale has 2ⁿ concepts and the ordinal scale
// n; the meet of any two concepts has the intersection of their extents and the join the intersection of their intents;
// the covering relation equals the one found by brute force; the stem base is sound, complete and minimal, and has one
// implication for each pseudo-intent found by its recursive definition.

#include <morphe/check.hpp>
#include <morphe/combinatorics/concepts.hpp>

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <set>
#include <string>
#include <vector>

using namespace morphe::combinatorics;

namespace {

std::uint64_t state = 0x9e3779b97f4a7c15ull;
std::uint64_t next() {
  state ^= state << 13, state ^= state >> 7, state ^= state << 17;
  return state;
}

Context random_context(std::size_t g, std::size_t m, double density) {
  Context c{g, m, std::vector<std::vector<char>>(g, std::vector<char>(m, 0))};
  for (auto& row : c.incidence)
    for (auto& x : row) x = (next() % 1000) < density * 1000;
  return c;
}

Set of_mask(std::uint32_t mask, std::size_t m) {
  Set s;
  for (std::size_t i = 0; i < m; i++)
    if (mask >> i & 1) s.push_back(i);
  return s;
}

// Lectic order: A < B when the least element of their symmetric difference lies in B.
bool lectic_less(const Set& a, const Set& b, std::size_t m) {
  std::vector<char> x(m, 0), y(m, 0);
  for (auto i : a) x[i] = 1;
  for (auto i : b) y[i] = 1;
  for (std::size_t i = 0; i < m; i++)
    if (x[i] != y[i]) return y[i];
  return false;
}

}  // namespace

int main() {
  morphe::Suite s;

  std::size_t contexts = 0, concept_faults = 0, order_faults = 0, lattice_faults = 0, cover_faults = 0, base_faults = 0, total_concepts = 0, total_implications = 0;
  for (int trial = 0; trial < 120; trial++) {
    const std::size_t g = next() % 8, m = 1 + next() % 7;
    const Context c = random_context(g, m, 0.25 + 0.5 * (next() % 100) / 100.0);
    contexts++;
    const auto cs = concepts(c);
    total_concepts += cs.size();
    // Brute force: every subset B with B″ = B.
    std::set<Set> intents, mine;
    for (std::uint32_t mask = 0; mask < (1u << m); mask++) {
      const Set b = of_mask(mask, m);
      if (closure(c, b) == b) intents.insert(b);
    }
    for (const auto& k : cs) {
      mine.insert(k.intent);
      if (common_objects(c, k.intent) != k.extent || common_attributes(c, k.extent) != k.intent) concept_faults++;
    }
    if (mine != intents || mine.size() != cs.size()) concept_faults++;
    for (std::size_t i = 1; i < cs.size(); i++)
      if (!lectic_less(cs[i - 1].intent, cs[i].intent, m)) order_faults++;

    // The lattice: meets and joins by intersection, and the covers.
    const auto l = lattice(c);
    std::set<Set> extents;
    for (const auto& k : l.concepts) extents.insert(k.extent);
    for (const auto& a : l.concepts)
      for (const auto& b : l.concepts) {
        Set meet, join;
        std::set_intersection(a.extent.begin(), a.extent.end(), b.extent.begin(), b.extent.end(), std::back_inserter(meet));
        std::set_intersection(a.intent.begin(), a.intent.end(), b.intent.begin(), b.intent.end(), std::back_inserter(join));
        if (!extents.count(meet) || !intents.count(join)) lattice_faults++;
      }
    if (l.concepts[l.top].extent.size() != g || l.concepts[l.bottom].intent.size() != m) lattice_faults++;
    std::set<std::pair<std::size_t, std::size_t>> covers, brute;
    for (const auto& a : l.covers) covers.insert({a.first, a.second});
    const auto proper = [&](std::size_t i, std::size_t j) {
      const auto& x = l.concepts[i].extent;
      const auto& y = l.concepts[j].extent;
      return x.size() < y.size() && std::includes(y.begin(), y.end(), x.begin(), x.end());
    };
    for (std::size_t i = 0; i < l.concepts.size(); i++)
      for (std::size_t j = 0; j < l.concepts.size(); j++) {
        if (!proper(i, j)) continue;
        bool between = false;
        for (std::size_t k = 0; k < l.concepts.size() && !between; k++) between = proper(i, k) && proper(k, j);
        if (!between) brute.insert({i, j});
      }
    if (covers != brute) cover_faults++;

    // The stem base: sound, complete, minimal, one implication for each pseudo-intent.
    const auto base = stem_base(c);
    total_implications += base.size();
    for (const auto& imp : base)
      if (!holds(c, imp)) base_faults++;
    const auto complete = [&](const std::vector<Implication>& b) {
      for (std::uint32_t mask = 0; mask < (1u << m); mask++) {
        const Set x = of_mask(mask, m);
        if (closure_under(b, x, m) != closure(c, x)) return false;
      }
      return true;
    };
    if (!complete(base)) base_faults++;
    for (std::size_t k = 0; k < base.size(); k++) {
      auto fewer = base;
      fewer.erase(fewer.begin() + static_cast<long>(k));
      if (complete(fewer)) base_faults++;
    }
    // Pseudo-intents by their definition, smallest first: P ≠ P″, and Q″ ⊆ P for every pseudo-intent Q ⊊ P.
    std::vector<Set> pseudo;
    for (std::size_t size = 0; size <= m; size++)
      for (std::uint32_t mask = 0; mask < (1u << m); mask++) {
        const Set p = of_mask(mask, m);
        if (p.size() != size || closure(c, p) == p) continue;
        bool ok = true;
        for (const auto& q : pseudo) {
          if (q.size() < p.size() && std::includes(p.begin(), p.end(), q.begin(), q.end())) {
            const Set qq = closure(c, q);
            if (!std::includes(p.begin(), p.end(), qq.begin(), qq.end())) ok = false;
          }
        }
        if (ok) pseudo.push_back(p);
      }
    std::set<Set> premises, defined(pseudo.begin(), pseudo.end());
    for (const auto& imp : base) premises.insert(imp.premise);
    if (premises != defined) base_faults++;
  }
  s.require("on " + std::to_string(contexts) + " random contexts the concepts are the closed sets found among every subset of attributes, each with its extent (" + std::to_string(total_concepts) + " concepts)", concept_faults == 0, std::to_string(concept_faults));
  s.require("listed in lectic order, each once", order_faults == 0);
  s.require("the meet of any two concepts has the intersection of their extents, the join the intersection of their intents", lattice_faults == 0);
  s.require("the covering relation is the one found by brute force", cover_faults == 0);
  s.require("the stem base is sound, complete and minimal, one implication for each pseudo-intent (" + std::to_string(total_implications) + " implications)", base_faults == 0, std::to_string(base_faults));

  // ---- Scales with closed forms -------------------------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 6;
    Context nominal{n, n, std::vector<std::vector<char>>(n, std::vector<char>(n, 0))};
    Context ordinal = nominal;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) nominal.incidence[i][j] = i != j, ordinal.incidence[i][j] = i <= j;
    s.require("the contranominal scale of 6 has 2⁶ = 64 concepts, a Boolean lattice, and an empty stem base", concepts(nominal).size() == 64 && stem_base(nominal).empty());
    const auto chain = lattice(ordinal);
    s.require("the ordinal scale (6, 6, ≤) has 6 concepts, a chain of 5 covers: every object has the last attribute, so ∅ is not closed", chain.concepts.size() == 6 && chain.covers.size() == 5);
    const Context none{0, 3, {}};
    const auto only = concepts(none);
    s.require("a context with no objects has one concept, (∅, M)", only.size() == 1 && only[0].extent.empty() && only[0].intent.size() == 3);
  }
  bool refused = false;
  try { (void)concepts(Context{2, 2, {{1, 0}}}); } catch (const std::invalid_argument&) { refused = true; }
  s.require("an incidence of the wrong shape is refused", refused);
  return s.done();
}
