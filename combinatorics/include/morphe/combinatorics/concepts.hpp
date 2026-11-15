// Formal contexts, their concepts and their implications: the Galois connection between objects and attributes.
//
// A formal context is a relation I between a set G of objects and a set M of attributes (Wille 1982; Ganter and Wille
// 1999, Formal Concept Analysis, ch. 1). For a set A of objects, A′ is the set of attributes they all have; for a set B of
// attributes, B′ the set of objects having them all. The two maps reverse inclusion, and A ⊆ B′ exactly when B ⊆ A′: they
// are a Galois connection, so that B ↦ B″ is a closure operator. A concept is a pair (A, B) with A′ = B and B′ = A: its
// extent A and its intent B determine each other. Ordered by inclusion of extents (equivalently, reverse inclusion of
// intents) the concepts form a complete lattice, in which the meet of two concepts has the intersection of their extents
// and the join the intersection of their intents (the basic theorem, Ganter and Wille 1999, Theorem 3).
//
//   - concepts: every concept, by Ganter's NextClosure, which lists the closed sets of attributes in lectic order, each once,
//     with one closure per candidate (Ganter 1984; Ganter 2010, Two basic algorithms in concept analysis, LNCS 5986).
//   - lattice: the concepts with their covering relation (this library's reduction of the strict order of extents), the
//     greatest concept (all objects) and the least (all attributes).
//   - stem_base: the implications B → C between sets of attributes that hold in the context (every object with B has C),
//     in their canonical basis: one implication P → P″ for each pseudo-intent P, a set that is not closed but contains the
//     closure of every pseudo-intent properly inside it. The basis is sound and complete (every implication that holds
//     follows from it), and no complete set of implications has fewer (Guigues and Duquenne 1986, Math. Sci. Humaines 95,
//     5–18). It is found by NextClosure under the closure of the implications found so far (Ganter 2010, Algorithm 2).
//   - closure_under: the least superset of a set closed under a list of implications, by forward chaining.
//
// Objects and attributes are numbered from 0; a set is a sorted list of numbers. A context with no objects has one concept
// (∅, M); one with no attributes one concept (G, ∅). Costs: each concept O(|M|·|G|·|M|); the lattice adds O(c²) comparisons
// and the reduction. Limits: `most` bounds the concepts or implications listed, and a longer list is refused, never cut.

#pragma once

#include <morphe/combinatorics/acyclic.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::combinatorics {

struct Context {
  std::size_t objects = 0, attributes = 0;
  std::vector<std::vector<char>> incidence;  // [object][attribute]: nonzero when the object has the attribute
};

using Set = std::vector<std::size_t>;

struct FormalConcept {
  Set extent, intent;
};

struct ConceptLattice {
  std::vector<FormalConcept> concepts;  // in lectic order of their intents
  std::vector<Arc> covers;              // (lower, upper): the extent of `lower` is a maximal proper subset of that of `upper`
  std::size_t top = 0, bottom = 0;      // the concept of all objects, and the concept of all attributes
};

struct Implication {
  Set premise, conclusion;  // conclusion ⊇ premise in the stem base: the premise's closure
};

namespace detail {

inline void check_context(const Context& c) {
  if (c.incidence.size() != c.objects) throw std::invalid_argument("morphe::combinatorics: the context's incidence has " + std::to_string(c.incidence.size()) + " rows for " + std::to_string(c.objects) + " objects");
  for (const auto& row : c.incidence)
    if (row.size() != c.attributes) throw std::invalid_argument("morphe::combinatorics: a row of the incidence has " + std::to_string(row.size()) + " entries for " + std::to_string(c.attributes) + " attributes");
}

inline std::vector<char> mask(const Set& s, std::size_t n, const char* what) {
  std::vector<char> m(n, 0);
  for (const auto x : s) {
    if (x >= n) throw std::invalid_argument(std::string("morphe::combinatorics: ") + what + " " + std::to_string(x) + " is out of range");
    m[x] = 1;
  }
  return m;
}

inline Set set_of(const std::vector<char>& m) {
  Set s;
  for (std::size_t i = 0; i < m.size(); i++)
    if (m[i]) s.push_back(i);
  return s;
}

// B′ as a mask of objects, B a mask of attributes.
inline std::vector<char> objects_of(const Context& c, const std::vector<char>& b) {
  std::vector<char> a(c.objects, 1);
  for (std::size_t g = 0; g < c.objects; g++)
    for (std::size_t m = 0; m < c.attributes && a[g]; m++)
      if (b[m] && !c.incidence[g][m]) a[g] = 0;
  return a;
}

// A′ as a mask of attributes, A a mask of objects.
inline std::vector<char> attributes_of(const Context& c, const std::vector<char>& a) {
  std::vector<char> b(c.attributes, 1);
  for (std::size_t g = 0; g < c.objects; g++)
    if (a[g])
      for (std::size_t m = 0; m < c.attributes; m++)
        if (!c.incidence[g][m]) b[m] = 0;
  return b;
}

// Ganter's next closed set after `a` in the lectic order, under the closure `close`; false when `a` is the last.
template <class Close>
bool next_closure(std::vector<char>& a, std::size_t n, Close&& close) {
  for (std::size_t i = n; i-- > 0;) {
    if (a[i]) {
      a[i] = 0;
      continue;
    }
    std::vector<char> candidate = a;  // a ∩ {0, …, i−1}, since the elements above i were cleared on the way down
    candidate[i] = 1;
    std::vector<char> b = close(candidate);
    bool same_below = true;
    for (std::size_t j = 0; j < i && same_below; j++)
      if (b[j] != a[j]) same_below = false;
    if (same_below) {
      a = std::move(b);
      return true;
    }
  }
  return false;
}

}  // namespace detail

// A′: the attributes every object of A has.
inline Set common_attributes(const Context& c, const Set& objects) {
  detail::check_context(c);
  return detail::set_of(detail::attributes_of(c, detail::mask(objects, c.objects, "object")));
}

// B′: the objects that have every attribute of B.
inline Set common_objects(const Context& c, const Set& attributes) {
  detail::check_context(c);
  return detail::set_of(detail::objects_of(c, detail::mask(attributes, c.attributes, "attribute")));
}

// B″, the closure of a set of attributes.
inline Set closure(const Context& c, const Set& attributes) {
  detail::check_context(c);
  return detail::set_of(detail::attributes_of(c, detail::objects_of(c, detail::mask(attributes, c.attributes, "attribute"))));
}

// Every concept, in lectic order of the intents.
inline std::vector<FormalConcept> concepts(const Context& c, std::size_t most = 1'000'000) {
  detail::check_context(c);
  const auto close = [&](const std::vector<char>& b) { return detail::attributes_of(c, detail::objects_of(c, b)); };
  std::vector<FormalConcept> out;
  std::vector<char> a = close(std::vector<char>(c.attributes, 0));
  do {
    if (out.size() == most) throw std::length_error("morphe::combinatorics::concepts: the context has more than " + std::to_string(most) + " concepts");
    out.push_back({detail::set_of(detail::objects_of(c, a)), detail::set_of(a)});
  } while (detail::next_closure(a, c.attributes, close));
  return out;
}

// The concepts with their covering relation.
inline ConceptLattice lattice(const Context& c, std::size_t most = 100'000) {
  ConceptLattice l;
  l.concepts = concepts(c, most);
  const std::size_t n = l.concepts.size();
  std::vector<Arc> below;  // (i, j) when the extent of i is a proper subset of that of j
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (i == j) continue;
      const auto& a = l.concepts[i].extent;
      const auto& b = l.concepts[j].extent;
      if (a.size() < b.size() && std::includes(b.begin(), b.end(), a.begin(), a.end())) below.push_back({i, j});
    }
  l.covers = reduction(n, below);
  for (std::size_t i = 0; i < n; i++) {
    if (l.concepts[i].extent.size() == c.objects) l.top = i;
    if (l.concepts[i].intent.size() == c.attributes) l.bottom = i;
  }
  return l;
}

// The least superset of `s` closed under the implications: whenever a premise lies inside, its conclusion is added.
inline Set closure_under(const std::vector<Implication>& implications, const Set& s, std::size_t attributes) {
  std::vector<char> m = detail::mask(s, attributes, "attribute");
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& imp : implications) {
      bool inside = true;
      for (const auto x : imp.premise) inside = inside && x < attributes && m[x];
      if (!inside) continue;
      for (const auto x : imp.conclusion)
        if (!m.at(x)) m[x] = 1, changed = true;
    }
  }
  return detail::set_of(m);
}

// Whether an implication holds in the context: every object with the premise has the conclusion.
inline bool holds(const Context& c, const Implication& imp) {
  const Set objects = common_objects(c, imp.premise);
  const Set has = common_objects(c, imp.conclusion);
  return std::includes(has.begin(), has.end(), objects.begin(), objects.end());
}

// The canonical basis of the implications that hold in the context, one for each pseudo-intent, in lectic order.
inline std::vector<Implication> stem_base(const Context& c, std::size_t most = 100'000) {
  detail::check_context(c);
  const std::size_t m = c.attributes;
  const auto context_close = [&](const std::vector<char>& b) { return detail::attributes_of(c, detail::objects_of(c, b)); };
  std::vector<Implication> basis;
  std::vector<std::vector<char>> premises, conclusions;  // as masks, for the pseudo-closure
  // The pseudo-closure under the basis so far: add the conclusion of every premise properly inside, until nothing changes.
  const auto pseudo_close = [&](const std::vector<char>& b) {
    std::vector<char> x = b;
    for (bool changed = true; changed;) {
      changed = false;
      for (std::size_t k = 0; k < premises.size(); k++) {
        bool inside = true, proper = false;
        for (std::size_t j = 0; j < m && inside; j++) {
          if (premises[k][j] && !x[j]) inside = false;
          if (!premises[k][j] && x[j]) proper = true;
        }
        if (!inside || !proper) continue;
        for (std::size_t j = 0; j < m; j++)
          if (conclusions[k][j] && !x[j]) x[j] = 1, changed = true;
      }
    }
    return x;
  };
  std::vector<char> a(m, 0);
  for (;;) {
    const auto closed = context_close(a);
    if (closed != a) {
      if (basis.size() == most) throw std::length_error("morphe::combinatorics::stem_base: more than " + std::to_string(most) + " implications");
      basis.push_back({detail::set_of(a), detail::set_of(closed)});
      premises.push_back(a);
      conclusions.push_back(closed);
    }
    if (!detail::next_closure(a, m, pseudo_close)) break;
  }
  return basis;
}

}  // namespace morphe::combinatorics
