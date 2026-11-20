// Checks of traces of records (traces.hpp).
//
// - Dependence: writing what another reads or writes; two readers of one item are independent.
// - A chain of dependent events has one record; n independent events have n!; a diamond (a before b and c,
//   both before d) has two, with Foata steps {a}, {b, c}, {d}.
// - On 300 random records of up to 7 events with random footprints: the count of linearisations equal to
//   the number of the n! orderings that keep every dependent pair in order (counted by enumeration); every
//   such ordering respected and every other refused; the Foata steps independent within, each event after
//   one of the step before.
// - Words over an alphabet (decision 0021): on 40 random independence relations of two to four letters, every
//   word of a length up to five classed by exchanging adjacent independent letters (breadth-first): the
//   lexicographic normal form constant on each class, distinct between classes, and the least word of its class;
//   the Foata normal form constant and distinct likewise; the projection lemma agreeing with the classes; each
//   class as large as the linearisations of its causal order. Independence from footprints; refusals.

#include <morphe/check.hpp>
#include <morphe/formal/traces.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace tr = morphe::formal::traces;

int main() {
  morphe::Suite s;
  const tr::Footprint r{{"x"}, {}}, w{{}, {"x"}}, rw{{"x"}, {"x"}};
  s.require("two readers are independent", !tr::dependent(r, r));
  s.require("a writer and a reader are dependent, either way", tr::dependent(w, r) && tr::dependent(r, w));
  s.require("two writers are dependent", tr::dependent(w, w));
  s.require("disjoint footprints are independent", !tr::dependent(tr::Footprint{{"a"}, {"b"}}, tr::Footprint{{"c"}, {"d"}}));

  {
    const auto P = tr::causal_order({rw, rw, rw, rw});
    s.equal("a chain of four has one record", tr::linearisations(P).str(), std::string("1"));
    std::vector<tr::Footprint> free;
    for (int i = 0; i < 6; i++) free.push_back({{"r"}, {"w" + std::to_string(i)}});
    s.equal("six independent events have 6! records", tr::linearisations(tr::causal_order(free)).str(), std::string("720"));
    const auto D = tr::causal_order({{{}, {"a"}}, {{"a"}, {"b"}}, {{"a"}, {"c"}}, {{"b", "c"}, {"d"}}});
    s.equal("a diamond has two", tr::linearisations(D).str(), std::string("2"));
    s.require("its Foata steps: {a}, {b, c}, {d}", D.steps.size() == 3 && D.steps[1] == std::vector<std::size_t>{1, 2} && D.covers.size() == 4);
  }

  std::mt19937_64 rng(3);
  bool counts = true, respect = true, foata = true;
  for (int trial = 0; trial < 300; trial++) {
    const std::size_t n = 1 + rng() % 7;
    std::vector<tr::Footprint> rec(n);
    for (auto& f : rec) {
      for (int k = 0; k < 2; k++) f.reads.insert(std::string(1, char('a' + rng() % 6)));
      if (rng() % 3) f.writes.insert(std::string(1, char('a' + rng() % 6)));
    }
    const auto P = tr::causal_order(rec);
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    long long good = 0;
    do {
      bool ok = true;
      for (std::size_t i = 0; i < n && ok; i++)
        for (std::size_t j = i + 1; j < n && ok; j++)
          if (order[j] < order[i] && tr::dependent(rec[order[i]], rec[order[j]])) ok = false;
      good += ok;
      respect = respect && P.respects(order) == ok;
    } while (std::next_permutation(order.begin(), order.end()));
    counts = counts && tr::linearisations(P).str() == std::to_string(good);
    for (std::size_t k = 0; k < P.steps.size(); k++)
      for (std::size_t a : P.steps[k]) {
        for (std::size_t b : P.steps[k]) foata = foata && (a == b || P.concurrent(a, b));
        if (k > 0) {
          bool after = false;
          for (std::size_t b : P.steps[k - 1]) after = after || P.before[b][a];
          foata = foata && after;
        }
      }
  }
  s.require("on 300 random records, the linearisations are the orderings that keep dependent pairs in order", counts);
  s.require("and exactly those orderings are respected", respect);
  s.require("Foata steps are independent within, each event after one of the step before", foata);

  // ---- words over an alphabet ----
  {
    std::mt19937_64 g(21);
    bool constant = true, distinct = true, least = true, foata_words = true, projection = true, counted = true;
    std::size_t words = 0, classes = 0;
    for (int trial = 0; trial < 40; trial++) {
      const std::size_t k = 2 + g() % 3, len = 1 + g() % 5;
      std::vector<std::pair<std::size_t, std::size_t>> pairs;
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t b = a + 1; b < k; b++)
          if (g() % 2) pairs.push_back({a, b});
      const auto I = tr::independence(k, pairs);
      std::vector<std::vector<std::size_t>> all;
      std::size_t total = 1;
      for (std::size_t i = 0; i < len; i++) total *= k;
      for (std::size_t code = 0; code < total; code++) {
        std::vector<std::size_t> w(len);
        std::size_t c = code;
        for (std::size_t i = 0; i < len; i++) w[i] = c % k, c /= k;
        all.push_back(w);
      }
      // Classes by exchanging adjacent independent letters, breadth-first.
      std::map<std::vector<std::size_t>, int> cls;
      std::vector<std::vector<std::vector<std::size_t>>> members;
      for (const auto& w : all) {
        if (cls.count(w)) continue;
        const int id = static_cast<int>(members.size());
        members.push_back({});
        std::vector<std::vector<std::size_t>> queue{w};
        cls[w] = id;
        for (std::size_t q = 0; q < queue.size(); q++) {
          const auto u = queue[q];
          members[id].push_back(u);
          for (std::size_t i = 0; i + 1 < len; i++)
            if (I.commute(u[i], u[i + 1])) {
              auto v = u;
              std::swap(v[i], v[i + 1]);
              if (!cls.count(v)) cls[v] = id, queue.push_back(v);
            }
        }
      }
      words += all.size(), classes += members.size();
      std::set<std::vector<std::size_t>> forms;
      std::set<std::vector<std::vector<std::size_t>>> foatas;
      for (const auto& m : members) {
        const auto lnf = tr::lexicographic_normal_form(m[0], I);
        const auto ff = tr::foata_normal_form(m[0], I);
        const auto smallest = *std::min_element(m.begin(), m.end());
        least = least && lnf == smallest;
        for (const auto& u : m) {
          constant = constant && tr::lexicographic_normal_form(u, I) == lnf;
          foata_words = foata_words && tr::foata_normal_form(u, I) == ff;
        }
        forms.insert(lnf), foatas.insert(ff);
        counted = counted && tr::linearisations(tr::causal_order(m[0], I)).str() == std::to_string(m.size());
      }
      distinct = distinct && forms.size() == members.size();
      foata_words = foata_words && foatas.size() == members.size();
      for (int t = 0; t < 300; t++) {
        const auto& u = all[g() % all.size()];
        const auto& v = t % 2 ? all[g() % all.size()] : members[cls[u]][g() % members[cls[u]].size()];
        projection = projection && tr::projections_agree(u, v, I) == (cls[u] == cls[v]) && tr::equivalent(u, v, I) == (cls[u] == cls[v]);
      }
    }
    const std::string tally = std::to_string(words) + " words in " + std::to_string(classes) + " traces";
    s.require("the lexicographic normal form is constant on each trace (" + tally + ")", constant);
    s.require("and distinct between traces", distinct);
    s.require("and the least word of its trace", least);
    s.require("the Foata normal form is constant on each trace and distinct between them", foata_words);
    s.require("the projection lemma and the normal forms agree with the traces found by exchanges", projection);
    s.require("each trace has as many words as its causal order has linearisations", counted);
  }
  {
    const std::vector<tr::Footprint> grants = {{{"A"}, {"A"}}, {{"B"}, {"B"}}, {{"A"}, {}}};
    const auto I = tr::independence(grants);
    s.require("from footprints: grants from two stocks commute; a reader and a writer of one stock do not", I.commute(0, 1) && !I.commute(0, 2) && I.commute(1, 2));
  }
  {
    bool self = false, range = false, asymmetric = false;
    try { tr::independence(3, {{1, 1}}); } catch (const std::invalid_argument&) { self = true; }
    try { tr::lexicographic_normal_form({0, 5}, tr::independence(3, {})); } catch (const std::invalid_argument&) { range = true; }
    tr::Independence bad = tr::independence(2, {});
    bad.independent[0][1] = 1;
    try { tr::causal_order(std::vector<std::size_t>{0, 1}, bad); } catch (const std::invalid_argument&) { asymmetric = true; }
    s.require("a letter commuting with itself, a letter out of range and an asymmetric independence are refused", self && range && asymmetric);
  }
  return s.done();
}
