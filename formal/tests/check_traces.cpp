// Checks of traces of records (traces.hpp).
//
// - Dependence: writing what another reads or writes; two readers of one item are independent.
// - A chain of dependent events has one record; n independent events have n!; a diamond (a before b and c,
//   both before d) has two, with Foata steps {a}, {b, c}, {d}.
// - On 300 random records of up to 7 events with random footprints: the count of linearisations equal to
//   the number of the n! orderings that keep every dependent pair in order (counted by enumeration); every
//   such ordering respected and every other refused; the Foata steps independent within, each event after
//   one of the step before.

#include <morphe/check.hpp>
#include <morphe/formal/traces.hpp>

#include <algorithm>
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
  return s.done();
}
