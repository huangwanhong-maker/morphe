// Checks of the worldline complex of a record of encounters and its refinement (worldlines.hpp).
//
// Known by hand: three participants meeting pairwise at two events make three lenses pinched at
// two vertices (the theta-sphere, Betti 1 0 1) and a refined prism, also a sphere; a circle of
// encounters met once each is a graph with one loop; a worldline's path runs along its segments in
// time. On random records: both complexes valid (every word a closed path, ∂∂ = 0), their cell
// counts those of the construction, every lens word as long as its two paths (and two rungs in
// the refinement), every triangle three rungs; faulty records refused.

#include <morphe/check.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/worldlines.hpp>

#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace mt = morphe::topology;

namespace {

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

std::vector<std::size_t> ranks(const mt::CWComplex& X) {
  std::vector<std::size_t> b;
  for (const auto& h : mt::homology(X, false)) b.push_back(h.rank);
  return b;
}

bool valid(const mt::CWComplex& X) {
  try {
    X.validate();
    return true;
  } catch (const std::exception&) {
    return false;
  }
}

}  // namespace

int main() {
  morphe::Suite s;

  // The theta: a, b, c at two events, every pair meeting at both.
  {
    mt::Encounters R;
    R.participants = 3;
    R.present = {{0, 1, 2}, {0, 1, 2}};
    R.pairs = {{0, 1, {0, 1}}, {1, 2, {0, 1}}, {0, 2, {0, 1}}};
    const auto W = mt::worldline_complex(R);
    s.require("the theta: 2 events, 3 segments, 3 lenses", W.cw.vertices == 2 && W.cw.edges.size() == 3 && W.cw.faces.size() == 3);
    s.require("the theta is a sphere's worth of shells (Betti 1 0 1)", ranks(W.cw) == std::vector<std::size_t>{1, 0, 1}, text(ranks(W.cw)));
    const auto F = mt::refined_worldlines(R, W);
    s.require("its refinement: 6 presences, 3 + 6 edges, 3 lenses and 2 triangles", F.cw.vertices == 6 && F.temporal == 3 && F.cw.edges.size() == 9 && F.cw.faces.size() == 5 && F.triangles.size() == 2);
    s.require("the refined theta is a prism, a sphere (Betti 1 0 1)", ranks(F.cw) == std::vector<std::size_t>{1, 0, 1}, text(ranks(F.cw)));
  }
  // A circle of encounters, each pair meeting once: no lens, one loop.
  {
    mt::Encounters R;
    R.participants = 3;
    R.present = {{0, 1}, {1, 2}, {2, 0}};
    R.pairs = {{0, 1, {0}}, {1, 2, {1}}, {2, 0, {2}}};
    const auto W = mt::worldline_complex(R);
    s.require("a circle of encounters: 3 events, 3 segments, no lens, b₁ = 1", W.cw.faces.empty() && W.cw.edges.size() == 3 && ranks(W.cw) == std::vector<std::size_t>{1, 1}, text(ranks(W.cw)));
    s.require("a worldline's path runs along its segments in time", W.path(0, 0, 2) == std::vector<std::size_t>{W.segments_of[0][0]});
    bool refused = false;
    try {
      (void)W.path(0, 2, 0);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("... and a path against time is refused", refused);
  }
  // Random records.
  {
    std::mt19937_64 rng(20261024);
    bool all_valid = true, counts = true, words = true;
    for (int trial = 0; trial < 60; trial++) {
      mt::Encounters R;
      R.participants = 3 + rng() % 6;
      const std::size_t events = 2 + rng() % 12;
      for (std::size_t v = 0; v < events; v++) {
        std::vector<std::size_t> here;
        for (std::size_t x = 0; x < R.participants; x++)
          if (rng() % 2) here.push_back(x);
        R.present.push_back(here);
      }
      for (std::size_t a = 0; a < R.participants; a++)
        for (std::size_t b = a + 1; b < R.participants; b++) {
          if (rng() % 3 == 0) continue;
          mt::Encounters::Pair p{a, b, {}};
          for (std::size_t v = 0; v < events; v++) {
            bool ha = false, hb = false;
            for (std::size_t x : R.present[v]) { ha = ha || x == a; hb = hb || x == b; }
            if (ha && hb && rng() % 4) p.meetings.push_back(v);
          }
          R.pairs.push_back(p);
        }
      const auto W = mt::worldline_complex(R);
      const auto F = mt::refined_worldlines(R, W);
      all_valid = all_valid && valid(W.cw) && valid(F.cw);
      std::size_t presences = 0, segments = 0, meetings = 0, lenses = 0;
      for (const auto& h : R.present) presences += h.size();
      for (const auto& w : W.worldline) segments += w.empty() ? 0 : w.size() - 1;
      for (const auto& p : R.pairs) { meetings += p.meetings.size(); lenses += p.meetings.empty() ? 0 : p.meetings.size() - 1; }
      counts = counts && W.cw.vertices == events && W.cw.edges.size() == segments && W.cw.faces.size() == lenses && F.cw.vertices == presences &&
               F.cw.edges.size() == segments + meetings && F.cw.faces.size() == lenses + F.triangles.size();
      for (std::size_t f = 0; f < W.lenses.size(); f++) {
        const auto& L = W.lenses[f];
        const auto& p = R.pairs[L.pair];
        const std::size_t len = W.path(p.a, L.from, L.to).size() + W.path(p.b, L.from, L.to).size();
        words = words && W.cw.faces[f].boundary.size() == len && F.cw.faces[f].boundary.size() == len + 2;
      }
      for (std::size_t t = 0; t < F.triangles.size(); t++) words = words && F.cw.faces[W.lenses.size() + t].boundary.size() == 3;
    }
    s.require("sixty random records: both complexes valid (closed words, ∂∂ = 0)", all_valid);
    s.require("... with the cells the construction makes", counts);
    s.require("... every lens word its two paths (and two rungs refined), every triangle three rungs", words);
  }
  // Faulty records.
  {
    mt::Encounters twice;
    twice.participants = 2;
    twice.present = {{0, 0, 1}};
    bool refused = false;
    try {
      (void)mt::worldline_complex(twice);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("a participant present twice at one event is refused", refused);
    mt::Encounters backwards;
    backwards.participants = 2;
    backwards.present = {{0, 1}, {0, 1}};
    backwards.pairs = {{0, 1, {1, 0}}};
    refused = false;
    try {
      (void)mt::worldline_complex(backwards);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("meetings out of order are refused", refused);
  }
  return s.done();
}
