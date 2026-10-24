// The worldline complex of a record of encounters, and its refinement.
//
// A record of encounters is a sequence of events in time, with the participants present at each,
// and pairs of participants with the events at which they meet (a subset of the events both
// attend, in order). From it, with no further choice:
//
//   the worldline complex W: its vertices the events; its edges the segments, each participant
//     between two consecutive events it attends; its 2-cells the lenses, one for each pair and each
//     two consecutive meetings p, q, attached along a[p, q] · b[p, q]⁻¹, the two worldlines between
//     the meetings, pinched at them. A 2-cell's cyclic order is time.
//   the refined complex Ŵ: its vertices the presences (x, v), a participant at an event; its edges
//     first the temporal edges, one for each segment, then the rungs, one for each pair at each
//     meeting, from (a, v) to (b, v); its 2-cells first the lenses as squares
//     a[p, q] · r(q) · b[p, q]⁻¹ · r(p)⁻¹, then the triangles of rungs at an event, one for each
//     three pairs meeting there that join three participants pairwise, attached x → y → z → x.
//
// Collapsing the rungs of each event of Ŵ gives back W. The construction is that of the
// infrastructure's history complex (its theory note §2), whose occasions, beings and ties are
// the events, participants and pairs; the economics domain's spacetime complex of transactions is
// another. Each pair's meetings must be events both attend, in increasing order; a participant
// attends an event at most once. O(cells) but for the triangles, O(Σ_v m_v³) for m_v pairs meeting
// at event v.

#pragma once

#include <morphe/topology/cw_complex.hpp>

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace morphe::topology {

struct Encounters {
  std::size_t participants = 0;
  std::vector<std::vector<std::size_t>> present;  // event (in time) -> the participants there, in order
  struct Pair {
    std::size_t a = 0, b = 0;
    std::vector<std::size_t> meetings;  // events, increasing, each attended by both
  };
  std::vector<Pair> pairs;
};

struct WorldlineComplex {
  struct Segment {
    std::size_t participant = 0, from = 0, to = 0;  // events
  };
  struct Lens {
    std::size_t pair = 0, from = 0, to = 0;  // events: the meetings k and k + 1
    std::size_t k = 0;
  };
  CWComplex cw;                                       // vertex v = event v; edge e = segments[e]; 2-cell f = lenses[f]
  std::vector<std::vector<std::size_t>> worldline;    // participant -> its events, in time
  std::vector<Segment> segments;
  std::vector<std::vector<std::size_t>> segments_of;  // participant -> its segments, in time
  std::vector<Lens> lenses;
  std::vector<std::vector<std::size_t>> lenses_of;    // pair -> its lenses, in time

  // The segments of a participant's worldline from event p to event q (p before q, both attended).
  std::vector<std::size_t> path(std::size_t participant, std::size_t p, std::size_t q) const {
    const auto& w = worldline.at(participant);
    std::size_t i = w.size(), j = w.size();
    for (std::size_t k = 0; k < w.size(); k++) {
      if (w[k] == p) i = k;
      if (w[k] == q) j = k;
    }
    if (i == w.size() || j == w.size() || i > j) throw std::invalid_argument("morphe::topology::WorldlineComplex::path: not on the worldline in this order");
    return {segments_of[participant].begin() + static_cast<std::ptrdiff_t>(i), segments_of[participant].begin() + static_cast<std::ptrdiff_t>(j)};
  }
};

inline WorldlineComplex worldline_complex(const Encounters& R) {
  WorldlineComplex W;
  const std::size_t V = R.present.size(), P = R.participants;
  W.worldline.assign(P, {});
  for (std::size_t v = 0; v < V; v++)
    for (std::size_t x : R.present[v]) {
      if (x >= P) throw std::invalid_argument("morphe::topology::worldline_complex: a participant out of range");
      if (!W.worldline[x].empty() && W.worldline[x].back() == v) throw std::invalid_argument("morphe::topology::worldline_complex: a participant present twice at one event");
      W.worldline[x].push_back(v);
    }
  W.segments_of.assign(P, {});
  for (std::size_t x = 0; x < P; x++) {
    const auto& w = W.worldline[x];
    for (std::size_t k = 0; k + 1 < w.size(); k++) {
      W.segments_of[x].push_back(W.segments.size());
      W.segments.push_back({x, w[k], w[k + 1]});
    }
  }
  W.lenses_of.assign(R.pairs.size(), {});
  for (std::size_t t = 0; t < R.pairs.size(); t++) {
    const auto& m = R.pairs[t].meetings;
    for (std::size_t k = 0; k + 1 < m.size(); k++) {
      if (m[k] >= m[k + 1]) throw std::invalid_argument("morphe::topology::worldline_complex: meetings out of order");
      W.lenses_of[t].push_back(W.lenses.size());
      W.lenses.push_back({t, m[k], m[k + 1], k});
    }
  }
  auto& X = W.cw;
  X.vertices = V;
  for (const auto& s : W.segments) X.edges.push_back({s.from, s.to});
  for (const auto& L : W.lenses) {
    const auto& pr = R.pairs[L.pair];
    Face f;
    f.base = L.from;
    for (std::size_t e : W.path(pr.a, L.from, L.to)) f.boundary.push_back({e, 1});
    const auto back = W.path(pr.b, L.from, L.to);
    for (auto it = back.rbegin(); it != back.rend(); ++it) f.boundary.push_back({*it, -1});
    X.faces.push_back(std::move(f));
  }
  return W;
}

struct RefinedWorldlines {
  struct Triangle {
    std::size_t event = 0;
    std::array<std::size_t, 3> participants{};  // x < y < z
    std::array<std::size_t, 3> pairs{};         // x–y, y–z, x–z
    std::array<std::size_t, 3> meetings{};      // the index of each pair's meeting at the event
  };
  CWComplex cw;
  std::vector<std::array<std::size_t, 2>> presences;  // vertex -> (participant, event)
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> presence_of;
  std::size_t temporal = 0;                           // edges [0, temporal) are the segments
  std::vector<std::vector<std::size_t>> rung_edge;    // pair -> meeting index -> edge
  std::vector<std::array<std::size_t, 2>> rung_of_edge;  // edge - temporal -> (pair, meeting index)
  std::vector<Triangle> triangles;                    // 2-cells [lenses, lenses + triangles)
};

inline RefinedWorldlines refined_worldlines(const Encounters& R, const WorldlineComplex& W) {
  RefinedWorldlines F;
  const std::size_t V = R.present.size();
  for (std::size_t v = 0; v < V; v++)
    for (std::size_t x : R.present[v]) {
      F.presence_of[{x, v}] = F.presences.size();
      F.presences.push_back({x, v});
    }
  auto& X = F.cw;
  X.vertices = F.presences.size();
  for (const auto& s : W.segments) X.edges.push_back({F.presence_of.at({s.participant, s.from}), F.presence_of.at({s.participant, s.to})});
  F.temporal = X.edges.size();
  F.rung_edge.assign(R.pairs.size(), {});
  for (std::size_t t = 0; t < R.pairs.size(); t++) {
    const auto& pr = R.pairs[t];
    for (std::size_t k = 0; k < pr.meetings.size(); k++) {
      const std::size_t v = pr.meetings[k];
      F.rung_edge[t].push_back(X.edges.size());
      F.rung_of_edge.push_back({t, k});
      X.edges.push_back({F.presence_of.at({pr.a, v}), F.presence_of.at({pr.b, v})});
    }
  }
  // The lenses as squares.
  for (const auto& L : W.lenses) {
    const auto& pr = R.pairs[L.pair];
    Face f;
    f.base = F.presence_of.at({pr.a, L.from});
    for (std::size_t e : W.path(pr.a, L.from, L.to)) f.boundary.push_back({e, 1});
    f.boundary.push_back({F.rung_edge[L.pair][L.k + 1], 1});
    const auto back = W.path(pr.b, L.from, L.to);
    for (auto it = back.rbegin(); it != back.rend(); ++it) f.boundary.push_back({*it, -1});
    f.boundary.push_back({F.rung_edge[L.pair][L.k], -1});
    X.faces.push_back(std::move(f));
  }
  // The triangles of rungs at each event: every three pairs meeting there that join three
  // participants pairwise.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> at(V);  // event -> (pair, meeting index)
  for (std::size_t t = 0; t < R.pairs.size(); t++)
    for (std::size_t k = 0; k < R.pairs[t].meetings.size(); k++) at[R.pairs[t].meetings[k]].push_back({t, k});
  for (std::size_t v = 0; v < V; v++) {
    const auto& here = at[v];
    for (std::size_t i = 0; i < here.size(); i++)
      for (std::size_t j = i + 1; j < here.size(); j++)
        for (std::size_t k = j + 1; k < here.size(); k++) {
          const std::array<std::pair<std::size_t, std::size_t>, 3> three{here[i], here[j], here[k]};
          std::map<std::size_t, int> degree;
          for (const auto& [t, g] : three) {
            degree[R.pairs[t].a]++;
            degree[R.pairs[t].b]++;
          }
          if (degree.size() != 3) continue;
          bool cycle = true;
          for (const auto& d : degree) cycle = cycle && d.second == 2;
          if (!cycle) continue;
          RefinedWorldlines::Triangle T;
          T.event = v;
          std::size_t q = 0;
          for (const auto& d : degree) T.participants[q++] = d.first;
          const auto joins = [&](std::size_t t, std::size_t u, std::size_t w) {
            return (R.pairs[t].a == u && R.pairs[t].b == w) || (R.pairs[t].a == w && R.pairs[t].b == u);
          };
          const std::array<std::array<std::size_t, 2>, 3> sides{{{T.participants[0], T.participants[1]}, {T.participants[1], T.participants[2]}, {T.participants[0], T.participants[2]}}};
          for (std::size_t n = 0; n < 3; n++)
            for (const auto& [t, g] : three)
              if (joins(t, sides[n][0], sides[n][1])) { T.pairs[n] = t; T.meetings[n] = g; }
          // The boundary x → y → z → x.
          Face f;
          f.base = F.presence_of.at({T.participants[0], v});
          const std::array<std::array<std::size_t, 2>, 3> walk{{{T.participants[0], T.participants[1]}, {T.participants[1], T.participants[2]}, {T.participants[2], T.participants[0]}}};
          for (std::size_t n = 0; n < 3; n++) {
            const std::size_t t = T.pairs[n];
            f.boundary.push_back({F.rung_edge[t][T.meetings[n]], R.pairs[t].a == walk[n][0] ? 1 : -1});
          }
          X.faces.push_back(std::move(f));
          F.triangles.push_back(T);
        }
  }
  return F;
}

}  // namespace morphe::topology
