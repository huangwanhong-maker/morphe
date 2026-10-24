// Checks of morphe/formal/polytope.hpp: hulls in canonical form (interior and collinear points
// dropped, segments and points kept); the semiring laws of the polytope algebra on sample
// polygons; the Minkowski sum by merging edges against the hull of all pairwise sums; evaluation
// of the support function a homomorphism to (max, +) at several directions; and, on short random
// words in every mode, the alignment polygon of polytope propagation equal to the hull of the
// vectors of every enumerated alignment, every vertex attained by its witness and optimal exactly
// in its normal cone, the support function at random directions equal to the tropical optimum
// with those weights; for global alignment the parity of the gaps, the vertex (0, m + n) of all
// gaps and the vertex (0, m + n − 2 LCS) of no mismatches.

#include <morphe/check.hpp>
#include <morphe/formal/alignment.hpp>
#include <morphe/formal/polytope.hpp>
#include <morphe/formal/semiring.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace mf = morphe::formal;
using mf::LatticePoint;
using mf::Polygon;

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
  std::int64_t in(std::int64_t lo, std::int64_t hi) { return lo + static_cast<std::int64_t>(below(static_cast<std::size_t>(hi - lo + 1))); }
};

std::vector<std::uint8_t> word(Random& r, std::size_t n, std::size_t sigma) {
  std::vector<std::uint8_t> w(n);
  for (auto& c : w) c = static_cast<std::uint8_t>(r.below(sigma));
  return w;
}

Polygon random_polygon(Random& r, std::size_t points, std::int64_t range) {
  std::vector<LatticePoint> p;
  for (std::size_t k = 0; k < points; k++) p.push_back({r.in(-range, range), r.in(-range, range)});
  return Polygon::hull(p);
}

std::size_t lcs(const std::vector<std::uint8_t>& x, const std::vector<std::uint8_t>& y) {
  std::vector<std::vector<std::size_t>> L(x.size() + 1, std::vector<std::size_t>(y.size() + 1, 0));
  for (std::size_t i = 1; i <= x.size(); i++)
    for (std::size_t j = 1; j <= y.size(); j++) L[i][j] = x[i - 1] == y[j - 1] ? L[i - 1][j - 1] + 1 : std::max(L[i - 1][j], L[i][j - 1]);
  return L[x.size()][y.size()];
}

}  // namespace

int main() {
  morphe::Suite s;
  Random r{0x706f6c79};

  // ---- hulls ------------------------------------------------------------------------------------------
  {
    const Polygon sq = Polygon::hull({{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}, {1, 0}, {2, 1}, {0, 0}});
    s.equal("hull: the square with interior and collinear points is its four corners, counterclockwise from the least", sq.str(), std::string("{(0, 0), (2, 0), (2, 2), (0, 2)}"));
    s.equal("hull: collinear points are a segment, its two ends", Polygon::hull({{3, 3}, {1, 1}, {2, 2}, {0, 0}}).str(), std::string("{(0, 0), (3, 3)}"));
    s.equal("hull: one point repeated is a point", Polygon::hull({{5, -1}, {5, -1}}).str(), std::string("{(5, -1)}"));
    s.require("hull: of nothing, the empty polygon", Polygon::hull({}).empty());
    s.require("hull: the square contains its centre and its edges' points, not (3, 1)", sq.contains({1, 1}) && sq.contains({2, 1}) && !sq.contains({3, 1}));
  }

  // ---- the polytope algebra --------------------------------------------------------------------------------
  std::vector<Polygon> samples = {Polygon{}, mf::PolytopeAlgebra::one(), Polygon::point({2, -1}), Polygon::hull({{0, 0}, {3, 1}})};
  for (int k = 0; k < 4; k++) samples.push_back(random_polygon(r, 2 + r.below(7), 5));
  {
    const auto v = mf::violated_law<mf::PolytopeAlgebra>(samples);
    s.require("polytope algebra: (polygons, conv ∪, Minkowski +) satisfies every semiring law, idempotent, on 8 samples", !v, v ? *v : "every law holds");
  }
  {
    bool ok = true;
    for (int rep = 0; rep < 300; rep++) {
      const Polygon P = random_polygon(r, 1 + r.below(9), 20), Q = random_polygon(r, 1 + r.below(9), 20);
      std::vector<LatticePoint> sums;
      for (const auto& p : P.vertices())
        for (const auto& q : Q.vertices()) sums.push_back(p + q);
      ok = ok && mf::minkowski_sum(P, Q) == Polygon::hull(sums);
    }
    s.require("polytope algebra: the Minkowski sum by merged edges is the hull of all pairwise sums, on 300 random pairs", ok);
  }
  {
    bool ok = true;
    for (LatticePoint w : {LatticePoint{1, 0}, LatticePoint{-3, 2}, LatticePoint{0, -1}, LatticePoint{7, 5}}) {
      const auto v = mf::violated_morphism_law(mf::Evaluate{w}, samples);
      ok = ok && !v;
    }
    s.require("evaluation: P ↦ h_P(w) is a homomorphism to (max, +) at w = (1, 0), (−3, 2), (0, −1), (7, 5)", ok);
  }

  // ---- alignment polygons ------------------------------------------------------------------------------------
  const mf::Mode modes[] = {mf::Mode::global, mf::Mode::local, mf::Mode::semiglobal, mf::Mode::overlap};
  const char* names[] = {"global", "local", "semiglobal", "overlap"};
  const auto same = [](std::uint8_t a, std::uint8_t b) { return a == b; };
  for (int k = 0; k < 4; k++) {
    const mf::Mode mode = modes[k];
    bool hull = true, witness = true, cones = true, support = true;
    std::size_t pairs = 0, vertices = 0;
    for (int rep = 0; rep < 60; rep++) {
      const auto x = word(r, r.below(7), 2 + r.below(3)), y = word(r, r.below(7), 2 + r.below(3));
      const auto A = mf::alignment_polygon(x, y, same, mode);
      std::vector<LatticePoint> v;
      mf::enumerate_paths(x.size(), y.size(), mode, [&](const mf::EditPath& p) { v.push_back(mf::alignment_vector(p, x, y, same)); });
      hull = hull && A.polygon == Polygon::hull(v);
      pairs++;
      vertices += A.polygon.size();
      for (std::size_t q = 0; q < A.polygon.size(); q++) {
        const LatticePoint vertex = A.polygon.vertices()[q];
        witness = witness && mf::alignment_vector(A.witness[q], x, y, same) == vertex && mf::is_source(mode, A.witness[q].x_begin, A.witness[q].y_begin);
        // The cone of the vertex is exactly where it is optimal: at random directions, the vertex is
        // a maximiser iff the direction lies in its closed cone.
        for (int t = 0; t < 12; t++) {
          const LatticePoint w{r.in(-6, 6), r.in(-6, 6)};
          const bool optimal = mf::dot(vertex, w) == static_cast<__int128>(mf::support(A.polygon, w));
          cones = cones && optimal == A.fan[q].contains(w);
        }
      }
      for (int t = 0; t < 8; t++) {
        const LatticePoint w{r.in(-9, 9), r.in(-9, 9)};
        const auto tw = mf::linear<mf::MaxPlus<>>([&](std::uint8_t a, std::uint8_t b) { return a == b ? std::int64_t(0) : w.x; }, [&](std::uint8_t) { return w.y; },
                                                  [&](std::uint8_t) { return w.y; });
        support = support && mf::support(A.polygon, w) == mf::path_sum(x, y, tw, mode);
      }
    }
    const std::string tag = std::string("alignment polygon, ") + names[k] + ": ";
    s.require(tag + "polytope propagation gives the hull of the vectors (mismatches, gaps) of every enumerated alignment", hull,
              std::to_string(pairs) + " pairs of words of length ≤ 6, " + std::to_string(vertices) + " vertices");
    s.require(tag + "every vertex is the vector of its witness, which begins where the mode allows", witness);
    s.require(tag + "a vertex is optimal at a direction exactly when the direction is in its normal cone", cones);
    s.require(tag + "the support function at a direction is the tropical optimum with those weights (evaluation commutes with propagation)", support);
  }
  {
    bool parity = true, all_gaps = true, no_mismatch = true;
    std::string sizes;
    for (std::size_t n : {20, 60, 120}) {
      const auto x = word(r, n, 4), y = word(r, n + n / 10, 4);
      const auto A = mf::alignment_polygon(x, y, same);
      const std::int64_t total = static_cast<std::int64_t>(x.size() + y.size());
      for (const auto& v : A.polygon.vertices()) parity = parity && (v.y - total) % 2 == 0 && v.x >= 0 && v.y >= 0;
      const auto& V = A.polygon.vertices();
      all_gaps = all_gaps && std::find(V.begin(), V.end(), LatticePoint{0, total}) != V.end();
      no_mismatch = no_mismatch && V.front() == LatticePoint{0, total - 2 * static_cast<std::int64_t>(lcs(x, y))};
      sizes += (sizes.empty() ? "" : ", ") + std::to_string(V.size()) + " vertices at n = " + std::to_string(n);
    }
    s.require("global polygons: every vertex has gaps ≡ m + n (mod 2), since 2 (pairs) + gaps = m + n", parity, sizes);
    s.require("global polygons: (0, m + n), all gaps, is a vertex", all_gaps);
    s.require("global polygons: the least vertex is (0, m + n − 2 LCS), the alignment of a longest common subsequence", no_mismatch);
  }
  return s.done();
}
