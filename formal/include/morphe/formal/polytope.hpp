// The polytope algebra of lattice polygons, and polytope propagation over the edit graph: the
// alignment polygon of two words, its vertices with an alignment attaining each, and its normal
// fan, which is parametric alignment.
//
// The polytope algebra. Let P₂ be the set of convex lattice polygons in ℤ² (the empty set, points
// and segments included). With
//
//   P ⊕ Q = conv(P ∪ Q),   P ⊗ Q = P + Q = {p + q} (the Minkowski sum),   0 = ∅,   1 = {(0, 0)},
//
// (P₂, ⊕, ⊗) is a commutative idempotent semiring (Pachter and Sturmfels 2004, PNAS 101,
// 16138–16143, §2; Pachter and Sturmfels 2005, Algebraic Statistics for Computational Biology,
// §2.3): distributivity is conv(P ∪ Q) + R = conv((P + R) ∪ (Q + R)). It is not selective, so it
// has no traceback of its own; its sums carry every alignment's vector at once.
//
// Polytope propagation. Weight each edge of the edit graph by a point: a pair of equal letters
// (0, 0), a mismatch (1, 0), a gap column (0, 1). The weight of a path is the point
// (mismatches, gap columns) of the alignment, and path_sum in the polytope algebra (alignment.hpp,
// the same dynamic programme as every other semiring) is the alignment polygon: the convex hull
// of the vectors of all alignments the mode allows. For global alignment the matches are then
// determined, 2 matches + 2 mismatches + gaps = m + n, so this plane projection of the triple
// (matches, mismatches, gaps) loses nothing, and the score with match 1, mismatch −α and gap −β
// is (m + n)/2 − (1 + α) mismatches − (1/2 + β) gaps, an affine function of the point (theory
// note alignment.md §6).
//
// The support function and tropical evaluation. For w ∈ ℤ², h_P(w) = max_{p ∈ P} ⟨w, p⟩ (−∞ on the
// empty polygon) satisfies h_{P⊕Q} = max(h_P, h_Q), h_{P⊗Q} = h_P + h_Q, h_∅ = −∞, h_{{0}} = 0:
// evaluation at w is a homomorphism from the polytope algebra to the tropical semiring MaxPlus
// (Evaluate). Applied to the alignment polygon it gives the optimal score of the tropical
// alignment with weights ⟨w, ·⟩ on the edges, so one polygon answers every parameter at once.
//
// The normal fan. The (outer) normal cone of a vertex v is the set of w for which v maximises
// ⟨w, ·⟩; the cones of the vertices cover ℤ² and form the normal fan. For w strictly inside the
// cone of v the maximiser is v alone, so the tropical dynamic programme at such a w returns an
// alignment whose vector is v: its witness (Pachter and Sturmfels 2004, §3; the parametric
// regions of Waterman, Eggert and Lander 1992, PNAS 89, 6090–6093, and of Gusfield,
// Balasubramanian and Naor 1994, Algorithmica 12, 312–326).
//
// Exactness and limits. Points are pairs of 64-bit integers with |coordinate| < 2^62, so that a
// difference fits 64 bits and a cross product is exact in 128. Hulls are Andrew's monotone chain
// (Andrew 1979, Information Processing Letters 9, 216–219) with exact cross products, kept in a
// canonical form: the vertices counterclockwise from the lexicographically least, no three
// collinear, so that equal polygons are equal vectors. Minkowski sums merge the edges by angle.
//
// Complexity. A hull of k points is O(k log k); the Minkowski sum of polygons of p and q vertices
// O((p + q) log(p + q)), a translation O(p). A convex lattice polygon inside an N × N box has
// O(N^{2/3}) vertices (Andrews 1963, Trans. Amer. Math. Soc. 106, 270–279; the constant
// 12/(4π²)^{1/3} of Acketa and Žunić 1995, J. Combin. Theory A 69, 358–368), so every polygon of
// the propagation has V = O((m + n)^{2/3}) vertices, the bound of Gusfield et al. 1994 on the
// number of parametric regions; propagation costs O(mn V log V) time and O(n V) memory, and the
// witnesses O(V mn) more.

#pragma once

#include <morphe/formal/alignment.hpp>
#include <morphe/formal/semiring.hpp>

#include <algorithm>
#include <compare>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace morphe::formal {

struct LatticePoint {
  std::int64_t x = 0, y = 0;
  friend auto operator<=>(const LatticePoint&, const LatticePoint&) = default;
  friend LatticePoint operator+(LatticePoint a, LatticePoint b) { return {a.x + b.x, a.y + b.y}; }
  friend LatticePoint operator-(LatticePoint a, LatticePoint b) { return {a.x - b.x, a.y - b.y}; }
  std::string str() const { return "(" + std::to_string(x) + ", " + std::to_string(y) + ")"; }
};

// The cross product a × b and the scalar product ⟨a, b⟩, exactly.
inline __int128 cross(LatticePoint a, LatticePoint b) { return static_cast<__int128>(a.x) * b.y - static_cast<__int128>(a.y) * b.x; }
inline __int128 dot(LatticePoint a, LatticePoint b) { return static_cast<__int128>(a.x) * b.x + static_cast<__int128>(a.y) * b.y; }

class Polygon {
 public:
  Polygon() = default;  // the empty polygon, the zero of the algebra

  static Polygon point(LatticePoint p) {
    Polygon P;
    check(p);
    P.v_.push_back(p);
    return P;
  }
  // The convex hull of the points, in canonical form.
  static Polygon hull(std::vector<LatticePoint> p) {
    for (const auto& q : p) check(q);
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    Polygon P;
    if (p.size() <= 1) {
      P.v_ = std::move(p);
      return P;
    }
    std::vector<LatticePoint> h(2 * p.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < p.size(); i++) {
      while (k >= 2 && cross(h[k - 1] - h[k - 2], p[i] - h[k - 2]) <= 0) k--;
      h[k++] = p[i];
    }
    for (std::size_t i = p.size() - 1, t = k + 1; i-- > 0;) {
      while (k >= t && cross(h[k - 1] - h[k - 2], p[i] - h[k - 2]) <= 0) k--;
      h[k++] = p[i];
    }
    h.resize(k - 1);
    P.v_ = std::move(h);
    return P;
  }

  // The vertices, counterclockwise from the lexicographically least.
  const std::vector<LatticePoint>& vertices() const { return v_; }
  std::size_t size() const { return v_.size(); }
  bool empty() const { return v_.empty(); }
  bool contains(LatticePoint q) const {
    if (v_.empty()) return false;
    if (v_.size() == 1) return q == v_[0];
    if (v_.size() == 2) return cross(v_[1] - v_[0], q - v_[0]) == 0 && dot(q - v_[0], q - v_[1]) <= 0;
    for (std::size_t k = 0; k < v_.size(); k++)
      if (cross(v_[(k + 1) % v_.size()] - v_[k], q - v_[k]) < 0) return false;
    return true;
  }
  // The polygon moved by d.
  Polygon translated(LatticePoint d) const {
    Polygon P = *this;
    for (auto& p : P.v_) {
      p = p + d;
      check(p);
    }
    return P;
  }
  friend bool operator==(const Polygon&, const Polygon&) = default;
  std::string str() const {
    std::string s = "{";
    for (std::size_t k = 0; k < v_.size(); k++) s += (k ? ", " : "") + v_[k].str();
    return s + "}";
  }

 private:
  static void check(LatticePoint p) {
    constexpr std::int64_t limit = std::int64_t(1) << 62;
    if (p.x >= limit || p.x <= -limit || p.y >= limit || p.y <= -limit) throw std::invalid_argument("morphe::formal::Polygon: a coordinate beyond 2^62");
  }
  std::vector<LatticePoint> v_;
};

// conv(P ∪ Q).
inline Polygon convex_union(const Polygon& P, const Polygon& Q) {
  if (P.empty()) return Q;
  if (Q.empty()) return P;
  std::vector<LatticePoint> p = P.vertices();
  p.insert(p.end(), Q.vertices().begin(), Q.vertices().end());
  return Polygon::hull(std::move(p));
}

// P + Q, by merging the edges of the two polygons in the order of their angles (de Berg et al.
// 2008, Computational Geometry, 3rd ed., §13.3).
inline Polygon minkowski_sum(const Polygon& P, const Polygon& Q) {
  if (P.empty() || Q.empty()) return {};
  if (P.size() == 1) return Q.translated(P.vertices()[0]);
  if (Q.size() == 1) return P.translated(Q.vertices()[0]);
  // Both counterclockwise from their lowest vertex (least y, then least x), so that the angles of
  // the edges increase through [0, 2π).
  const auto from_lowest = [](const Polygon& R) {
    std::vector<LatticePoint> v = R.vertices();
    const auto lowest = std::min_element(v.begin(), v.end(), [](LatticePoint a, LatticePoint b) { return a.y < b.y || (a.y == b.y && a.x < b.x); });
    std::rotate(v.begin(), lowest, v.end());
    return v;
  };
  std::vector<LatticePoint> a = from_lowest(P), b = from_lowest(Q), out;
  const std::size_t na = a.size(), nb = b.size();
  out.reserve(na + nb);
  std::size_t i = 0, j = 0;
  // When one polygon's edges are spent, its next edge is its first again, at angle near 0 ≡ 2π,
  // and the other's remaining edges turn clockwise from it, so the loop always advances.
  while (i < na || j < nb) {
    if (out.size() > na + nb) throw std::logic_error("morphe::formal::minkowski_sum: the merge of edges did not advance");
    out.push_back(a[i % na] + b[j % nb]);
    const __int128 c = cross(a[(i + 1) % na] - a[i % na], b[(j + 1) % nb] - b[j % nb]);
    if (c >= 0 && i < na) i++;
    if (c <= 0 && j < nb) j++;
  }
  return Polygon::hull(std::move(out));
}

struct PolytopeAlgebra {
  using value_type = Polygon;
  static constexpr bool idempotent = true, selective = false, commutative = true;
  static Polygon zero() { return {}; }
  static Polygon one() { return Polygon::point({0, 0}); }
  static Polygon plus(const Polygon& a, const Polygon& b) { return convex_union(a, b); }
  static Polygon times(const Polygon& a, const Polygon& b) { return minkowski_sum(a, b); }
  static bool equal(const Polygon& a, const Polygon& b) { return a == b; }
};

// h_P(w) = max over P of ⟨w, p⟩; MaxPlus's −∞ for the empty polygon. Throws std::overflow_error
// when the value leaves 64 bits.
inline std::int64_t support(const Polygon& P, LatticePoint w) {
  if (P.empty()) return MaxPlus<std::int64_t>::zero();
  __int128 best = dot(P.vertices()[0], w);
  for (const auto& v : P.vertices()) best = std::max(best, dot(v, w));
  if (best >= std::numeric_limits<std::int64_t>::max() || best <= std::numeric_limits<std::int64_t>::min())
    throw std::overflow_error("morphe::formal::support: beyond 64 bits");
  return static_cast<std::int64_t>(best);
}

// Evaluation at w: the homomorphism P ↦ h_P(w) from the polytope algebra to MaxPlus.
struct Evaluate {
  using source = PolytopeAlgebra;
  using target = MaxPlus<std::int64_t>;
  LatticePoint w;
  std::int64_t operator()(const Polygon& P) const { return support(P, w); }
};

// ---- the normal fan -------------------------------------------------------------------------------
// The outer normal cone of a vertex: the w for which it maximises ⟨w, ·⟩ over the polygon.
struct NormalCone {
  enum class Kind { plane, half_plane, sector };
  Kind kind = Kind::plane;
  // sector: the closed cone from `first` counterclockwise to `second` (less than a half-turn),
  // the outer normals of the incoming and the outgoing edge; half-plane: {w : ⟨w, first⟩ ≥ 0}.
  LatticePoint first, second;
  LatticePoint interior;  // a direction strictly inside
  bool contains(LatticePoint w) const {
    switch (kind) {
      case Kind::plane: return true;
      case Kind::half_plane: return dot(w, first) >= 0;
      case Kind::sector: return cross(first, w) >= 0 && cross(w, second) >= 0;
    }
    return false;
  }
  bool strictly_contains(LatticePoint w) const {
    switch (kind) {
      case Kind::plane: return true;
      case Kind::half_plane: return dot(w, first) > 0;
      case Kind::sector: return cross(first, w) > 0 && cross(w, second) > 0;
    }
    return false;
  }
};

// The normal cone of each vertex, in the order of the vertices.
inline std::vector<NormalCone> normal_fan(const Polygon& P) {
  std::vector<NormalCone> fan;
  const auto& v = P.vertices();
  if (v.size() == 1) {
    fan.push_back({NormalCone::Kind::plane, {0, 0}, {0, 0}, {0, 0}});
  } else if (v.size() == 2) {
    for (std::size_t k = 0; k < 2; k++) {
      const LatticePoint d = v[k] - v[1 - k];
      fan.push_back({NormalCone::Kind::half_plane, d, {0, 0}, d});
    }
  } else {
    const std::size_t n = v.size();
    // The outer normal of a counterclockwise edge (dx, dy) is (dy, −dx).
    const auto normal = [](LatticePoint e) { return LatticePoint{e.y, -e.x}; };
    for (std::size_t k = 0; k < n; k++) {
      const LatticePoint in = normal(v[k] - v[(k + n - 1) % n]), out = normal(v[(k + 1) % n] - v[k]);
      fan.push_back({NormalCone::Kind::sector, in, out, in + out});
    }
  }
  return fan;
}

// ---- the alignment polygon --------------------------------------------------------------------------
// The vector (mismatches, gap columns) of a path, letters being equal when same(a, b).
template <class X, class Y, class Same>
LatticePoint alignment_vector(const EditPath& p, const X& x, const Y& y, Same same) {
  LatticePoint v;
  std::size_t i = p.x_begin, j = p.y_begin;
  for (char op : p.ops) {
    if (op == 'M') {
      if (i >= x.size() || j >= y.size()) throw std::out_of_range("morphe::formal::alignment_vector: the path leaves the words");
      v.x += same(x[i], y[j]) ? 0 : 1;
      i++;
      j++;
    } else {
      v.y++;
      (op == 'I' ? i : j)++;
    }
  }
  return v;
}

struct AlignmentPolygon {
  Polygon polygon;                // the hull of the vectors (mismatches, gap columns) of the paths
  std::vector<EditPath> witness;  // witness[k] is a path whose vector is vertex k
  std::vector<NormalCone> fan;    // fan[k] is the outer normal cone of vertex k
};

// Polytope propagation: the alignment polygon of x and y in the mode, its vertices each with an
// alignment attaining it (the tropical optimum at an interior direction of its cone), and its
// normal fan.
template <class X, class Y, class Same>
AlignmentPolygon alignment_polygon(const X& x, const Y& y, Same same, Mode mode = Mode::global) {
  const Polygon mismatch = Polygon::point({1, 0}), gap = Polygon::point({0, 1}), match = PolytopeAlgebra::one();
  const auto w = linear<PolytopeAlgebra>([&](const auto& a, const auto& b) -> const Polygon& { return same(a, b) ? match : mismatch; },
                                         [&](const auto&) -> const Polygon& { return gap; }, [&](const auto&) -> const Polygon& { return gap; });
  AlignmentPolygon out;
  out.polygon = path_sum(x, y, w, mode);
  out.fan = normal_fan(out.polygon);
  for (std::size_t k = 0; k < out.polygon.size(); k++) {
    const LatticePoint d = out.fan[k].interior;
    const auto t = linear<MaxPlus<std::int64_t>>([&](const auto& a, const auto& b) { return same(a, b) ? std::int64_t(0) : d.x; },
                                                 [&](const auto&) { return d.y; }, [&](const auto&) { return d.y; });
    const auto best = optimal_path(x, y, t, mode);
    if (!best.exists || alignment_vector(best.path, x, y, same) != out.polygon.vertices()[k])
      throw std::logic_error("morphe::formal::alignment_polygon: the optimum inside a normal cone is not its vertex");
    out.witness.push_back(best.path);
  }
  return out;
}

}  // namespace morphe::formal
