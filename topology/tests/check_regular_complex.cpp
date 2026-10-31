// Checks of regular CW complexes held by their face posets (regular_complex.hpp).
//
// Known by hand:
// - a square 2-cell with its edges in mixed directions: regular, oriented (+, +, −, −) round its boundary
//   from a given sign, ∂∂ = 0, a point's homology, the same as a CWComplex; a conflicting sign refused;
// - cubical grids: the cell counts of the 3 × 2 × 2 grid, Euler characteristic 1, a ball's homology, the
//   product orientation satisfying ∂∂ = 0 and kept by orient(), and recovered by it from no signs at all;
//   the boundary of the cube a 2-sphere (Betti 1 0 1);
// - the filled 3-simplex is a ball and its 2-skeleton a sphere: dropping the top cell changes b₂;
// - the 3 × 3 torus is regular with Betti 1 2 1 as a regular complex and as a CWComplex;
// - the 4-simplex and its boundary: the 4-cell admitted by homology alone, marked inexact;
// - refusals, each with its reason: a loop, a 2-cell bounded by two circles or by a figure eight, a 3-cell
//   bounded by a torus, by two spheres, by two spheres pinched at a vertex, by an open cube;
// - the face poset is Eulerian: μ(σ, τ) = (−1)^{dim τ − dim σ} on every interval of the cube (Philip
//   Hall's theorem, the intervals being spheres), by morphe's Möbius function.

#include <morphe/check.hpp>
#include <morphe/combinatorics/mobius.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/regular_complex.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace mt = morphe::topology;

namespace {

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

std::vector<std::size_t> betti(const mt::RegularComplex& X) { return morphe::algebra::betti(mt::chain_complex(X)); }

// The simplex on n + 1 vertices as a regular complex (every nonempty face a cell), optionally without its
// top cell.
mt::RegularComplex simplex(std::size_t n, bool with_top = true) {
  mt::RegularComplex X;
  std::map<unsigned, std::size_t> cell;
  for (std::size_t k = 0; k <= n; k++)
    for (unsigned S = 1; S < (1u << (n + 1)); S++) {
      if (static_cast<std::size_t>(__builtin_popcount(S)) != k + 1) continue;
      if (k == n && !with_top) continue;
      std::vector<std::size_t> facets;
      if (k > 0)
        for (unsigned b = 0; b <= n; b++)
          if (S & (1u << b)) facets.push_back(cell.at(S ^ (1u << b)));
      cell[S] = X.add(k, "", facets);
    }
  mt::orient(X);
  return X;
}

// The torus from a 3 × 3 grid of squares with opposite sides identified.
mt::RegularComplex torus3() {
  mt::RegularComplex X;
  std::size_t v[3][3], h[3][3], w[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) v[i][j] = X.add_vertex();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      h[i][j] = X.add_edge(v[i][j], v[(i + 1) % 3][j]);
      w[i][j] = X.add_edge(v[i][j], v[i][(j + 1) % 3]);
    }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) X.add(2, "square", {h[i][j], w[(i + 1) % 3][j], h[i][(j + 1) % 3], w[i][j]});
  mt::orient(X);
  return X;
}

}  // namespace

int main() {
  morphe::Suite s;

  // ---- a square ----
  {
    mt::RegularComplex X;
    const auto v0 = X.add_vertex(), v1 = X.add_vertex(), v2 = X.add_vertex(), v3 = X.add_vertex();
    const auto a = X.add_edge(v0, v1), b = X.add_edge(v1, v2), c = X.add_edge(v3, v2), d = X.add_edge(v0, v3);
    X.add(2, "block", {a, b, c, d}, {1, 0, 0, 0});
    const auto adm = mt::admissible(X);
    s.require("a square 2-cell is regular", adm.regular && adm.exact, adm.reason);
    mt::orient(X);
    const auto& f = X.cells.back();
    s.require("oriented round its boundary: + + − −", f.signs == std::vector<int>{1, 1, -1, -1});
    bool dd = true;
    try { mt::chain_complex(X).validate(); } catch (...) { dd = false; }
    s.require("∂∂ = 0", dd);
    s.equal("a point's homology", text(betti(X)), std::string("1 0 0"));
    s.equal("Euler characteristic 1", X.euler_characteristic(), 1LL);
    const auto Y = mt::cw(X);
    bool valid = true;
    try { Y.validate(); } catch (...) { valid = false; }
    s.require("as a CWComplex its attaching word is a closed path", valid && Y.faces.size() == 1 && Y.faces[0].boundary.size() == 4);
    s.equal("and its homology is the same", text(mt::betti(Y)), std::string("1 0 0"));
    s.require("the ends of an edge read from its signs", X.ends(c)[0] == v3 && X.ends(c)[1] == v2);
    mt::RegularComplex Z = X;
    Z.cells.back().signs = {1, 1, 1, 0};
    bool refused = false;
    try { mt::orient(Z); } catch (const std::exception&) { refused = true; }
    s.require("signs that break ∂∂ = 0 are refused", refused);
    Z.cells.back().signs = {-1, 0, 0, 0};
    mt::orient(Z);
    s.require("from the opposite sign, the opposite orientation", Z.cells.back().signs == std::vector<int>{-1, -1, 1, 1});
  }

  // ---- cubical grids ----
  {
    const auto G = mt::cubical_grid(3, 2, 2);
    const auto& X = G.complex;
    s.equal("the 3 × 2 × 2 grid: 36 vertices, 75 edges, 52 squares, 12 cubes", text(X.counts()), std::string("36 75 52 12"));
    const auto adm = mt::admissible(X);
    s.require("is regular", adm.regular && adm.exact, adm.reason);
    bool dd = true;
    try { mt::chain_complex(X).validate(); } catch (...) { dd = false; }
    s.require("its product orientation has ∂∂ = 0", dd);
    s.equal("a ball's homology", text(betti(X)), std::string("1 0 0 0"));
    mt::RegularComplex Y = X;
    mt::orient(Y);
    bool same = true;
    for (std::size_t i = 0; i < X.size(); i++) same = same && X.cells[i].signs == Y.cells[i].signs;
    s.require("orient() keeps the signs given", same);
    mt::RegularComplex Z = X;
    for (auto& c : Z.cells) std::fill(c.signs.begin(), c.signs.end(), 0);
    mt::orient(Z);
    dd = true;
    try { mt::chain_complex(Z).validate(); } catch (...) { dd = false; }
    s.require("and finds an orientation from no signs at all", dd && text(betti(Z)) == "1 0 0 0");
    s.require("cells are found by corner and span", G.at({2, 1, 0}, {1, 1, 0}) != mt::no_cell && X.cells[G.at({2, 1, 0}, {1, 1, 0})].type == "xy" && G.at({3, 0, 0}, {1, 0, 0}) == mt::no_cell);

    auto C = mt::cubical_grid(1, 1, 1).complex;
    s.equal("the solid cube: Betti 1 0 0 0", text(betti(C)), std::string("1 0 0 0"));
    C.cells.pop_back();
    s.equal("its boundary is a 2-sphere: Betti 1 0 1", text(betti(C)), std::string("1 0 1"));
  }

  // ---- simplices ----
  {
    const auto T = simplex(3), T2 = simplex(3, false);
    s.require("the filled 3-simplex is regular", mt::admissible(T).regular);
    s.equal("and a ball", text(betti(T)), std::string("1 0 0 0"));
    s.equal("its 2-skeleton is a sphere: b₂ = 1, so dropping the 3-cell loses what it filled", text(betti(T2)), std::string("1 0 1"));
    const auto P = simplex(4), dP = simplex(4, false);
    const auto a = mt::admissible(P), b = mt::admissible(dP);
    s.require("the 4-simplex is admitted, by homology alone (inexact)", a.regular && !a.exact, a.reason);
    s.require("its boundary, of dimension 3, exactly", b.regular && b.exact, b.reason);
    s.equal("the 4-simplex: Betti 1 0 0 0 0", text(betti(P)), std::string("1 0 0 0 0"));
    s.equal("its boundary: Betti 1 0 0 1", text(betti(dP)), std::string("1 0 0 1"));
  }

  // ---- the torus ----
  {
    const auto X = torus3();
    s.require("the 3 × 3 torus is regular", mt::admissible(X).regular);
    s.equal("Betti 1 2 1 as a regular complex", text(betti(X)), std::string("1 2 1"));
    s.equal("and as a CWComplex", text(mt::betti(mt::cw(X))), std::string("1 2 1"));
  }

  // ---- refusals ----
  {
    const auto refused = [&](const std::string& name, const mt::RegularComplex& X, const std::string& word) {
      const auto a = mt::admissible(X);
      s.require(name, !a.regular && a.reason.find(word) != std::string::npos, a.reason);
    };
    {
      mt::RegularComplex X;
      const auto v = X.add_vertex();
      X.add(1, "", {v, v});
      refused("a loop is not a regular 1-cell", X, "twice");
    }
    const auto triangle = [](mt::RegularComplex& X) {
      const auto a = X.add_vertex(), b = X.add_vertex(), c = X.add_vertex();
      return std::vector<std::size_t>{X.add_edge(a, b), X.add_edge(b, c), X.add_edge(c, a)};
    };
    {
      mt::RegularComplex X;
      auto e = triangle(X), f = triangle(X);
      e.insert(e.end(), f.begin(), f.end());
      X.add(2, "", e);
      refused("a 2-cell bounded by two circles", X, "not connected");
    }
    {
      mt::RegularComplex X;
      const auto o = X.add_vertex(), a = X.add_vertex(), b = X.add_vertex(), c = X.add_vertex(), d = X.add_vertex();
      X.add(2, "", {X.add_edge(o, a), X.add_edge(a, b), X.add_edge(b, o), X.add_edge(o, c), X.add_edge(c, d), X.add_edge(d, o)});
      refused("a 2-cell bounded by a figure eight", X, "on 4 of its edges");
    }
    {
      auto X = torus3();
      X.add(3, "", X.of_dimension(2));
      refused("a 3-cell bounded by a torus", X, "Euler characteristic 0");
    }
    {
      auto A = simplex(3, false);
      const auto n = A.size();
      auto B = simplex(3, false);
      for (auto c : B.cells) {
        for (auto& f : c.facets) f += n;
        A.add(c.dim, c.type, c.facets, c.signs);
      }
      A.add(3, "", A.of_dimension(2));
      refused("a 3-cell bounded by two spheres", A, "not connected");
    }
    {
      // Two tetrahedral boundaries sharing vertex 0.
      auto A = simplex(3, false);
      const auto B = simplex(3, false);
      std::map<std::size_t, std::size_t> at{{0, 0}};
      for (std::size_t i = 0; i < B.size(); i++) {
        const auto& c = B.cells[i];
        if (i == 0) continue;
        std::vector<std::size_t> facets;
        for (auto f : c.facets) facets.push_back(at.at(f));
        at[i] = A.add(c.dim, c.type, facets, c.signs);
      }
      A.add(3, "", A.of_dimension(2));
      refused("a 3-cell bounded by two spheres pinched at a vertex", A, "pinched");
    }
    {
      auto C = mt::cubical_grid(1, 1, 1).complex;
      auto faces = C.of_dimension(2);
      C.cells.pop_back();
      faces.pop_back();
      C.add(3, "", faces);
      refused("a 3-cell bounded by an open box", C, "on 1 of its 2-cells");
    }
    {
      mt::RegularComplex X;
      const auto a = X.add_vertex(), b = X.add_vertex();
      X.add(2, "", {a, b});
      refused("a 2-cell with vertices as facets", X, "dimension 0");
    }
  }

  // ---- the face poset is Eulerian ----
  {
    const auto X = mt::cubical_grid(1, 1, 1).complex;
    const std::size_t n = X.size();
    std::vector<std::vector<char>> leq(n, std::vector<char>(n, 0));
    for (std::size_t t = 0; t < n; t++)
      for (std::size_t q : X.closure(t)) leq[q][t] = 1;
    const auto mu = morphe::combinatorics::mobius<long long>(leq);
    bool eulerian = true;
    std::size_t intervals = 0;
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = 0; b < n; b++)
        if (leq[a][b]) {
          intervals++;
          const long long want = (X.cells[b].dim - X.cells[a].dim) % 2 ? -1 : 1;
          eulerian = eulerian && mu[a][b] == want;
        }
    s.require("μ(σ, τ) = (−1)^{dim τ − dim σ} on every interval of the cube (" + std::to_string(intervals) + " intervals)", eulerian);
  }

  return s.done();
}
