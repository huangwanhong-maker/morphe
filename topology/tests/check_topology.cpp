// Checks of morphe/topology against the classical computations: the homology and fundamental
// groups of surfaces given by words, of projective and lens spaces and spheres, of products
// (Künneth) and lattices; the classification of surfaces, orientability and boundary; the
// fundamental class of an orientable surface against its orientation; triangulations and
// their links, subdivisions, cones and suspensions; manifold recognition; collapses (the
// dunce hat contractible and not collapsible); dual 2-complexes; and coloured graphs: their
// pseudomanifolds, jackets, Gurau's face formula and the Fuss–Catalan count of melonic graphs.

#include <morphe/check.hpp>
#include <morphe/topology/collapse.hpp>
#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/dual.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/manifolds.hpp>
#include <morphe/topology/surfaces.hpp>
#include <morphe/algebra/finite_group.hpp>

#include <string>
#include <vector>

namespace mt = morphe::topology;
namespace ma = morphe::algebra;
using morphe::arithmetic::Integer;

namespace {

// The homology as text, "ℤ ; ℤ^2 ; ℤ".
std::string H(const std::vector<ma::HomologyGroup>& h) {
  std::string s;
  for (std::size_t k = 0; k < h.size(); k++) s += (k ? " ; " : "") + h[k].str();
  return s;
}
std::string H(const mt::CWComplex& X) { return H(mt::homology(X)); }
std::string H(const mt::SimplicialComplex& K) { return H(ma::homology(K.chain_complex())); }

std::size_t factorial(std::size_t n) { return n <= 1 ? 1 : n * factorial(n - 1); }

}  // namespace

int main() {
  morphe::Suite s;

  // ---- surfaces from words -----------------------------------------------------------------------------------------
  {
    struct Case { const char* word; const char* homology; const char* name; bool orientable; std::size_t genus, crosscaps, boundary; };
    const Case cases[] = {
        {"aA", "ℤ ; 0 ; ℤ", "sphere", true, 0, 0, 0},
        {"abAB", "ℤ ; ℤ^2 ; ℤ", "torus", true, 1, 0, 0},
        {"a1b1A1B1a2b2A2B2", "ℤ ; ℤ^4 ; ℤ", "orientable surface of genus 2", true, 2, 0, 0},
        {"aa", "ℤ ; ℤ/2 ; 0", "projective plane", false, 0, 1, 0},
        {"abAb", "ℤ ; ℤ ⊕ ℤ/2 ; 0", "Klein bottle", false, 0, 2, 0},
        {"c1c1c2c2c3c3", "ℤ ; ℤ^2 ⊕ ℤ/2 ; 0", "non-orientable surface with 3 crosscaps", false, 0, 3, 0},
        {"abac", "ℤ ; ℤ ; 0", "Möbius band", false, 0, 1, 1},
        {"abAc", "ℤ ; ℤ ; 0", "annulus", true, 0, 0, 2},
        {"abcd", "ℤ ; 0 ; 0", "disc", true, 0, 0, 1},
        {"abAB acAC", "ℤ ; ℤ^3 ; ℤ^2", "", false, 0, 0, 0},  // two tori sharing the circle a: not a surface
        {"abAB cdCD", "ℤ^2 ; ℤ^4 ; ℤ^2", "", false, 0, 0, 0},  // two disjoint tori: not connected
    };
    for (const auto& c : cases) {
      const auto X = mt::polygon_complex(c.word);
      const auto S = mt::classify_surface(X);
      const bool surface_expected = std::string(c.name) != "";
      bool ok = H(X) == c.homology;
      if (surface_expected) ok = ok && S.surface && S.name == c.name && S.orientable == c.orientable && S.genus == c.genus && S.crosscaps == c.crosscaps && S.boundary_components == c.boundary;
      else ok = ok && !S.surface;
      s.require(std::string("topology: the complex of the word ") + c.word + " has homology " + c.homology + (surface_expected ? " and is the " + std::string(c.name) : " and is not a surface"),
                ok, H(X) + " / " + (S.surface ? S.name : S.reason));
    }
    // The fundamental group of the torus and its loops.
    const auto T = mt::polygon_complex("abAB");
    const auto G = mt::fundamental_group(T);
    bool loops = G.loops.size() == 2;
    for (const auto& l : G.loops) loops = loops && T.walk(0, l) == 0;
    s.require("topology: π1(torus) = ⟨a, b | a b a⁻¹ b⁻¹⟩ with closed loops at the base", G.presentation.generators == 2 && G.presentation.relators.size() == 1 && loops,
              G.presentation.str());
    // The fundamental class of a closed orientable surface is Σ o_f f for its orientation o.
    const auto X = mt::torus_lattice({5, 4});
    const auto SX = mt::classify_surface(X);
    const auto HX = mt::homology(X);
    bool fundamental = SX.surface && SX.orientable && HX[2].free_generators.size() == 1 && HX[2].free_generators[0].size() == X.faces.size();
    if (fundamental) {
      const auto& z = HX[2].free_generators[0];
      const int sign = (z[0].value.sign() > 0 ? 1 : -1) * SX.orientation[0];
      for (const auto& e : z) fundamental = fundamental && e.value == Integer(sign * SX.orientation[e.index]);
    }
    s.require("topology: the generator of H2 of the 5 × 4 lattice torus is its fundamental class, ± Σ o_f f", fundamental);
  }
  // ---- CW spaces --------------------------------------------------------------------------------------------------------
  s.equal("topology: H(ℝP⁵) = ℤ, ℤ/2, 0, ℤ/2, 0, ℤ", H(mt::real_projective_space(5)), std::string("ℤ ; ℤ/2 ; 0 ; ℤ/2 ; 0 ; ℤ"));
  s.equal("topology: H(ℂP³) = ℤ, 0, ℤ, 0, ℤ, 0, ℤ", H(mt::complex_projective_space(3)), std::string("ℤ ; 0 ; ℤ ; 0 ; ℤ ; 0 ; ℤ"));
  s.equal("topology: H(S⁴) = ℤ, 0, 0, 0, ℤ", H(mt::sphere(4)), std::string("ℤ ; 0 ; 0 ; 0 ; ℤ"));
  {
    const auto L = mt::lens_space(5);
    const auto G = mt::fundamental_group(L);
    s.require("topology: the lens space L(5, q) has H = ℤ, ℤ/5, 0, ℤ and π1 = ⟨a | a⁵⟩", H(L) == "ℤ ; ℤ/5 ; 0 ; ℤ" && G.simplified.presentation.str() == "⟨a | a a a a a⟩",
              H(L) + " / " + G.simplified.presentation.str());
  }
  // ---- products and lattices -----------------------------------------------------------------------------------------
  {
    const auto P = mt::product(mt::real_projective_space(2), mt::circle(3));
    P.validate();
    s.equal("topology: H(ℝP² × S¹) = ℤ, ℤ ⊕ ℤ/2, ℤ/2, 0 (Künneth)", H(P), std::string("ℤ ; ℤ ⊕ ℤ/2 ; ℤ/2 ; 0"));
    const auto T3 = mt::torus_lattice({3, 4, 3});
    T3.validate();
    s.equal("topology: H of the 3 × 4 × 3 lattice 3-torus is ℤ, ℤ³, ℤ³, ℤ", H(T3), std::string("ℤ ; ℤ^3 ; ℤ^3 ; ℤ"));
    const auto T4 = mt::torus_lattice({3, 3, 3, 3});
    s.equal("topology: Betti numbers of the 3⁴ lattice 4-torus are 1, 4, 6, 4, 1", mt::betti(T4), std::vector<std::size_t>{1, 4, 6, 4, 1});
    const auto big = mt::torus_lattice({32, 32});
    const auto Hb = mt::homology(big);
    bool cycles = Hb[1].free_generators.size() == 2;
    const auto C = big.chain_complex();
    for (const auto& z : Hb[1].free_generators) cycles = cycles && C.boundary[1].apply(z).empty();
    s.require("topology: the 32 × 32 lattice torus: H1 = ℤ² with generators that are cycles", H(Hb) == "ℤ ; ℤ^2 ; ℤ" && cycles, H(Hb));
    const auto Q = mt::product(mt::interval(3), mt::interval(4));
    s.require("topology: the 3 × 4 grid square is contractible and collapsible", H(Q) == "ℤ ; 0 ; 0" && mt::collapse(Q).collapsible);
  }
  // ---- presentation complexes and homomorphisms -------------------------------------------------------------------------
  {
    ma::Presentation P;
    P.generators = 2;
    P.names = {"a", "b"};
    P.relators = {{1, 1}, {2, 2, 2}, {1, 2, 1, 2, 1, 2, 1, 2, 1, 2}};  // ⟨a, b | a², b³, (a b)⁵⟩ ≅ A5
    const auto X = mt::presentation_complex(P);
    const auto G = mt::fundamental_group(X);
    const auto homs = ma::count_homomorphisms(G.presentation, ma::alternating(5));
    s.require("topology: the presentation complex of ⟨a, b | a², b³, (ab)⁵⟩ has H1 = 0, H2 = ℤ, and its π1 has 1 + |Aut A5| = 121 maps to A5",
              H(X) == "ℤ ; 0 ; ℤ" && homs == 121, H(X) + " / " + std::to_string(homs));
  }
  // ---- triangulations ----------------------------------------------------------------------------------------------------
  {
    const auto T7 = mt::torus_7();
    const auto S7 = mt::classify_surface(T7.cw());
    const auto L = mt::renumbered(T7.link({0}));
    s.require("topology: the 7-vertex torus: χ = 0, H = ℤ, ℤ², ℤ, a torus, the link of a vertex a hexagon",
              T7.euler_characteristic() == 0 && H(T7) == "ℤ ; ℤ^2 ; ℤ" && S7.name == "torus" && L.count(0) == 6 && L.count(1) == 6);
    const auto P6 = mt::projective_plane_6();
    s.require("topology: the 6-vertex projective plane: H = ℤ, ℤ/2, 0, non-orientable", H(P6) == "ℤ ; ℤ/2 ; 0" && mt::classify_surface(P6.cw()).name == "projective plane");
    const auto B = mt::barycentric_subdivision(T7);
    s.require("topology: the barycentric subdivision of the 7-vertex torus has 42, 126, 84 simplices and the same homology",
              B.counts() == std::vector<std::size_t>{42, 126, 84} && H(B) == "ℤ ; ℤ^2 ; ℤ");
    const auto S3 = mt::simplex_boundary(4), S3b = mt::suspension(mt::simplex_boundary(3));
    const auto M1 = mt::manifold_type(S3), M2 = mt::manifold_type(S3b);
    s.require("topology: ∂Δ⁴ and the suspension of ∂Δ³ are closed combinatorial 3-manifolds with the homology of S³",
              H(S3) == "ℤ ; 0 ; 0 ; ℤ" && H(S3b) == "ℤ ; 0 ; 0 ; ℤ" && M1.manifold && M1.closed && M1.exact && M2.manifold && M2.closed);
    const auto C = mt::cone(mt::octahedron());
    s.require("topology: the cone on the octahedron is a contractible 3-ball", H(C) == "ℤ ; 0 ; 0 ; 0" && mt::manifold_type(C).manifold && !mt::manifold_type(C).closed);
    const auto pinch = mt::SimplicialComplex::from_facets({{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 2, 3}, {0, 4, 5}, {0, 5, 6}, {0, 6, 4}, {4, 5, 6}});
    const auto M3 = mt::manifold_type(pinch);
    s.require("topology: two spheres sharing a vertex are a pseudomanifold but not a manifold", M3.pseudomanifold && !M3.manifold, M3.reason);
    s.require("topology: the cone on a hexagon (a triangulated disc) collapses to a point", mt::collapse(mt::cone(mt::renumbered(L)).cw()).collapsible);
  }
  // ---- collapses -----------------------------------------------------------------------------------------------------------
  {
    const auto D = mt::dunce_hat();
    const auto G = mt::fundamental_group(D);
    const auto C = mt::collapse(D);
    s.require("topology: the dunce hat is contractible (H = ℤ, 0, 0, π1 trivial) and not collapsible",
              H(D) == "ℤ ; 0 ; 0" && G.simplified.presentation.generators == 0 && C.attempted && !C.collapsible && C.steps == 0);
    s.require("topology: the torus does not collapse", !mt::collapse(mt::polygon_complex("abAB")).collapsible);
  }
  // ---- dual 2-complexes ------------------------------------------------------------------------------------------------------
  {
    const auto D = mt::dual_two_complex(mt::simplex_boundary(4));
    const auto G = mt::fundamental_group(D.complex);
    s.require("topology: the dual 2-complex of ∂Δ⁴ has 5 vertices, 10 edges, 10 faces, π1 = 1, H2 = ℤ⁴",
              D.complex.counts() == std::vector<std::size_t>{5, 10, 10} && G.simplified.presentation.generators == 0 && H(D.complex) == "ℤ ; 0 ; ℤ^4",
              H(D.complex));
    const auto T = mt::dual_two_complex(mt::torus_7());
    s.require("topology: the dual of the 7-vertex torus is the hexagonal torus (14, 21, 7)", T.complex.counts() == std::vector<std::size_t>{14, 21, 7} &&
                                                                                          mt::classify_surface(T.complex).name == "torus");
  }
  // ---- coloured graphs -----------------------------------------------------------------------------------------------------
  {
    const auto M = mt::elementary_melon(3);
    const auto X = M.dual_complex();
    X.validate();
    const auto two = M.two_complex();
    s.require("topology: the elementary 4-coloured melon encodes S³; its 2-complex (2, 4, 6) is simply connected",
              H(X) == "ℤ ; 0 ; 0 ; ℤ" && two.counts() == std::vector<std::size_t>{2, 4, 6} && mt::fundamental_group(two).simplified.presentation.generators == 0 && M.degree() == 0);
    auto G = M;
    for (std::size_t k = 0; k < 6; k++) G = mt::insert_melon(G, k % 4, k % G.p);
    s.require("topology: melonic insertions keep degree 0 and the 3-sphere", G.degree() == 0 && H(G.dual_complex()) == "ℤ ; 0 ; 0 ; ℤ" && G.p == 7);
    int formula = 0, surfaces = 0, valid = 0;
    for (std::uint64_t seed = 1; seed <= 40; seed++) {
      const auto R = mt::random_coloured_graph(seed % 2 ? 3 : 4, 2 + seed % 6, seed);
      if (!R.connected()) continue;
      const std::size_t D = R.dimension();
      // F = D(D-1)p/2 + D - 2ω/(D-1)!
      const long long lhs = static_cast<long long>(R.face_count()) * static_cast<long long>(factorial(D - 1));
      const long long rhs = (static_cast<long long>(D * (D - 1) * R.p / 2 + D)) * static_cast<long long>(factorial(D - 1)) - 2 * static_cast<long long>(R.degree());
      if (lhs != rhs) formula++;
      try { R.dual_complex().validate(); } catch (...) { valid++; }
      const auto S2 = mt::random_coloured_graph(2, 2 + seed % 7, seed);
      if (S2.connected()) {
        const auto surf = mt::classify_surface(S2.dual_complex());
        if (!surf.surface || !surf.orientable || surf.genus != S2.degree()) surfaces++;
      }
    }
    s.equal("topology: Gurau's formula F = D(D-1)p/2 + D - 2ω/(D-1)! on random 4- and 5-coloured graphs (failures)", formula, 0);
    s.equal("topology: the pseudomanifolds of random coloured graphs are chain complexes, ∂∂ = 0 (failures)", valid, 0);
    s.equal("topology: a 3-coloured graph encodes an orientable surface whose genus is its degree (failures)", surfaces, 0);
    // Rooted melonic 4-coloured graphs with 2p vertices: the Fuss–Catalan numbers 1, 4, 22, 140.
    std::vector<std::size_t> melonic;
    for (std::size_t p = 1; p <= 4; p++) {
      std::vector<std::size_t> id(p);
      std::iota(id.begin(), id.end(), 0);
      std::vector<std::size_t> s1 = id, s2 = id, s3 = id;
      std::size_t count = 0;
      do {
        do {
          do {
            const mt::ColouredGraph Gr{4, p, {id, s1, s2, s3}};
            if (Gr.connected() && Gr.degree() == 0) count++;
          } while (std::next_permutation(s3.begin(), s3.end()));
        } while (std::next_permutation(s2.begin(), s2.end()));
      } while (std::next_permutation(s1.begin(), s1.end()));
      melonic.push_back(count / factorial(p - 1));
    }
    s.equal("topology: rooted melonic 4-coloured graphs with 2, 4, 6, 8 vertices are counted by Fuss–Catalan numbers", melonic, std::vector<std::size_t>{1, 4, 22, 140});
  }
  return s.done();
}
