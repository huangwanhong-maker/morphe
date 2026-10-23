// Checks of morphe/topology/bubbles.hpp: the automorphism groups of the classical bubbles
// ((k!)² for (φ̄φ)^k, k! for (φ̄·φ)^k, k rotations for tr (MM†)^k, two for a quartic melon,
// three for the wheel and the necklace, one for φ̄·φ); bubbles glued into Feynman graphs (their
// colours, connectedness, faces, and the faces of colours (0, i) of a single trace counted as
// closed index lines); bubbles refused when they are not bubbles.

#include <morphe/check.hpp>
#include <morphe/topology/bubbles.hpp>

#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace mt = morphe::topology;

int main() {
  morphe::Suite s;
  auto aut = [](const mt::ColouredGraph& B) { return mt::bubble_automorphisms(B).size(); };
  s.equal("bubbles: |Aut (φ̄φ)³| = (3!)² = 36", aut(mt::scalar_bubble(3)), std::size_t{36});
  s.equal("bubbles: |Aut (φ̄·φ)⁴| = 4! = 24", aut(mt::vector_bubble(4)), std::size_t{24});
  s.equal("bubbles: |Aut tr (MM†)⁵| = 5, its rotations", aut(mt::trace_bubble(5)), std::size_t{5});
  s.equal("bubbles: |Aut| of the quartic melon of colour 2 in rank 4 is 2", aut(mt::quartic_melonic_bubble(4, 1)), std::size_t{2});
  s.equal("bubbles: |Aut K₃,₃| = 3", aut(mt::wheel_bubble()), std::size_t{3});
  s.equal("bubbles: |Aut| of the necklace is 3", aut(mt::necklace_bubble()), std::size_t{3});
  s.equal("bubbles: φ̄·φ in rank 3 has no automorphism but the identity", aut(mt::quadratic_bubble(3)), std::size_t{1});
  {
    // Two quartic melons glued to themselves (σ₀ = id): two components, rank 3 + propagators.
    const auto B = mt::quartic_melonic_bubble(3, 0);
    const auto G = mt::glue({B, B}, {0, 1, 2, 3});
    std::size_t comps = 0;
    G.components(0xF, &comps);
    s.require("bubbles: two melons glued each to itself form 4-coloured graphs in 2 components", G.colours == 4 && G.p == 4 && comps == 2);
    const auto H = mt::glue({B, B}, {2, 3, 0, 1});
    s.require("bubbles: crossed, they form one connected graph", H.connected());
  }
  {
    // ⟨tr (MM†)^k⟩ with σ₀ = id: the index lines of colours (0, 1) and (0, 2), N^{k + 1}… the
    // planar term: one face of colour 1 runs around, k of colour 2 close at once.
    const std::size_t k = 4;
    const auto G = mt::glue({mt::trace_bubble(k)}, {0, 1, 2, 3});
    s.require("bubbles: tr (MM†)⁴ with each M joined to its own M†: k faces of colours (0, 2) and one of (0, 1)",
              G.faces(0, 2).size() == k && G.faces(0, 1).size() == 1, std::to_string(G.faces(0, 1).size()) + ", " + std::to_string(G.faces(0, 2).size()));
    const auto R0 = mt::glue({mt::scalar_bubble(2), mt::scalar_bubble(1)}, {2, 0, 1});
    s.require("bubbles: rank-0 bubbles glue into a graph of one colour, the propagators", R0.colours == 1 && R0.p == 3);
  }
  bool refused = false;
  try {
    mt::glue({mt::trace_bubble(2), mt::quartic_melonic_bubble(3, 0)}, {0, 1, 2, 3});
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  s.require("bubbles: bubbles of different ranks are not glued", refused);
  refused = false;
  try {
    mt::validate_bubble(mt::ColouredGraph{2, 2, {{0, 0}, {0, 1}}});
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  s.require("bubbles: a colour that is not a matching is refused", refused);
  return s.done();
}
