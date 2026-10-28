// A check that every header of morphe::formal compiles together in one translation unit (their
// names share the namespace morphe::formal, and a collision would break any program that uses two
// of them, as the engine of the Genome Atlas does), and that a few of their names coexist: a
// monoid morphism, a semiring morphism, an edit path and a polygon.

#include <morphe/check.hpp>
#include <morphe/formal/aho_corasick.hpp>
#include <morphe/formal/alignment.hpp>
#include <morphe/formal/breaking.hpp>
#include <morphe/formal/bwt.hpp>
#include <morphe/formal/earley.hpp>
#include <morphe/formal/equivariant.hpp>
#include <morphe/formal/fm_index.hpp>
#include <morphe/formal/grammar.hpp>
#include <morphe/formal/inside_outside.hpp>
#include <morphe/formal/pbwt.hpp>
#include <morphe/formal/polytope.hpp>
#include <morphe/formal/r_index.hpp>
#include <morphe/formal/semiring.hpp>
#include <morphe/formal/succinct.hpp>
#include <morphe/formal/suffix_array.hpp>
#include <morphe/formal/words.hpp>

#include <cstdint>
#include <vector>

namespace mf = morphe::formal;

int main() {
  morphe::Suite s;
  static_assert(mf::SemiringMorphism<mf::ForgetCount>, "ForgetCount is a semiring morphism");
  static_assert(mf::Semiring<mf::PolytopeAlgebra>, "the polytope algebra is a semiring");
  const std::vector<std::uint8_t> x = {0, 1, 1}, y = {0, 1};
  const auto one = [](auto...) { return morphe::arithmetic::Integer(1); };
  s.equal("formal: every header compiles in one translation unit, and the alignments of words of lengths 3 and 2 number D(3, 2) = 25",
          mf::path_sum(x, y, mf::linear<mf::Counting>(one, one, one), mf::Mode::global), mf::delannoy(3, 2));
  return s.done();
}
