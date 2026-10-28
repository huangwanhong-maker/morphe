// The group algebra K[G] of a finite group over an exact field, its central idempotents and its
// Wedderburn decomposition.
//
// K[G] has the elements of G as a basis; an element is its vector of coefficients (Σ_g a_g g
// held as a_g at index g), added coordinate by coordinate and multiplied by convolution,
//
//   (Σ_g a_g g)(Σ_h b_h h) = Σ_{g,h} a_g b_h (g h),   O(|G|²) for dense factors.
//
// Over ℂ the algebra is semisimple and
//
//   ℂ[G] ≅ ⊕_χ M_{χ(1)}(ℂ),   the sum over the irreducible characters χ,
//
// (Maschke and Wedderburn; Serre 1977, Linear Representations of Finite Groups, ch. 6), the
// block of χ being the two-sided ideal ℂ[G] e_χ cut out by the central primitive idempotent
//
//   e_χ = (χ(1) / |G|) Σ_g χ(g⁻¹) g,
//
// which acts on any representation as the projection onto its χ-isotypic component (Serre 1977,
// §2.6, Theorem 8). They satisfy e_χ² = e_χ, e_χ e_ψ = 0 for χ ≠ ψ, Σ_χ e_χ = 1, and each is
// central; the dimension of ℂ[G] e_χ is χ(1)². Here all of this is computed, not assumed: the
// idempotents from the character table (characters.hpp), their relations by exact
// multiplication, and each dimension as the rank of x ↦ x e_χ by exact elimination
// (field_matrix.hpp) over ℚ(ζ_e), which holds every coefficient of e_χ; a rank does not change
// when the field is extended, so it is the dimension over ℂ. Any f ∈ K[G] is the sum of its
// isotypic components f e_χ.
//
// Over ℚ the algebra splits less. The central primitive idempotents of ℚ[G] are the sums of the
// e_χ over the orbits of the Galois group of ℚ(ζ_e) on the characters. Proof: such a sum is
// fixed by every σ_k, so its coefficients are rational; it is a central idempotent; and a
// central idempotent of ℚ[G] is a sum of e_χ (it is one of ℂ[G]) over a set of characters that
// σ_k must preserve, so a union of orbits. The block of the orbit of χ has dimension
// |orbit| χ(1)² over ℚ, since tensoring with ℂ gives the sum of the blocks of the orbit. For a
// cyclic group C_n the orbits are the characters of each order d | n, and
//
//   ℚ[C_n] ≅ ⊕_{d | n} ℚ(ζ_d),
//
// the block of the characters of order d being ℚ[u], u = g e_d for a generator g, with e_d its
// unit. Proof: g e_χ = χ(g) e_χ for a linear χ, and χ(g) is a primitive d-th root of unity, so
// Φ_d(u) = Φ_d(g) e_d = 0; Φ_d is irreducible, so ℚ[u] ≅ ℚ[x]/(Φ_d) = ℚ(ζ_d), of dimension φ(d),
// which is the dimension of the block. This is the decomposition exact rational arithmetic sees;
// over ℂ each ℚ(ζ_d) splits further into φ(d) copies of ℂ.
//
// Costs: a product O(|G|²) field operations; an idempotent O(|G|); the dimension of a block
// O(|G|² χ(1)²) by elimination (O(|G|³) for all blocks together). Exact throughout.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/arithmetic/cyclotomic.hpp>
#include <morphe/arithmetic/rational.hpp>

#include "characters.hpp"
#include "field_matrix.hpp"
#include "finite_group.hpp"

namespace morphe::algebra {

template <class K>
class GroupAlgebra {
 public:
  using Element = std::vector<K>;

  explicit GroupAlgebra(FiniteGroup G) : G_(std::move(G)) {}
  const FiniteGroup& group() const { return G_; }
  std::size_t dimension() const { return G_.order(); }

  Element zero() const { return Element(G_.order(), K(0)); }
  Element one() const { return basis(0); }
  // The basis element g.
  Element basis(std::uint32_t g) const {
    Element x = zero();
    x.at(g) = K(1);
    return x;
  }
  Element add(const Element& x, const Element& y) const {
    check(x), check(y);
    Element z = x;
    for (std::size_t g = 0; g < z.size(); g++) z[g] += y[g];
    return z;
  }
  Element subtract(const Element& x, const Element& y) const {
    check(x), check(y);
    Element z = x;
    for (std::size_t g = 0; g < z.size(); g++) z[g] -= y[g];
    return z;
  }
  Element scale(const K& c, const Element& x) const {
    check(x);
    Element z = x;
    for (auto& v : z) v *= c;
    return z;
  }
  // The convolution product x y.
  Element multiply(const Element& x, const Element& y) const {
    check(x), check(y);
    const std::uint32_t n = static_cast<std::uint32_t>(G_.order());
    Element z = zero();
    for (std::uint32_t a = 0; a < n; a++) {
      if (x[a].is_zero()) continue;
      for (std::uint32_t b = 0; b < n; b++)
        if (!y[b].is_zero()) z[G_.mul(a, b)] += x[a] * y[b];
    }
    return z;
  }
  // x is central exactly when its coefficients are constant on conjugacy classes. O(|G|²).
  bool is_central(const Element& x) const {
    check(x);
    const std::uint32_t n = static_cast<std::uint32_t>(G_.order());
    for (std::uint32_t g = 0; g < n; g++)
      for (std::uint32_t h = 0; h < n; h++)
        if (!(x[h] == x[G_.mul(G_.mul(g, h), G_.inv(g))])) return false;
    return true;
  }
  // The matrix of x ↦ x e in the basis of the elements: column g holds g e.
  FieldMatrix<K> right_multiplication(const Element& e) const {
    check(e);
    const std::uint32_t n = static_cast<std::uint32_t>(G_.order());
    FieldMatrix<K> R(n, n);
    for (std::uint32_t g = 0; g < n; g++)
      for (std::uint32_t h = 0; h < n; h++)
        if (!e[h].is_zero()) R(G_.mul(g, h), g) = e[h];
    return R;
  }
  // dim K[G] e, the left ideal generated by e (two-sided when e is central).
  std::size_t ideal_dimension(const Element& e) const { return right_multiplication(e).rank(); }

 private:
  void check(const Element& x) const {
    if (x.size() != G_.order()) throw std::invalid_argument("morphe::algebra::GroupAlgebra: an element of the wrong length");
  }
  FiniteGroup G_;
};

// The central primitive idempotent e_χ = (χ(1)/|G|) Σ_g χ(g⁻¹) g of ℚ(ζ_e)[G].
inline std::vector<Cyclotomic> central_idempotent(const CharacterTable& T, std::size_t chi) {
  const Cyclotomic scale(Rational(arithmetic::Integer(static_cast<long long>(T.degree(chi))), arithmetic::Integer(static_cast<long long>(T.order))));
  std::vector<Cyclotomic> by_class;
  for (std::size_t j = 0; j < T.classes.size(); j++) by_class.push_back(scale * T.values[chi][T.inverse_class[j]]);
  std::vector<Cyclotomic> e(T.order);
  for (std::size_t g = 0; g < T.order; g++) e[g] = by_class[T.class_of[g]];
  return e;
}

// A block of ℂ[G] ≅ ⊕_χ M_{χ(1)}(ℂ): the ideal ℂ[G] e_χ, isomorphic to M_d(ℂ), d = χ(1).
struct WedderburnBlock {
  std::size_t character = 0;
  std::size_t degree = 0;     // d = χ(1)
  std::size_t dimension = 0;  // the rank of x ↦ x e_χ, measured; d² by the theorem
  std::vector<Cyclotomic> idempotent;
};

inline std::vector<WedderburnBlock> wedderburn(const GroupAlgebra<Cyclotomic>& A, const CharacterTable& T) {
  if (A.dimension() != T.order) throw std::invalid_argument("morphe::algebra::wedderburn: the table is of another group");
  std::vector<WedderburnBlock> blocks;
  for (std::size_t chi = 0; chi < T.size(); chi++) {
    WedderburnBlock B;
    B.character = chi;
    B.degree = T.degree(chi);
    B.idempotent = central_idempotent(T, chi);
    B.dimension = A.ideal_dimension(B.idempotent);
    blocks.push_back(std::move(B));
  }
  return blocks;
}

// The isotypic components f e_χ of f, one for each irreducible character; they sum to f.
inline std::vector<std::vector<Cyclotomic>> isotypic_components(const GroupAlgebra<Cyclotomic>& A, const CharacterTable& T,
                                                                const std::vector<Cyclotomic>& f) {
  std::vector<std::vector<Cyclotomic>> parts;
  for (std::size_t chi = 0; chi < T.size(); chi++) parts.push_back(A.multiply(f, central_idempotent(T, chi)));
  return parts;
}

// An element of ℚ[G] read in ℚ(ζ_e)[G].
inline std::vector<Cyclotomic> to_cyclotomic(const std::vector<Rational>& f) { return {f.begin(), f.end()}; }

// A block of ℚ[G]: the Galois orbit of characters it gathers, its central primitive idempotent
// Σ_{χ in the orbit} e_χ (rational), and its dimension over ℚ, |orbit| χ(1)² by the theorem.
struct RationalBlock {
  std::vector<std::size_t> characters;
  std::vector<Rational> idempotent;
  std::size_t dimension = 0;
};

inline std::vector<RationalBlock> rational_wedderburn(const GroupAlgebra<Rational>& A, const CharacterTable& T) {
  if (A.dimension() != T.order) throw std::invalid_argument("morphe::algebra::rational_wedderburn: the table is of another group");
  const long long e = static_cast<long long>(T.exponent);
  std::vector<char> seen(T.size(), 0);
  std::vector<RationalBlock> blocks;
  for (std::size_t chi = 0; chi < T.size(); chi++) {
    if (seen[chi]) continue;
    RationalBlock B;
    for (long long k = 1; k <= e; k++) {
      if (std::gcd(k, e) != 1) continue;
      const std::size_t psi = galois_conjugate(T, chi, k);
      if (!seen[psi]) seen[psi] = 1, B.characters.push_back(psi);
    }
    std::sort(B.characters.begin(), B.characters.end());
    std::vector<Cyclotomic> sum(T.order, Cyclotomic(0));
    for (std::size_t psi : B.characters) {
      const auto e_psi = central_idempotent(T, psi);
      for (std::size_t g = 0; g < T.order; g++) sum[g] += e_psi[g];
    }
    for (const auto& x : sum) B.idempotent.push_back(x.to_rational());
    B.dimension = A.ideal_dimension(B.idempotent);
    blocks.push_back(std::move(B));
  }
  return blocks;
}

// ℚ[C_n] ≅ ⊕_{d | n} ℚ(ζ_d): for each divisor d, the block of the characters of order d.
struct CyclicBlock {
  std::size_t divisor = 1;  // d: the block is isomorphic to ℚ(ζ_d)
  RationalBlock block;
};

inline std::vector<CyclicBlock> cyclic_rational_decomposition(std::size_t n) {
  const GroupAlgebra<Rational> A(cyclic(n));
  const CharacterTable T = character_table(A.group());
  std::vector<CyclicBlock> out;
  for (auto& B : rational_wedderburn(A, T)) {
    // The order of χ(g) as a root of unity, g = 1 the generator (the identity when n = 1).
    const Cyclotomic x = T.value(B.characters[0], n > 1 ? 1u : 0u);
    std::size_t d = 1;
    for (Cyclotomic y = x; !(y == Cyclotomic(1)); y *= x) d++;
    out.push_back({d, std::move(B)});
  }
  std::sort(out.begin(), out.end(), [](const CyclicBlock& a, const CyclicBlock& b) { return a.divisor < b.divisor; });
  return out;
}

}  // namespace morphe::algebra
