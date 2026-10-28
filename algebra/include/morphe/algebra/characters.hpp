// The complex character table of a finite group, exactly, by the Dixon–Schneider method.
//
// The irreducible complex characters of a finite group G take their values in ℤ[ζ_e], e the
// exponent of G, and there are as many as there are conjugacy classes. This header computes
// them all from G's multiplication table, as cyclotomic numbers (arithmetic/cyclotomic.hpp),
// with no value typed or approximated.
//
// The method (Dixon 1967, "High speed computation of group characters", Numer. Math. 10,
// 446–450; Schneider 1990, "Dixon's character table algorithm revisited", J. Symbolic Comput.
// 9, 601–606). Let K_1 = {1}, K_2, …, K_r be the classes, h_j = |K_j|, g_j ∈ K_j, and C_j the
// class sums, which span the centre of ℂ[G]. With
//
//   a_{jkl} = #{(x, y) ∈ K_j × K_k : x y = g_l},   so that C_j C_k = Σ_l a_{jkl} C_l,
//
// the central character ω_χ(C_j) = h_j χ(g_j) / χ(1) of each irreducible χ is an algebra
// homomorphism of the centre, so ω_χ(C_j) ω_χ(C_k) = Σ_l a_{jkl} ω_χ(C_l): the vector
// w_χ = (ω_χ(C_1), …, ω_χ(C_r)), whose first coordinate is 1, is a common eigenvector of the
// class matrices M_j = (a_{jkl})_{k,l}, with eigenvalue ω_χ(C_j) for M_j. The r vectors w_χ are
// linearly independent, so they are, up to scalars, the only common eigenvectors.
//
// Dixon computes them modulo a prime p ≡ 1 (mod e). Then p does not divide |G| (every prime
// dividing |G| divides e), 𝔽_p contains a primitive e-th root of unity z, and ζ_e ↦ z defines a
// ring homomorphism θ: ℤ[ζ_e] → 𝔽_p, defined on every ω_χ(C_j), an algebraic integer of ℚ(ζ_e)
// (Isaacs 1976, ch. 3), whose ring of integers is ℤ[ζ_e] (Washington 1997, ch. 2). The images
// θ(w_χ) remain linearly independent: the reduced
// table satisfies the column orthogonality relations, whose right-hand sides |C_G(g_j)| are
// units modulo p, and θ(χ(1)) ≠ 0 since χ(1) < p. So the common eigenspaces of the M_j over 𝔽_p
// are the r lines θ(w_χ). They are found as Schneider does, one class matrix at a time: each
// space already found is split into the eigenspaces of the next M_j restricted to it (its
// characteristic polynomial by reduction to Hessenberg form, its roots by trying every element
// of 𝔽_p), until every space is a line. The line normalised to first coordinate 1 is θ(w_χ).
//
// From θ(w_χ): Σ_j h_j χ(g_j) χ(g_j⁻¹) = |G| gives Σ_j ω_χ(C_j) ω_χ(C_{j*}) / h_j = |G| / χ(1)²
// (K_{j*} the class of g_j⁻¹), so χ(1)² modulo p, and χ(1) is the unique root in [1, √|G|]
// when p > 2√|G| (two roots x, y ≤ √|G| with x² ≡ y² would have p | (x − y)(x + y), x + y < p).
// Then θ(χ(g_j)) = θ(ω_χ(C_j)) χ(1) / h_j. To lift a value to ℤ[ζ_e], let o be the order of g = g_j:
// the eigenvalues of ρ(g) are o-th roots of unity, ζ_o^k with multiplicity m_k, so
//
//   χ(g^i) = Σ_k m_k ζ_o^{ik},   m_k = (1/o) Σ_{i<o} χ(g^i) ζ_o^{−ik},
//
// and the m_k are integers in [0, χ(1)], below p: they are read exactly from their residues
// θ(m_k) = (1/o) Σ_i θ(χ(g^i)) z_o^{−ik}, z_o = z^{e/o}, and χ(g) = Σ_k m_k ζ_o^k exactly. The
// prime is the least p ≡ 1 (mod e) with p² > 4|G|. A different z would permute the characters by
// a Galois automorphism; the table, as a set, does not depend on it.
//
// Termination and cost: it terminates on every finite group. With n = |G| and r classes:
// the classes O(n²) (conjugacy_classes); the coefficients a_{jkl} O(n r); the splitting at most
// r − 1 class matrices, each O(r³) to restrict, and in all at most 2r spaces split, each O(d³)
// for its characteristic polynomial and for each of its eigenspaces plus O(p d) to find the
// roots, O(r⁴ + p r²) together; the lifting O(r Σ_j o(g_j)²). The prime exists by Dirichlet's
// theorem on primes in arithmetic progressions and is found by trying 1 + k e in turn; it is a
// small multiple of e for the groups met here (61 for S_5 and S_6). Everything is
// exact; a failure of an internal consistency check (a space that does not split, a multiplicity
// outside [0, χ(1)]) would be a defect and throws std::logic_error.
//
// The characters are ordered by degree, and among equal degrees by the multiplicity vectors
// (m_k) class by class, in decreasing lexicographic order, which puts the trivial character
// first. The classes are those of conjugacy_classes: the identity's first, each sorted, its
// representative its least element.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/arithmetic/cyclotomic.hpp>
#include <morphe/arithmetic/modular.hpp>
#include <morphe/arithmetic/rational.hpp>

#include "finite_group.hpp"

namespace morphe::algebra {

using arithmetic::Cyclotomic;
using arithmetic::Rational;

struct CharacterTable {
  std::size_t order = 1;     // |G|
  std::size_t exponent = 1;  // e: every value lies in ℚ(ζ_e)
  std::vector<std::vector<std::uint32_t>> classes;  // the conjugacy classes, the identity's first
  std::vector<std::size_t> class_sizes;
  std::vector<std::uint32_t> representatives;        // g_j, the least element of K_j
  std::vector<std::size_t> representative_orders;    // o(g_j)
  std::vector<std::size_t> class_of;                 // the class of each element
  std::vector<std::size_t> inverse_class;            // j*: the class of g_j⁻¹
  std::vector<std::vector<Cyclotomic>> values;       // values[χ][j] = χ(g_j), each in ℚ(ζ_{o(g_j)})
  std::uint64_t prime = 0;  // the prime p of the modular computation
  std::uint64_t root = 0;   // the primitive e-th root of unity modulo p that ζ_e was sent to

  // The number of irreducible characters, which is the number of classes.
  std::size_t size() const { return values.size(); }
  const Cyclotomic& operator()(std::size_t chi, std::size_t j) const { return values[chi][j]; }
  // χ(g) for an element g of the group.
  const Cyclotomic& value(std::size_t chi, std::uint32_t g) const { return values[chi][class_of[g]]; }
  std::size_t degree(std::size_t chi) const {
    return static_cast<std::size_t>(values[chi][0].to_rational().numerator().to_int64());
  }
};

// a_{jkl} = #{(x, y) ∈ K_j × K_k : x y = g_l} at index (j r + k) r + l, for the classes of T
// (only the classes, class_of and representatives are read). O(|G| r).
inline std::vector<std::uint64_t> class_multiplication_coefficients(const FiniteGroup& G, const CharacterTable& T) {
  const std::size_t r = T.classes.size();
  std::vector<std::uint64_t> a(r * r * r, 0);
  for (std::size_t l = 0; l < r; l++) {
    const std::uint32_t z = T.representatives[l];
    for (std::uint32_t x = 0; x < G.n; x++) {
      const std::uint32_t y = G.mul(G.inv(x), z);
      a[(T.class_of[x] * r + T.class_of[y]) * r + l]++;
    }
  }
  return a;
}

// The class matrix M_j = (a_{jkl})_{k,l}, row-major: a_{jkl} = #{x ∈ K_j : x⁻¹ g_l ∈ K_k}.
// O(h_j r) operations and O(r²) memory, so that the splitting never holds all r³ coefficients.
inline std::vector<std::uint64_t> class_matrix(const FiniteGroup& G, const CharacterTable& T, std::size_t j) {
  const std::size_t r = T.classes.size();
  std::vector<std::uint64_t> M(r * r, 0);
  for (std::size_t l = 0; l < r; l++)
    for (auto x : T.classes[j]) M[T.class_of[G.mul(G.inv(x), T.representatives[l])] * r + l]++;
  return M;
}

namespace detail {

using arithmetic::mod_add;
using arithmetic::mod_mul;
using arithmetic::mod_pow;
using arithmetic::mod_sub;
using arithmetic::Residue;
using FpRows = std::vector<std::vector<Residue>>;

inline Residue fp_inverse(Residue a, Residue p) { return mod_pow(a, p - 2, p); }

// Reduced row echelon form of the rows over 𝔽_p, in place, zero rows dropped; returns the
// pivot columns. O(rows · cols · rank).
inline std::vector<std::size_t> fp_rref(FpRows& rows, Residue p) {
  std::vector<std::size_t> pivots;
  const std::size_t cols = rows.empty() ? 0 : rows[0].size();
  std::size_t rank = 0;
  for (std::size_t c = 0; c < cols && rank < rows.size(); c++) {
    std::size_t i = rank;
    while (i < rows.size() && rows[i][c] == 0) i++;
    if (i == rows.size()) continue;
    std::swap(rows[i], rows[rank]);
    const Residue inv = fp_inverse(rows[rank][c], p);
    for (auto& x : rows[rank]) x = mod_mul(x, inv, p);
    for (std::size_t k = 0; k < rows.size(); k++) {
      if (k == rank || rows[k][c] == 0) continue;
      const Residue f = rows[k][c];
      for (std::size_t j = c; j < cols; j++) rows[k][j] = mod_sub(rows[k][j], mod_mul(f, rows[rank][j], p), p);
    }
    pivots.push_back(c);
    rank++;
  }
  rows.resize(rank);
  return pivots;
}

// A basis of the kernel {u : A u = 0} of a matrix given by its rows.
inline FpRows fp_kernel(FpRows A, std::size_t cols, Residue p) {
  const auto pivots = fp_rref(A, p);
  std::vector<char> is_pivot(cols, 0);
  for (auto c : pivots) is_pivot[c] = 1;
  FpRows kernel;
  for (std::size_t f = 0; f < cols; f++) {
    if (is_pivot[f]) continue;
    std::vector<Residue> u(cols, 0);
    u[f] = 1;
    for (std::size_t i = 0; i < pivots.size(); i++) u[pivots[i]] = (p - A[i][f]) % p;
    kernel.push_back(std::move(u));
  }
  return kernel;
}

// The characteristic polynomial det(x I − A) of a square matrix over 𝔽_p, low degree first.
// A is brought to upper Hessenberg form H by similarity (row k ← row k − u row m together with
// column m ← column m + u column k), and then, with p_0 = 1,
//   p_m = (x − h_{mm}) p_{m−1} − Σ_{i<m} h_{im} (h_{i+1,i} h_{i+2,i+1} ⋯ h_{m,m−1}) p_{i−1}
// (Cohen 1993, A Course in Computational Algebraic Number Theory, §2.2). O(n³).
inline std::vector<Residue> fp_charpoly(FpRows H, Residue p) {
  const std::size_t n = H.size();
  for (std::size_t m = 1; m + 1 < n; m++) {
    std::size_t i = m;
    while (i < n && H[i][m - 1] == 0) i++;
    if (i == n) continue;
    if (i != m) {
      std::swap(H[i], H[m]);
      for (auto& row : H) std::swap(row[i], row[m]);
    }
    const Residue inv = fp_inverse(H[m][m - 1], p);
    for (std::size_t k = m + 1; k < n; k++) {
      const Residue u = mod_mul(H[k][m - 1], inv, p);
      if (u == 0) continue;
      for (std::size_t c = 0; c < n; c++) H[k][c] = mod_sub(H[k][c], mod_mul(u, H[m][c], p), p);
      for (std::size_t c = 0; c < n; c++) H[c][m] = mod_add(H[c][m], mod_mul(u, H[c][k], p), p);
    }
  }
  // 1-based in the recurrence: h_{ab} = H[a−1][b−1].
  std::vector<std::vector<Residue>> P(n + 1);
  P[0] = {1};
  for (std::size_t m = 1; m <= n; m++) {
    std::vector<Residue> q(m + 1, 0);
    for (std::size_t k = 0; k < m; k++) {
      q[k + 1] = mod_add(q[k + 1], P[m - 1][k], p);
      q[k] = mod_sub(q[k], mod_mul(H[m - 1][m - 1], P[m - 1][k], p), p);
    }
    Residue t = 1;
    for (std::size_t i = m - 1; i >= 1; i--) {
      t = mod_mul(t, H[i][i - 1], p);
      const Residue c = mod_mul(H[i - 1][m - 1], t, p);
      if (c)
        for (std::size_t k = 0; k < P[i - 1].size(); k++) q[k] = mod_sub(q[k], mod_mul(c, P[i - 1][k], p), p);
    }
    P[m] = std::move(q);
  }
  return P[n];
}

inline Residue fp_evaluate(const std::vector<Residue>& poly, Residue x, Residue p) {
  Residue v = 0;
  for (std::size_t k = poly.size(); k-- > 0;) v = mod_add(mod_mul(v, x, p), poly[k], p);
  return v;
}

inline std::vector<std::uint64_t> prime_factors(std::uint64_t e) {
  std::vector<std::uint64_t> f;
  for (std::uint64_t q = 2; q * q <= e; q++)
    if (e % q == 0) {
      f.push_back(q);
      while (e % q == 0) e /= q;
    }
  if (e > 1) f.push_back(e);
  return f;
}

}  // namespace detail

inline CharacterTable character_table(const FiniteGroup& G) {
  using arithmetic::Residue;
  using detail::FpRows;
  CharacterTable T;
  const std::size_t n = G.order();
  T.order = n;
  T.classes = conjugacy_classes(G);
  const std::size_t r = T.classes.size();
  T.class_of.assign(n, 0);
  std::size_t e = 1;
  for (std::size_t j = 0; j < r; j++) {
    for (auto g : T.classes[j]) T.class_of[g] = j;
    T.class_sizes.push_back(T.classes[j].size());
    T.representatives.push_back(T.classes[j][0]);
    T.representative_orders.push_back(G.element_order(T.classes[j][0]));
    e = std::lcm(e, T.representative_orders.back());
  }
  for (std::size_t j = 0; j < r; j++) T.inverse_class.push_back(T.class_of[G.inv(T.representatives[j])]);
  T.exponent = e;

  // The least prime p ≡ 1 (mod e) with p > 2√|G|, and a primitive e-th root of unity z modulo p.
  Residue p = e + 1;
  while (!(arithmetic::is_prime(p) && p * p > 4 * static_cast<Residue>(n))) p += e;
  const auto qs = detail::prime_factors(e);
  // z = a^{(p−1)/e} has order e unless z^{e/q} = 1 for a prime q | e; some a in [1, p) succeeds.
  Residue z = 0;
  for (Residue a = 1; a < p && z == 0; a++) {
    const Residue c = arithmetic::mod_pow(a, (p - 1) / e, p);
    bool primitive = true;
    for (auto q : qs)
      if (arithmetic::mod_pow(c, e / q, p) == 1) primitive = false;
    if (primitive) z = c;
  }
  if (z == 0) throw std::logic_error("morphe::algebra::character_table: no primitive e-th root of unity modulo p");
  T.prime = p;
  T.root = z;

  // The common eigenspaces, each a basis in reduced row echelon form with its pivots.
  struct Space {
    FpRows basis;
    std::vector<std::size_t> pivots;
  };
  std::vector<Space> lines, todo;
  {
    Space all;
    for (std::size_t i = 0; i < r; i++) {
      all.basis.emplace_back(r, 0);
      all.basis.back()[i] = 1;
      all.pivots.push_back(i);
    }
    (r == 1 ? lines : todo).push_back(std::move(all));
  }
  // Small classes first: their matrices are the sparsest.
  std::vector<std::size_t> order(r > 0 ? r - 1 : 0);
  std::iota(order.begin(), order.end(), std::size_t(1));
  std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return T.class_sizes[x] < T.class_sizes[y]; });
  for (std::size_t j : order) {
    if (todo.empty()) break;
    std::vector<Residue> Mj = class_matrix(G, T, j);
    for (auto& x : Mj) x %= p;
    std::vector<Space> next;
    for (Space& S : todo) {
      const std::size_t d = S.basis.size();
      // M_j maps the space to itself; A is its matrix in the basis: column c holds the
      // coordinates of M_j v_c, which are its entries at the pivots.
      FpRows A(d, std::vector<Residue>(d, 0));
      for (std::size_t c = 0; c < d; c++)
        for (std::size_t i = 0; i < d; i++) {
          const std::size_t k = S.pivots[i];
          Residue s = 0;
          for (std::size_t l = 0; l < r; l++)
            if (S.basis[c][l]) s = arithmetic::mod_add(s, arithmetic::mod_mul(Mj[k * r + l], S.basis[c][l], p), p);
          A[i][c] = s;
        }
      const auto chi = detail::fp_charpoly(A, p);
      std::size_t found = 0;
      std::vector<Space> parts;
      for (Residue lambda = 0; lambda < p && found < d; lambda++) {
        if (detail::fp_evaluate(chi, lambda, p) != 0) continue;
        FpRows B = A;
        for (std::size_t i = 0; i < d; i++) B[i][i] = arithmetic::mod_sub(B[i][i], lambda, p);
        const FpRows U = detail::fp_kernel(B, d, p);
        Space part;
        for (const auto& u : U) {
          std::vector<Residue> v(r, 0);
          for (std::size_t c = 0; c < d; c++)
            if (u[c])
              for (std::size_t l = 0; l < r; l++) v[l] = arithmetic::mod_add(v[l], arithmetic::mod_mul(u[c], S.basis[c][l], p), p);
          part.basis.push_back(std::move(v));
        }
        part.pivots = detail::fp_rref(part.basis, p);
        found += part.basis.size();
        parts.push_back(std::move(part));
      }
      if (found != d) throw std::logic_error("morphe::algebra::character_table: a class matrix is not diagonalisable modulo p");
      for (auto& part : parts) (part.basis.size() == 1 ? lines : next).push_back(std::move(part));
    }
    todo = std::move(next);
  }
  if (!todo.empty() || lines.size() != r) throw std::logic_error("morphe::algebra::character_table: the class matrices do not separate the characters");

  // For each class, the classes of the powers g_j^i, i < o(g_j).
  std::vector<std::vector<std::size_t>> powers(r);
  for (std::size_t j = 0; j < r; j++) {
    std::uint32_t h = 0;
    for (std::size_t i = 0; i < T.representative_orders[j]; i++) {
      powers[j].push_back(T.class_of[h]);
      h = G.mul(h, T.representatives[j]);
    }
  }
  std::vector<Residue> size_inv(r);
  for (std::size_t j = 0; j < r; j++) size_inv[j] = detail::fp_inverse(T.class_sizes[j] % p, p);

  struct Character {
    std::size_t degree;
    std::vector<std::vector<long long>> multiplicities;
    std::vector<Cyclotomic> values;
  };
  std::vector<Character> chars;
  for (const Space& L : lines) {
    const auto& v = L.basis[0];
    if (v[0] == 0) throw std::logic_error("morphe::algebra::character_table: an eigenvector vanishes at the identity");
    const Residue v0 = detail::fp_inverse(v[0], p);
    std::vector<Residue> w(r);
    for (std::size_t j = 0; j < r; j++) w[j] = arithmetic::mod_mul(v[j], v0, p);
    // χ(1)² = |G| / Σ_j ω_j ω_{j*} / h_j.
    Residue s = 0;
    for (std::size_t j = 0; j < r; j++)
      s = arithmetic::mod_add(s, arithmetic::mod_mul(arithmetic::mod_mul(w[j], w[T.inverse_class[j]], p), size_inv[j], p), p);
    if (s == 0) throw std::logic_error("morphe::algebra::character_table: a vanishing norm modulo p");
    const Residue square = arithmetic::mod_mul(n % p, detail::fp_inverse(s, p), p);
    std::size_t degree = 0;
    for (std::size_t x = 1; x * x <= n; x++)
      if (arithmetic::mod_mul(x, x, p) == square) degree = x;
    if (degree == 0) throw std::logic_error("morphe::algebra::character_table: no degree in [1, √|G|]");
    std::vector<Residue> value(r);
    for (std::size_t j = 0; j < r; j++) value[j] = arithmetic::mod_mul(arithmetic::mod_mul(w[j], degree % p, p), size_inv[j], p);
    Character X;
    X.degree = degree;
    for (std::size_t j = 0; j < r; j++) {
      const std::size_t o = T.representative_orders[j];
      const Residue zo_inv = arithmetic::mod_pow(z, (e / o) * (o - 1), p);  // z_o⁻¹
      const Residue o_inv = detail::fp_inverse(o % p, p);
      std::vector<long long> m(o);
      std::vector<Rational> coeff(o);
      long long total = 0;
      for (std::size_t k = 0; k < o; k++) {
        const Residue step = arithmetic::mod_pow(zo_inv, k, p);
        Residue acc = 0, cur = 1;
        for (std::size_t i = 0; i < o; i++) {
          acc = arithmetic::mod_add(acc, arithmetic::mod_mul(value[powers[j][i]], cur, p), p);
          cur = arithmetic::mod_mul(cur, step, p);
        }
        const Residue mk = arithmetic::mod_mul(acc, o_inv, p);
        if (mk > degree) throw std::logic_error("morphe::algebra::character_table: a multiplicity outside [0, χ(1)]");
        m[k] = static_cast<long long>(mk);
        coeff[k] = Rational(m[k]);
        total += m[k];
      }
      if (total != static_cast<long long>(degree)) throw std::logic_error("morphe::algebra::character_table: the multiplicities do not sum to χ(1)");
      X.values.push_back(Cyclotomic::from_powers(static_cast<unsigned>(o), coeff));
      X.multiplicities.push_back(std::move(m));
    }
    chars.push_back(std::move(X));
  }
  std::sort(chars.begin(), chars.end(), [](const Character& x, const Character& y) {
    if (x.degree != y.degree) return x.degree < y.degree;
    return x.multiplicities > y.multiplicities;
  });
  for (auto& X : chars) T.values.push_back(std::move(X.values));
  return T;
}

// The classes of g_j^k for every class j, k any integer.
inline std::vector<std::size_t> power_map(const FiniteGroup& G, const CharacterTable& T, long long k) {
  std::vector<std::size_t> out;
  for (std::size_t j = 0; j < T.classes.size(); j++) {
    const long long o = static_cast<long long>(T.representative_orders[j]);
    const long long kk = ((k % o) + o) % o;
    std::uint32_t h = 0;
    for (long long i = 0; i < kk; i++) h = G.mul(h, T.representatives[j]);
    out.push_back(T.class_of[h]);
  }
  return out;
}

// ⟨f, g⟩ = (1/|G|) Σ_j h_j f(g_j) conj(g(g_j)) for class functions given on the classes.
inline Cyclotomic inner_product(const CharacterTable& T, const std::vector<Cyclotomic>& f, const std::vector<Cyclotomic>& g) {
  Cyclotomic s = 0;
  for (std::size_t j = 0; j < T.classes.size(); j++) s += Cyclotomic(static_cast<long long>(T.class_sizes[j])) * f[j] * g[j].conj();
  return s / Cyclotomic(static_cast<long long>(T.order));
}

// The multiplicities ⟨f, χ⟩ of the irreducible characters in a class function f.
inline std::vector<Cyclotomic> decompose(const CharacterTable& T, const std::vector<Cyclotomic>& f) {
  std::vector<Cyclotomic> m;
  for (std::size_t chi = 0; chi < T.size(); chi++) m.push_back(inner_product(T, f, T.values[chi]));
  return m;
}

// The Frobenius–Schur indicator (1/|G|) Σ_g χ(g²), which is 1, 0 or −1 as χ is the character of
// a real representation, is not real-valued, or is real-valued without a real representation
// (Isaacs 1976, Character Theory of Finite Groups, ch. 4). Throws if the sum were not rational.
inline Rational frobenius_schur(const FiniteGroup& G, const CharacterTable& T, std::size_t chi) {
  const auto sq = power_map(G, T, 2);
  Cyclotomic s = 0;
  for (std::size_t j = 0; j < T.classes.size(); j++) s += Cyclotomic(static_cast<long long>(T.class_sizes[j])) * T.values[chi][sq[j]];
  return (s / Cyclotomic(static_cast<long long>(T.order))).to_rational();
}

// The index of σ_k ∘ χ in the table, gcd(k, e) = 1: the Galois group of ℚ(ζ_e) permutes the
// irreducible characters. Throws std::logic_error if σ_k ∘ χ were not a row.
inline std::size_t galois_conjugate(const CharacterTable& T, std::size_t chi, long long k) {
  std::vector<Cyclotomic> row;
  for (const auto& x : T.values[chi]) row.push_back(x.galois(k));
  for (std::size_t psi = 0; psi < T.size(); psi++)
    if (T.values[psi] == row) return psi;
  throw std::logic_error("morphe::algebra::galois_conjugate: σ_k ∘ χ is not in the table");
}

}  // namespace morphe::algebra
