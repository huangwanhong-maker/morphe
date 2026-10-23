// Checks of morphe/algebra: the Smith normal form (U A V = D with unimodular transforms and
// their inverses, the divisibility chain, invariant factors against gcds of minors, a
// textbook example); ranks modulo primes; homology of chain complexes disguised by random
// unimodular changes of basis, with generators that are cycles, independent modulo
// boundaries and of the stated orders; presentations (Tietze simplification, abelianisation);
// finite groups (orders, classes) and homomorphisms from surface groups, counted against the
// formulas of Frobenius and Mednykh.

#include <morphe/algebra/chain_complex.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/algebra/presentation.hpp>
#include <morphe/algebra/smith.hpp>
#include <morphe/check.hpp>

#include <random>
#include <string>

namespace ma = morphe::algebra;
using morphe::arithmetic::Integer;

namespace {

ma::Matrix random_matrix(std::mt19937_64& g, std::size_t r, std::size_t c, int spread) {
  ma::Matrix A(r, c);
  for (auto& x : A.a) x = static_cast<long long>(g() % (2 * spread + 1)) - spread;
  return A;
}

// A random unimodular n × n matrix and its inverse, as a product of elementary operations.
std::pair<ma::Matrix, ma::Matrix> random_unimodular(std::mt19937_64& g, std::size_t n, int steps) {
  ma::Matrix B = ma::Matrix::identity(n), Bi = ma::Matrix::identity(n);
  if (n < 2) return {B, Bi};
  for (int s = 0; s < steps; s++) {
    const std::size_t i = g() % n, j = (i + 1 + g() % (n - 1)) % n;
    const long long c = static_cast<long long>(g() % 5) - 2;
    // B <- E B with E = I + c e_i e_j^T; Bi <- Bi E^-1.
    for (std::size_t k = 0; k < n; k++) B(i, k) += Integer(c) * B(j, k);
    for (std::size_t k = 0; k < n; k++) Bi(k, j) -= Integer(c) * Bi(k, i);
  }
  return {B, Bi};
}

Integer gcd_of_minors2(const ma::Matrix& A) {
  Integer g = 0;
  for (std::size_t i = 0; i < A.rows; i++)
    for (std::size_t k = i + 1; k < A.rows; k++)
      for (std::size_t j = 0; j < A.cols; j++)
        for (std::size_t l = j + 1; l < A.cols; l++) g = morphe::arithmetic::gcd(g, A(i, j) * A(k, l) - A(i, l) * A(k, j));
  return g;
}

ma::Word word(std::initializer_list<int> w) { return ma::Word(w); }

// ⟨a1, b1, ..., ag, bg | [a1, b1] ... [ag, bg]⟩.
ma::Presentation surface_group(int g) {
  ma::Presentation P;
  P.generators = static_cast<std::size_t>(2 * g);
  ma::Word r;
  for (int i = 0; i < g; i++) {
    const int a = 2 * i + 1, b = 2 * i + 2;
    r.insert(r.end(), {a, b, -a, -b});
  }
  P.relators = {r};
  return P;
}
// ⟨a1..ak | a1 a1 ... ak ak⟩.
ma::Presentation crosscap_group(int k) {
  ma::Presentation P;
  P.generators = static_cast<std::size_t>(k);
  ma::Word r;
  for (int i = 1; i <= k; i++) r.insert(r.end(), {i, i});
  P.relators = {r};
  return P;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(17);

  // ---- the Smith normal form -----------------------------------------------------------------------------------------
  {
    const ma::Matrix A = [] {
      ma::Matrix m(3, 3);
      const long long v[9] = {2, 4, 4, -6, 6, 12, 10, -4, -16};
      for (int i = 0; i < 9; i++) m.a[static_cast<std::size_t>(i)] = v[i];
      return m;
    }();
    const auto S = ma::smith(A);
    s.require("algebra: the Smith form of [[2,4,4],[-6,6,12],[10,-4,-16]] is diag(2, 6, 12)",
              S.diagonal == std::vector<Integer>{2, 6, 12});
    int bad = 0, chain = 0, minors = 0;
    for (int t = 0; t < 300; t++) {
      const std::size_t r = 1 + g() % 6, c = 1 + g() % 6;
      const ma::Matrix B = random_matrix(g, r, c, t % 3 == 0 ? 1 : 9);
      const auto T = ma::smith(B, {.U = true, .V = true, .U_inverse = true, .V_inverse = true});
      if (T.U * B * T.V != T.D) bad++;
      if (T.U * T.U_inverse != ma::Matrix::identity(r) || T.V * T.V_inverse != ma::Matrix::identity(c)) bad++;
      for (std::size_t i = 0; i < r; i++)
        for (std::size_t j = 0; j < c; j++)
          if (!T.D(i, j).is_zero() && (i != j || i >= T.rank())) bad++;
      for (std::size_t i = 0; i + 1 < T.rank(); i++)
        if (!(T.diagonal[i + 1] % T.diagonal[i]).is_zero()) chain++;
      // d1 = gcd of the entries, d1 d2 = gcd of the 2 × 2 minors.
      Integer g1 = 0;
      for (const auto& x : B.a) g1 = morphe::arithmetic::gcd(g1, x);
      const Integer d1 = T.rank() >= 1 ? T.diagonal[0] : Integer(0);
      const Integer d12 = T.rank() >= 2 ? T.diagonal[0] * T.diagonal[1] : Integer(0);
      if (d1 != g1 || d12 != gcd_of_minors2(B)) minors++;
    }
    s.equal("algebra: U A V = D, U U^-1 = I, V V^-1 = I, D diagonal, on 300 random matrices up to 6 × 6 (failures)", bad, 0);
    s.equal("algebra: each invariant factor divides the next (failures)", chain, 0);
    s.equal("algebra: d1 is the gcd of the entries and d1 d2 the gcd of the 2 × 2 minors (failures)", minors, 0);
    ma::Matrix six(2, 2);
    six(0, 0) = 2, six(1, 1) = 3;
    s.require("algebra: the rank of diag(2, 3) is 1 mod 2, 1 mod 3, 2 mod 5 and over ℚ",
              ma::rank_mod(six, 2) == 1 && ma::rank_mod(six, 3) == 1 && ma::rank_mod(six, 5) == 2 && ma::rank_mod(six, 0) == 2);
  }
  // ---- homology of disguised complexes --------------------------------------------------------------------------------
  {
    // A direct sum of elementary complexes ℤ --d--> ℤ and free ℤ's, then a random change of
    // basis in every degree: ∂'_k = B_{k-1} ∂_k B_k^{-1}. Homology is known by construction.
    int wrong = 0, not_cycles = 0, dependent = 0;
    for (int t = 0; t < 40; t++) {
      const std::size_t d = 1 + g() % 3;
      std::vector<std::size_t> rank(d + 1, 0);
      std::vector<std::vector<Integer>> tors(d + 1);
      std::vector<std::size_t> n(d + 1, 0);
      struct Block { std::size_t k; long long factor; };
      std::vector<Block> blocks;
      for (std::size_t k = 0; k <= d; k++) {
        const std::size_t free = g() % 3;
        rank[k] += free;
        n[k] += free;
        for (std::size_t i = 0; i < free; i++) blocks.push_back({k, 0});
        if (k >= 1)
          for (int b = 0; b < static_cast<int>(g() % 3); b++) {
            const long long f = std::vector<long long>{1, 1, 2, 3, 4, 6}[g() % 6];
            blocks.push_back({k, f});
            n[k]++, n[k - 1]++;
            if (f > 1) tors[k - 1].push_back(f);
          }
      }
      // Build the diagonal complex.
      std::vector<ma::Matrix> D(d + 1);
      for (std::size_t k = 1; k <= d; k++) D[k] = ma::Matrix(n[k - 1], n[k]);
      std::vector<std::size_t> used(d + 1, 0);
      for (const auto& b : blocks) {
        if (b.factor == 0) { used[b.k]++; continue; }
        D[b.k](used[b.k - 1]++, used[b.k]++) = b.factor;
      }
      std::vector<std::pair<ma::Matrix, ma::Matrix>> B(d + 1);
      for (std::size_t k = 0; k <= d; k++) B[k] = random_unimodular(g, n[k], 30);
      ma::ChainComplex C;
      C.n = n;
      C.boundary.resize(d + 1);
      for (std::size_t k = 1; k <= d; k++) C.boundary[k] = ma::SparseMatrix::from(B[k - 1].first * D[k] * B[k].second);
      C.validate();
      const auto H = ma::homology(C);
      for (std::size_t k = 0; k <= d; k++) {
        std::sort(tors[k].begin(), tors[k].end());
        // Invariant factors of ⊕ ℤ/t_i are not the t_i themselves; compare the group orders and
        // the number of cyclic factors of each prime power by comparing with the Smith form.
        ma::Matrix T(tors[k].size(), tors[k].size());
        for (std::size_t i = 0; i < tors[k].size(); i++) T(i, i) = tors[k][i];
        const auto want = ma::torsion(ma::smith(T));
        if (H[k].rank != rank[k] || H[k].torsion != want) wrong++;
        // Generators are cycles.
        for (const auto* gens : {&H[k].free_generators, &H[k].torsion_generators})
          for (const auto& z : *gens)
            if (k >= 1 && !C.boundary[k].apply(z).empty()) not_cycles++;
        // Free generators are independent modulo boundaries: rank over ℚ of [im ∂_{k+1} | gens].
        const std::size_t cols_im = k < d ? n[k + 1] : 0;
        ma::Matrix M(n[k], cols_im + H[k].free_generators.size());
        if (k < d)
          for (std::size_t j = 0; j < n[k + 1]; j++)
            for (const auto& e : C.boundary[k + 1].columns[j]) M(e.index, j) = e.value;
        for (std::size_t j = 0; j < H[k].free_generators.size(); j++)
          for (const auto& e : H[k].free_generators[j]) M(e.index, cols_im + j) = e.value;
        ma::Matrix Im(n[k], cols_im);
        for (std::size_t i = 0; i < n[k]; i++)
          for (std::size_t j = 0; j < cols_im; j++) Im(i, j) = M(i, j);
        if (ma::rank_mod(M, 0) != ma::rank_mod(Im, 0) + H[k].free_generators.size()) dependent++;
      }
      // Betti numbers mod 2 by the universal coefficient theorem: b_k + #even t in T_k and T_{k-1}.
      const auto b2 = ma::betti(C, 2);
      for (std::size_t k = 0; k <= d; k++) {
        std::size_t want = rank[k];
        for (const auto& t : H[k].torsion) want += (t % 2).is_zero();
        if (k >= 1)
          for (const auto& t : H[k - 1].torsion) want += (t % 2).is_zero();
        if (b2[k] != want) wrong++;
      }
    }
    s.equal("algebra: homology of 40 disguised complexes (ranks, torsion, Betti numbers mod 2) as constructed (failures)", wrong, 0);
    s.equal("algebra: every generator is a cycle (failures)", not_cycles, 0);
    s.equal("algebra: free generators are independent modulo boundaries (failures)", dependent, 0);
  }
  {
    // ℝP² as a CW complex: one cell in each dimension, ∂_2 = 2, ∂_1 = 0.
    ma::ChainComplex C;
    C.n = {1, 1, 1};
    C.boundary.resize(3);
    C.boundary[1] = ma::SparseMatrix(1, 1);
    C.boundary[2] = ma::SparseMatrix(1, 1);
    C.boundary[2].add(0, 0, 2);
    const auto H = ma::homology(C);
    s.require("algebra: H(ℝP²) = ℤ, ℤ/2, 0, the torsion generated by the 1-cell", H[0].str() == "ℤ" && H[1].str() == "ℤ/2" && H[2].str() == "0" &&
                                                                         H[1].torsion_generators.size() == 1 && H[1].torsion_generators[0] == ma::Chain{{0, 1}});
    const auto co = ma::cohomology(H);
    s.require("algebra: H^*(ℝP²) = ℤ, 0, ℤ/2 by universal coefficients", co[0].rank == 1 && co[1].rank == 0 && co[1].torsion.empty() && co[2].torsion == std::vector<Integer>{2});
  }
  // ---- presentations ---------------------------------------------------------------------------------------------------
  {
    const auto T = ma::simplify(surface_group(1));
    s.require("algebra: the torus group ⟨a, b | a b a⁻¹ b⁻¹⟩ is left as it is", T.presentation.generators == 2 && T.presentation.relators.size() == 1);
    ma::Presentation P;
    P.generators = 3;
    P.names = {"a", "b", "c"};
    P.relators = {word({1, 2, -3}), word({3, 3, 3}), word({2})};  // c = a b, c³ = 1, b = 1
    const auto S = ma::simplify(P);
    s.require("algebra: ⟨a, b, c | a b c⁻¹, c³, b⟩ simplifies to one generator and one relator of length 3",
              S.presentation.generators == 1 && S.presentation.relators.size() == 1 && S.presentation.relators[0].size() == 3, S.presentation.str());
    s.equal("algebra: ... and abelianises to ℤ/3", ma::abelianisation(S.presentation).str(), std::string("ℤ/3"));
    s.equal("algebra: the abelianised genus-3 surface group is ℤ^6", ma::abelianisation(surface_group(3)).str(), std::string("ℤ^6"));
    s.equal("algebra: the abelianised Klein bottle group ⟨a, b | a a b b⟩ is ℤ ⊕ ℤ/2", ma::abelianisation(crosscap_group(2)).str(), std::string("ℤ ⊕ ℤ/2"));
    s.equal("algebra: free reduction of a b b⁻¹ a⁻¹ c", ma::reduce(word({1, 2, -2, -1, 3})), word({3}));
    s.equal("algebra: cyclic reduction of a⁻¹ b a", ma::cyclic_reduce(word({-1, 2, 1})), word({2}));
  }
  // ---- finite groups ---------------------------------------------------------------------------------------------------
  {
    struct Known { ma::FiniteGroup G; std::size_t order, classes; std::vector<int> dims; std::vector<int> fs; };
    // Degrees of the irreducible representations and their Frobenius–Schur indicators (James and Liebeck 2001).
    std::vector<Known> groups = {
        {ma::cyclic(5), 5, 5, {1, 1, 1, 1, 1}, {1, 0, 0, 0, 0}},
        {ma::symmetric(3), 6, 3, {1, 1, 2}, {1, 1, 1}},
        {ma::quaternion(), 8, 5, {1, 1, 1, 1, 2}, {1, 1, 1, 1, -1}},
        {ma::dihedral(4), 8, 5, {1, 1, 1, 1, 2}, {1, 1, 1, 1, 1}},
        {ma::alternating(4), 12, 4, {1, 1, 1, 3}, {1, 0, 0, 1}},
        {ma::symmetric(4), 24, 5, {1, 1, 2, 3, 3}, {1, 1, 1, 1, 1}},
        {ma::alternating(5), 60, 5, {1, 3, 3, 4, 5}, {1, 1, 1, 1, 1}},
    };
    int bad = 0, hom = 0;
    std::string detail;
    for (const auto& k : groups) {
      if (k.G.order() != k.order || ma::conjugacy_classes(k.G).size() != k.classes) bad++, detail += k.G.name + " ";
      // Orientable surfaces: |Hom(π1(Σ_g), G)| = |G|^{2g-1} Σ_ρ d_ρ^{2-2g}, g = 1, 2.
      for (int genus = 1; genus <= 2; genus++) {
        if (genus == 2 && k.order > 24) continue;
        // Σ_ρ |G|^{2g-1} / d^{2g-2}, exact since d divides |G|.
        Integer want = 0;
        for (int dim : k.dims)
          want += morphe::arithmetic::exact_div(morphe::arithmetic::pow(Integer(static_cast<long long>(k.order)), 2 * genus - 1),
                                                morphe::arithmetic::pow(Integer(dim), 2 * genus - 2));
        if (Integer(ma::count_homomorphisms(surface_group(genus), k.G)) != want) hom++, detail += k.G.name + "/g" + std::to_string(genus) + " ";
      }
      // Non-orientable surfaces N_1 = ℝP², N_2 = Klein bottle: |G|^{k-1} Σ_ρ ν_ρ^k d_ρ^{2-k}.
      for (int cc = 1; cc <= 2; cc++) {
        long long want = 0;
        for (std::size_t r = 0; r < k.dims.size(); r++) {
          const long long nu = k.fs[r], dk = cc == 1 ? k.dims[r] : 1;
          want += (cc == 1 ? nu : nu * nu) * dk;
        }
        if (cc == 2) want *= static_cast<long long>(k.order);
        if (ma::count_homomorphisms(crosscap_group(cc), k.G) != static_cast<std::uint64_t>(want)) hom++, detail += k.G.name + "/N" + std::to_string(cc) + " ";
      }
    }
    s.require("algebra: orders and class numbers of Z5, S3, Q8, D4, A4, S4, A5", bad == 0, detail);
    s.require("algebra: homomorphisms from surface groups (torus, genus 2, ℝP², Klein bottle) agree with Frobenius and Mednykh", hom == 0, detail);
    s.equal("algebra: Z2 × Z3 is abelian of order 6 with elements of order 6", ma::direct_product(ma::cyclic(2), ma::cyclic(3)).element_order(4), 6u);
  }
  return s.done();
}
