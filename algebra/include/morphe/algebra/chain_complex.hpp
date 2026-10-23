// Chain complexes of finitely generated free abelian groups, and their homology.
//
// A chain complex C is a sequence of free abelian groups C_k = ℤ^{n_k}, k = 0..d, with
// boundary maps ∂_k: C_k → C_{k-1} satisfying ∂_{k-1} ∂_k = 0, given by sparse integer
// matrices (rows: the basis of C_{k-1}; columns: that of C_k). Its homology is
// H_k = ker ∂_k / im ∂_{k+1}, a finitely generated abelian group ℤ^{b_k} ⊕ ⊕ ℤ/t_i.
//
// The computation has two phases.
//
// 1. Algebraic reduction (Kaczynski, Mischaikow and Mrozek 2004, Computational Homology,
//    §4.4). A pair (a, b), a a basis element of C_{k-1} and b of C_k with ⟨∂b, a⟩ = ±1, can
//    be removed: the complex with ∂'c = ∂c - (⟨∂c, a⟩ / ⟨∂b, a⟩) ∂b on the other cells of
//    dimension k, and the other maps restricted, is chain equivalent to C. The inclusion of
//    the reduced complex sends a k-cell c to c - (⟨∂c, a⟩ / ⟨∂b, a⟩) b, so the lift of each
//    surviving cell to an original chain is kept as the reduction proceeds, and cycles of
//    the reduced complex lift to cycles of C. Pivots are chosen to merge the shorter row into
//    the longer, which for ∂_1 contracts a spanning forest and for ∂_2 removes free faces.
//
// 2. The reduced complex, usually small, is made dense. With U_k ∂_k V_k = S_k in Smith
//    normal form of rank r_k, the last n_k - r_k columns of V_k are a basis of ker ∂_k; the
//    image of ∂_{k+1} in that basis is M = (V_k^{-1} ∂_{k+1}) restricted to its last rows,
//    and with P M Q = diag(e_1, ..., e_s), H_k ≅ ⊕ ℤ/e_i ⊕ ℤ^{n_k - r_k - s}. The columns of
//    K P^{-1} (K the kernel basis) are the generators: those with e_i > 1 generate the
//    torsion, those beyond s the free part (Munkres 1984, Elements of Algebraic Topology,
//    §11; Dumas, Saunders and Villard 2001, J. Symb. Comput. 32, 71–99).
//
// Cohomology follows from homology by the universal coefficient theorem,
// H^k ≅ ℤ^{b_k} ⊕ T_{k-1} (T the torsion), and Betti numbers over 𝔽_p from ranks mod p of
// the reduced complex, since unit pivots are units mod every p.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "matrix.hpp"
#include "smith.hpp"

namespace morphe::algebra {

using Chain = Vector;

struct ChainComplex {
  // n[k]: the rank of C_k, k = 0..d.
  std::vector<std::size_t> n;
  // boundary[k]: ∂_k, an n[k-1] × n[k] matrix, for k = 1..d; boundary[0] is unused.
  std::vector<SparseMatrix> boundary;

  std::size_t dimension() const { return n.empty() ? 0 : n.size() - 1; }
  // Throws std::invalid_argument unless the shapes agree and ∂_{k-1} ∂_k = 0.
  void validate() const {
    if (boundary.size() != n.size()) throw std::invalid_argument("morphe::algebra::ChainComplex: one boundary map per dimension");
    for (std::size_t k = 1; k < n.size(); k++) {
      if (boundary[k].rows != n[k - 1] || boundary[k].cols() != n[k])
        throw std::invalid_argument("morphe::algebra::ChainComplex: ∂_" + std::to_string(k) + " has the wrong shape");
      for (const auto& col : boundary[k].columns)
        for (const auto& e : col)
          if (e.index >= n[k - 1]) throw std::invalid_argument("morphe::algebra::ChainComplex: an entry outside C_" + std::to_string(k - 1));
      if (k >= 2) {
        const SparseMatrix dd = boundary[k - 1] * boundary[k];
        for (const auto& col : dd.columns)
          if (!col.empty()) throw std::invalid_argument("morphe::algebra::ChainComplex: ∂_" + std::to_string(k - 1) + " ∂_" + std::to_string(k) + " is not zero");
      }
    }
  }
  // The Euler characteristic Σ (-1)^k n_k.
  long long euler_characteristic() const {
    long long chi = 0;
    for (std::size_t k = 0; k < n.size(); k++) chi += (k % 2 ? -1 : 1) * static_cast<long long>(n[k]);
    return chi;
  }
};

struct HomologyGroup {
  std::size_t dimension = 0;
  std::size_t rank = 0;              // the Betti number b_k
  std::vector<Integer> torsion;      // the invariant factors > 1, each dividing the next
  std::vector<Chain> free_generators;     // cycles whose classes form a basis of the free part
  std::vector<Chain> torsion_generators;  // a cycle of order torsion[i] for each i
  std::string str() const {
    std::string s;
    if (rank) s = rank == 1 ? "ℤ" : "ℤ^" + std::to_string(rank);
    for (const auto& t : torsion) s += (s.empty() ? "" : " ⊕ ") + std::string("ℤ/") + t.str();
    return s.empty() ? "0" : s;
  }
};

struct CohomologyGroup {
  std::size_t dimension = 0, rank = 0;
  std::vector<Integer> torsion;
};

// The complex after algebraic reduction: dense boundaries on the surviving cells, and the
// lift of each surviving cell to a chain of the original complex.
struct ReducedComplex {
  std::vector<std::vector<std::size_t>> cells;  // cells[k]: the surviving original k-cells
  std::vector<Matrix> boundary;                  // boundary[k]: ∂'_k on the survivors, k >= 1
  std::vector<std::vector<Chain>> lift;          // lift[k][i]: survivor i of dimension k as an original chain
  std::size_t pairs = 0;                         // reduction pairs removed
};

inline ReducedComplex reduce(const ChainComplex& C, bool lifts = true) {
  const std::size_t d1 = C.n.size();
  std::vector<std::vector<Vector>> col(d1);                  // col[k][c]: ∂_k of cell c (current)
  std::vector<std::vector<std::vector<std::size_t>>> row(d1);  // row[k][r]: columns of ∂_k that may hold r
  std::vector<std::vector<char>> alive(d1);
  std::vector<std::vector<Chain>> lift(d1);
  for (std::size_t k = 0; k < d1; k++) {
    alive[k].assign(C.n[k], 1);
    if (lifts) {
      lift[k].resize(C.n[k]);
      for (std::size_t i = 0; i < C.n[k]; i++) lift[k][i] = {{i, 1}};
    }
    if (k == 0) continue;
    col[k] = C.boundary[k].columns;
    row[k].resize(C.n[k - 1]);
    for (std::size_t c = 0; c < C.n[k]; c++)
      for (const auto& e : col[k][c]) row[k][e.index].push_back(c);
  }
  ReducedComplex R;
  for (std::size_t k = 1; k < d1; k++) {
    bool progress = true;
    while (progress) {
      progress = false;
      for (std::size_t b = 0; b < C.n[k]; b++) {
        if (!alive[k][b]) continue;
        // A unit entry of column b in the row with the fewest columns.
        std::size_t a = SIZE_MAX, best = SIZE_MAX;
        Integer pivot;
        for (const auto& e : col[k][b])
          if ((e.value.is_one() || e.value == Integer(-1)) && row[k][e.index].size() < best) {
            a = e.index, best = row[k][e.index].size(), pivot = e.value;
          }
        if (a == SIZE_MAX) continue;
        // Clear row a by column operations c <- c - (⟨∂c, a⟩ / pivot) b.
        const Vector cb = col[k][b];
        std::vector<std::size_t> holders;
        holders.swap(row[k][a]);
        for (std::size_t c : holders) {
          if (c == b || !alive[k][c]) continue;
          const Integer v = coefficient(col[k][c], a);
          if (v.is_zero()) continue;
          const Integer lambda = v * pivot;  // v / pivot, pivot = ±1
          const Vector before = std::move(col[k][c]);
          col[k][c] = axpy(before, -lambda, cb);
          for (const auto& e : col[k][c])
            if (e.index != a && coefficient(before, e.index).is_zero()) row[k][e.index].push_back(c);
          if (lifts) lift[k][c] = axpy(lift[k][c], -lambda, lift[k][b]);
        }
        // Remove the pair: b from C_k, a from C_{k-1}.
        alive[k][b] = 0;
        alive[k - 1][a] = 0;
        for (const auto& e : cb)
          if (e.index != a) {
            auto& r = row[k][e.index];
            r.erase(std::remove(r.begin(), r.end(), b), r.end());
          }
        if (k + 1 < d1) {
          // b is no longer a basis element of C_k: delete its row from ∂_{k+1}.
          for (std::size_t c : row[k + 1][b])
            if (alive[k + 1][c]) {
              auto& v = col[k + 1][c];
              v.erase(std::remove_if(v.begin(), v.end(), [&](const Entry& e) { return e.index == b; }), v.end());
            }
          row[k + 1][b].clear();
        }
        R.pairs++;
        progress = true;
      }
    }
  }
  // The survivors, renumbered, and their dense boundaries.
  R.cells.resize(d1);
  R.boundary.resize(d1);
  R.lift.resize(d1);
  std::vector<std::vector<std::size_t>> position(d1);
  for (std::size_t k = 0; k < d1; k++) {
    position[k].assign(C.n[k], SIZE_MAX);
    for (std::size_t i = 0; i < C.n[k]; i++)
      if (alive[k][i]) {
        position[k][i] = R.cells[k].size();
        R.cells[k].push_back(i);
        if (lifts) R.lift[k].push_back(std::move(lift[k][i]));
      }
  }
  for (std::size_t k = 1; k < d1; k++) {
    R.boundary[k] = Matrix(R.cells[k - 1].size(), R.cells[k].size());
    for (std::size_t j = 0; j < R.cells[k].size(); j++)
      for (const auto& e : col[k][R.cells[k][j]]) {
        if (position[k - 1][e.index] == SIZE_MAX) throw std::logic_error("morphe::algebra::reduce: a boundary meets a removed cell");
        R.boundary[k](position[k - 1][e.index], j) = e.value;
      }
  }
  return R;
}

// The homology groups H_0 .. H_d, with generators when asked.
inline std::vector<HomologyGroup> homology(const ChainComplex& C, bool generators = true) {
  const std::size_t d1 = C.n.size();
  const ReducedComplex R = reduce(C, generators);
  std::vector<HomologyGroup> H(d1);
  // S[k]: the Smith form of ∂'_k, with V and V^-1 (kernel basis and coordinates in it).
  std::vector<Smith> S(d1 + 1);
  for (std::size_t k = 1; k < d1; k++) S[k] = smith(R.boundary[k], {.V = generators, .V_inverse = true});
  for (std::size_t k = 0; k < d1; k++) {
    const std::size_t nk = R.cells[k].size();
    const std::size_t rk = k == 0 ? 0 : S[k].rank();
    H[k].dimension = k;
    // The image of ∂'_{k+1} in the coordinates of the basis V_k of C'_k.
    Matrix X = k + 1 < d1 ? R.boundary[k + 1] : Matrix(nk, 0);
    if (k >= 1) X = S[k].V_inverse * X;
    Matrix M(nk - rk, X.cols);
    for (std::size_t i = rk; i < nk; i++)
      for (std::size_t j = 0; j < X.cols; j++) M(i - rk, j) = X(i, j);
    const Smith T = smith(M, {.U_inverse = generators});
    const std::size_t s = T.rank();
    H[k].rank = nk - rk - s;
    for (std::size_t i = 0; i < s; i++)
      if (!T.diagonal[i].is_one()) H[k].torsion.push_back(T.diagonal[i]);
    if (!generators) continue;
    // Generators: the columns of K P^-1, K the last nk - rk columns of V_k (identity for k = 0).
    auto generator = [&](std::size_t col) {
      Chain z;
      for (std::size_t j = 0; j < nk; j++) {
        // (K P^-1)(j, col) = Σ_i K(j, i) P^-1(i, col)
        Integer v = 0;
        for (std::size_t i = 0; i < nk - rk; i++) {
          const Integer& p = T.U_inverse(i, col);
          if (p.is_zero()) continue;
          const Integer kji = k == 0 ? Integer(j == i ? 1 : 0) : S[k].V(j, rk + i);
          if (!kji.is_zero()) v += kji * p;
        }
        if (!v.is_zero()) z = axpy(z, v, R.lift[k][j]);
      }
      return z;
    };
    for (std::size_t i = 0; i < s; i++)
      if (!T.diagonal[i].is_one()) H[k].torsion_generators.push_back(generator(i));
    for (std::size_t i = s; i < nk - rk; i++) H[k].free_generators.push_back(generator(i));
  }
  return H;
}

// Betti numbers over 𝔽_p (p prime) or, for p = 0, over ℚ.
inline std::vector<std::size_t> betti(const ChainComplex& C, std::uint64_t p = 0) {
  const ReducedComplex R = reduce(C, false);
  const std::size_t d1 = C.n.size();
  std::vector<std::size_t> r(d1 + 1, 0), b(d1);
  for (std::size_t k = 1; k < d1; k++) r[k] = rank_mod(R.boundary[k], p);
  for (std::size_t k = 0; k < d1; k++) b[k] = R.cells[k].size() - r[k] - r[k + 1];
  return b;
}

// Cohomology from homology by the universal coefficient theorem: H^k ≅ ℤ^{b_k} ⊕ T_{k-1}.
inline std::vector<CohomologyGroup> cohomology(const std::vector<HomologyGroup>& H) {
  std::vector<CohomologyGroup> co(H.size());
  for (std::size_t k = 0; k < H.size(); k++) {
    co[k].dimension = k;
    co[k].rank = H[k].rank;
    if (k > 0) co[k].torsion = H[k - 1].torsion;
  }
  return co;
}

}  // namespace morphe::algebra
