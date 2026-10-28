# Glossary

*Status: draft, 2026-10-23.* The conventions morphe fixes once (CLAUDE.md §6), and the notation of its modules with the names in the code. Created with the first slice ([decision 0004](../development/decisions/0004-the-first-slice-cell-complexes.md)).

## 1. Conventions Fixed Once

| Convention | Choice | Where it shows |
|---|---|---|
| Order of composition | Applicative: g ∘ f applies f first. Permutations compose so: (p q)(i) = p(q(i)). | `permutation_group` |
| Words in groups | Read left to right: the word (1, 2) is x₁ x₂; letters are ±i, 1-based, −i the inverse of x_i. | `Word`, `Presentation` |
| Group actions | Left actions, unless a construction states otherwise. | — |
| Gauge transformations (physis) | U_e ↦ g_tail U_e g_head⁻¹, a holonomy multiplying along the word. | physis `gauge` |
| Orientation of a simplex | By increasing vertex label: ∂[v₀ … v_k] = Σ (−1)^i [v₀ … v̂_i … v_k]. | `SimplicialComplex` |
| Boundary of an edge | head − tail. | `CWComplex::chain_complex` |
| Strict or up to isomorphism | Constructions are strict (they return explicit cells); every comparison of results states its category (isomorphism of groups, homotopy, simple homotopy). | — |
| Division of integers | Truncated toward zero; `mod` and `floor_div` give the Euclidean pair. | `Integer` |
| Fourier transform on a finite group | f̂(ρ) = Σ_g f(g) ρ(g), with no conjugate (Diaconis 1988), so that the product of K[G] goes to the product of matrices; for C_{n_1} × ⋯ × C_{n_k}, f̂(y) = Σ_x f(x) ζ_N^{⟨x, y⟩}, the transform with kernel exp(−2πi x y/n) being f̂(−y). | `SnFourier`, `AbelianFourier` |
| Complex values of cyclotomic numbers | ζ_m = exp(2πi/m), for display only; exact work never leaves ℚ(ζ_m). | `Cyclotomic::to_complex` |

## 2. Notation

| Symbol | Meaning | Code |
|---|---|---|
| ℤ, ℚ, 𝔽_p | integers, rationals, the prime field | `Integer`, `Rational`, `Residue` |
| U A V = D | Smith normal form; d₁ ∣ d₂ ∣ … the invariant factors | `smith`, `Smith::diagonal` |
| C_k, ∂_k | chains of dimension k, the boundary map C_k → C_{k−1} | `ChainComplex::n`, `ChainComplex::boundary` |
| H_k, b_k, T_k | homology, Betti number, torsion | `HomologyGroup::rank`, `::torsion` |
| H^k | cohomology, ℤ^{b_k} ⊕ T_{k−1} | `cohomology` |
| χ | Euler characteristic Σ (−1)^k n_k | `euler_characteristic` |
| ⟨X ∣ R⟩ | a presentation | `Presentation` |
| G^{ab} | abelianisation | `abelianisation` |
| π₁(X, x) | fundamental group | `fundamental_group` |
| e^{±1} | a step along an edge | `Step{edge, sign}` |
| g, k, b | orientable genus, crosscaps, boundary circles | `SurfaceType` |
| K*, the dual 2-complex | of a pseudomanifold K | `dual_two_complex` |
| σ_c | the permutation of colour c of a coloured graph | `ColouredGraph::sigma` |
| g_J, ω | genus of a jacket, Gurau's degree Σ g_J | `jacket_genera`, `degree` |
| X_s, s(c) | the sub-complex at stage s of a filtration, the stage of a cell | `Filtration` |
| [b, d) | an interval of persistence: a class of H_k over ℤ/2 born at stage b and dying at d (or never) | `PersistenceInterval` |
| ℚ(ζ_m), Φ_m, φ(m) | the m-th cyclotomic field ℚ[x]/(Φ_m), the m-th cyclotomic polynomial, Euler's totient | `Cyclotomic`, `cyclotomic_polynomial`, `euler_phi` |
| σ_k | the automorphism ζ ↦ ζ^k of ℚ(ζ_m), gcd(k, m) = 1; σ_{−1} is complex conjugation | `Cyclotomic::galois`, `Cyclotomic::conj` |
| K_j, h_j, g_j, C_j | a conjugacy class, its size, its representative (its least element), its class sum | `CharacterTable::classes`, `::class_sizes`, `::representatives` |
| χ (of a group), χ(1) | an irreducible complex character, its degree (not the Euler characteristic) | `CharacterTable::values`, `::degree` |
| ω_χ, a_{jkl} | the central character h_j χ(g_j)/χ(1); the class coefficients, C_j C_k = Σ_l a_{jkl} C_l | `class_multiplication_coefficients` |
| ν(χ) | the Frobenius–Schur indicator (1/\|G\|) Σ_g χ(g²) | `frobenius_schur` |
| K[G], e_χ | the group algebra; the central primitive idempotent (χ(1)/\|G\|) Σ_g χ(g⁻¹) g | `GroupAlgebra`, `central_idempotent` |
| λ ⊢ n, ρ_λ, d_λ | a partition of n, the irreducible representation of S_n it indexes (Young's seminormal form), its dimension | `Partition`, `SeminormalRep`, `hook_length_dimension` |
| s_i, c_T(k) | the Coxeter generator (i i+1) of S_n, 0-based; the content (column − row) of the cell of k in a tableau T | `SeminormalRep::generator`, `::content` |
| f̂(λ), E^λ_{ST} | the Fourier transform Σ_σ f(σ) ρ_λ(σ); the matrix units of ℚ[S_n], inverse transforms of e_{ST} | `SnFourier::transform`, `::matrix_unit` |
| ⟨x, y⟩, χ_y | the pairing Σ_i x_i y_i N/n_i on C_{n_1} × ⋯ × C_{n_k}, N = lcm n_i; the character ζ_N^{⟨x, y⟩} | `AbelianFourier::pairing`, `::character` |
