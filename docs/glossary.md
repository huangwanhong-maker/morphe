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
