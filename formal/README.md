# formal: Words, Their Indexes, and Their Alignment over Semirings

*Status: draft (2026-10-24): the first slice, words and their indexes, semirings, the alignment of two words in any semiring, and the polytope algebra ([decision 0005](../development/decisions/0005-words-their-indexes-and-their-alignment.md)).* Automata beyond Aho–Corasick, grammars, transducers, rewriting systems, transition systems and bimorphisms are planned (morphe's README, §2).

## 1. Purpose

Formal structures: words and the structures built on them. This first slice serves the sequences of the sciences, the genomes of The Science Infrastructures' genomics domain first, and the texts of its linguistics next: the combinatorics of words, the indexes that make a long word searchable, and the alignment of two words, defined once as a sum over the paths of an edit graph in a semiring.

## 2. Concepts

- **Words and their morphisms.** A word is an element of the free monoid A* on an alphabet A = [0, σ); a morphism of monoids is fixed by the images of the letters (a substitution), and iterated it generates the classical infinite words, Thue–Morse, Fibonacci, period-doubling (Lothaire 1983, 2002). Borders and periods (Fine and Wilf), primitive roots (Lyndon and Schützenberger), conjugacy and least rotations are computed exactly, with the algorithms of Knuth, Morris and Pratt, Gusfield's Z-array and Booth.
- **Lyndon words.** A word strictly smaller than all its proper rotations; every word factors uniquely into a non-increasing product of Lyndon words (Chen, Fox and Lyndon 1958), found in linear time by Duval's algorithm, and their number is Witt's necklace polynomial, checked exactly against enumeration.
- **Suffix arrays and the Burrows–Wheeler transform.** The suffixes of a word sorted, by SA-IS in linear time (Nong, Zhang and Chan 2011); the longest common prefixes of neighbours (Kasai et al. 2001); the BWT, the last column of the sorted rotations, a permutation that groups equal letters and is inverted by the LF map. The bijective BWT (Gil and Scott, Kufleitner) needs no sentinel: it is the extended BWT of the word's Lyndon factors, a composite of three bijections, and the Gessel–Reutenauer correspondence reads multisets of necklaces off any word.
- **Indexes.** Bit vectors with rank in constant time and select; wavelet matrices; the FM-index, whose backward search counts the occurrences of a pattern in time proportional to its length (Ferragina and Manzini 2000); the Aho–Corasick automaton of a dictionary.
- **Semirings and the edit graph.** An alignment of two words is a path in their edit graph; a scoring scheme assigns its steps weights in a semiring S, and the sum over paths of the products of weights is, in (max, +), the optimal score; in ℕ, the number of alignments (the Delannoy numbers for global ones); in (max, +) with multiplicities, the optimum and the number of optimal alignments; in the polytope algebra, the Newton polygon whose vertices are every alignment optimal for some choice of parameters (Pachter and Sturmfels 2004). Semiring homomorphisms carry one reading to another: forgetting multiplicities, or evaluating a polygon's support function at a parameter, which is the tropical value there.
- **The polytope algebra.** Lattice polygons, with ⊕ the convex hull of a union and ⊗ the Minkowski sum, computed exactly in integer arithmetic; the normal fan of the alignment polygon divides the space of parameters into the cones where each vertex is optimal (Gusfield, Balasubramanian and Naor 1994; Waterman, Eggert and Lander 1992).

## 3. Content

| Header | Content |
|---|---|
| `words.hpp` | `Word`, `LetterRange`; the free monoid (`concat`, `power`, `rotation`) and its morphisms (`Morphism`, `compose` as g ∘ f, `iterate`, `fixed_point_prefix`; the Thue–Morse, Fibonacci and period-doubling words); `border_array` (the failure function of Knuth, Morris and Pratt), `borders`, `periods`, `z_array` (Gusfield), `find_all` (KMP); `primitive_root`, `is_primitive` (Lyndon–Schützenberger); `lyndon_factorisation` (Duval 1983), `is_lyndon`, `for_each_lyndon_word` (Duval 1988), `lyndon_count` (Witt's formula, exact); `least_rotation` (Booth 1980), `are_conjugate` |
| `suffix_array.hpp` | `suffix_array<Index>(w, σ)` by SA-IS, linear on integer alphabets (Nong, Zhang and Chan 2011), with `suffix_array_doubling` (Manber–Myers) and `suffix_array_naive` as references; `inverse_suffix_array`; `lcp_array` (Kasai et al. 2001); `distinct_substrings`, `factor_complexity`; `Concatenation`, `concatenate`, `generalised_suffix_array` for sets of words with distinct separators |
| `bwt.hpp` | `bwt` of w$ and `inverse_bwt` by LF (Burrows and Wheeler 1994); `lf_map`, `standard_permutation`, `sort_conjugates` in the ω-order; `extended_bwt` and its inverse (Mantaci, Restivo, Rosone and Sciortino 2007), `necklaces` (Gessel–Reutenauer); `rotation_bwt` and its inverse; `bijective_bwt` and its inverse (Gil and Scott 2012, Kufleitner 2009); `runs`, `move_to_front` |
| `succinct.hpp` | `BitVector` with rank in O(1) (Vigna's rank9) and select; `WaveletMatrix` with access, rank and select on integer sequences (Claude, Navarro and Ordóñez 2015) |
| `fm_index.hpp` | `FmIndex<Index>(text, σ, sample)` over any alphabet: backward search `interval` and `count`, `locate` by a sampled suffix array, `extract` by a sampled inverse (Ferragina and Manzini 2000) |
| `aho_corasick.hpp` | `AhoCorasick(patterns, σ)`: goto, failure and output functions, the deterministic automaton δ, and every occurrence in one pass (Aho and Corasick 1975) |
| `semiring.hpp` | Semirings as a C++20 concept; Boolean, (max, +) and (min, +) over integers (exact, or refused at their infinities), counting over ℕ on `Integer`, the tropical semiring with multiplicities, products; semiring homomorphisms (`SemiringMorphism`: forgetting multiplicities, the embedding of ℕ, negation, support, projections) and the laws checked on samples |
| `alignment.hpp` | The edit graph of two words and the sum over its paths in any semiring, global, local, semiglobal and overlap, linear and three-state affine schemes, traceback in selective semirings, enumeration of every path, Delannoy numbers |
| `polytope.hpp` | The polytope algebra of lattice polygons (exact hulls, Minkowski sums), support functions as homomorphisms to (max, +), normal fans, polytope propagation: the alignment polygon with a witness and a cone for each vertex |

## 4. Dependencies

`arithmetic` (`Integer`, for exact counts). No other library: morphe's invariant 1.

## 5. Checks

Every `tests/check_<topic>.cpp` is a check program `check_morphe_formal_<topic>`:

| Check | What it verifies |
|---|---|
| `check_words` | Morphisms are monoid homomorphisms; the closed forms of the Thue–Morse, Fibonacci and period-doubling words, and Thue–Morse overlap-free; borders, periods, Z-array and KMP against brute force; Fine and Wilf on every binary word up to length 14, with the bound attained; primitive roots and Lyndon–Schützenberger; Duval's factorisation against the greedy one; Duval's enumeration against Witt's formula and brute force, and Σ d L_k(d) = kⁿ exactly; Booth's least rotation and conjugacy |
| `check_suffix_array` | SA-IS against the naive sort (every binary word up to length 12; random words, σ 1 to 256, up to 10⁵ letters), unary closed forms, and prefix doubling on Fibonacci and Thue–Morse words of 10⁵ letters; Kasai's LCP; distinct factors and factor complexity, the Fibonacci word Sturmian and Thue–Morse matching Brlek's formula; sets of words with separators |
| `check_bwt` | The BWT of w$ against sorted rotations, and inverted; Burrows and Wheeler's example; the rotation BWT on every binary word up to length 12; the bijective BWT a bijection (binary up to 12, ternary up to 7); Gessel–Reutenauer on every binary word up to 10 and ternary up to 6; extended BWTs inverted; Mantaci, Restivo and Sciortino's characterisation of standard words in both directions; runs and move-to-front |
| `check_indexes` | Bit-vector rank and select; wavelet-matrix access, rank and select (σ 1 to 300); FM-index counts, positions and extracts against brute force (σ 1 to 256); Aho–Corasick against the paper's example and brute force on 300 dictionaries |
| `check_semiring` | The laws of every semiring on samples; the homomorphisms, and a map that is not one refused; exactness at the infinities |
| `check_alignment` | The generic programme against the enumeration of every path in five semirings, all modes and both gap models; Delannoy numbers up to D(60, 60) exactly; Durbin et al.'s example |
| `check_polytope` | Hulls, the laws of the algebra, Minkowski sums against pairwise sums, evaluation as a homomorphism; the alignment polygon equal to the hull of the enumerated vectors; witnesses and cones; the vertex of the longest common subsequence |
| `check_formal` | Every header of the module compiles in one translation unit |

## 6. References

- Aho, A. V. and Corasick, M. J. (1975), Efficient string matching, *Communications of the ACM* 18, 333–340.
- Burrows, M. and Wheeler, D. J. (1994), A block-sorting lossless data compression algorithm, DEC SRC research report 124.
- Claude, F., Navarro, G. and Ordóñez, A. (2015), The wavelet matrix, *Information Systems* 47, 15–32.
- Duval, J.-P. (1983), Factorizing words over an ordered alphabet, *Journal of Algorithms* 4, 363–381.
- Ferragina, P. and Manzini, G. (2000), Opportunistic data structures with applications, *FOCS*, 390–398.
- Gessel, I. M. and Reutenauer, C. (1993), Counting permutations with given cycle structure and descent set, *J. Combin. Theory A* 64, 189–215.
- Golan, J. S. (1999), *Semirings and their Applications*, Kluwer.
- Gusfield, D. (1997), *Algorithms on Strings, Trees, and Sequences*, Cambridge University Press.
- Gusfield, D., Balasubramanian, K. and Naor, D. (1994), Parametric optimization of sequence alignment, *Algorithmica* 12, 312–326.
- Kasai, T., Lee, G., Arimura, H., Arikawa, S. and Park, K. (2001), Linear-time longest-common-prefix computation in suffix arrays, *CPM*, LNCS 2089, 181–192.
- Lothaire, M. (1983), *Combinatorics on Words*, Addison-Wesley; (2002), *Algebraic Combinatorics on Words*, Cambridge University Press.
- Mantaci, S., Restivo, A., Rosone, G. and Sciortino, M. (2007), An extension of the Burrows–Wheeler transform, *Theoretical Computer Science* 387, 298–312.
- Mohri, M. (2002), Semiring frameworks and algorithms for shortest-distance problems, *J. Automata, Languages and Combinatorics* 7, 321–350.
- Nong, G., Zhang, S. and Chan, W. H. (2011), Two efficient algorithms for linear time suffix array construction, *IEEE Transactions on Computers* 60, 1471–1484.
- Pachter, L. and Sturmfels, B. (2004), Parametric inference for biological sequence analysis, *PNAS* 101, 16138–16143; (2005), *Algebraic Statistics for Computational Biology*, Cambridge University Press, ch. 2.
- Waterman, M. S., Eggert, M. and Lander, E. (1992), Parametric sequence comparisons, *PNAS* 89, 6090–6093.

## 7. Open Questions

1. **Three-dimensional polytopes**, for parametric alignment with affine gaps or in local mode.
2. **Automata and grammars** beyond Aho–Corasick, for linguistics: deterministic and weighted automata (whose shortest-distance problems are the same semiring sums), context-free grammars and their parsing.
3. **Large indexes.** The FM-index is in memory with 32- or 64-bit indices; external-memory construction would carry it to whole human genomes.
