# 0005. Words, Their Indexes, and Their Alignment over Semirings

- **Status:** Accepted
- **Date:** 2026-10-24
- **Decided by:** the user asked for a genome analysis infrastructure in The Science Infrastructures, with classical and advanced algorithms written at low level, and told another session that domains must reuse the foundations and never reimplement them. Claude placed the general mathematics of sequences here, as the first slice of the planned module `formal`, and wrote this record; the composition is The Science Infrastructures' [ADR 0043](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0043-the-genome-infrastructure.md).

## 1. Context

A genome is a word over four letters, read on two strands; a text, which linguistics will study next, is a word over many. Both need the combinatorics of words and the indexes that make long words searchable, and both need the alignment of two words. None of this names an object in nature, and the map of the mathematics (the infrastructure's `docs/mathematics.md`, §3 and §8) gives words, automata and grammars to morphe's `formal`, which was planned and empty.

The alignment of two words is, mathematically, a sum over the paths of a directed acyclic graph of products of weights in a semiring: the optimal score in (max, +), the number of alignments in ℕ, the number of optimal ones in (max, +) with multiplicities, the probability of a pair HMM in (+, ×), and every optimal alignment for every scoring in the polytope algebra (Pachter and Sturmfels 2004). The genomics domain's fast algorithms are schedules of that sum; its definition is mathematics, and exact.

## 2. Decision

1. **The module `formal` is created** with its first slice, `morphe::formal`, header-only, status `draft`, depending on `arithmetic`:
   - words and their morphisms, borders and periods, primitivity, Lyndon factorisation and enumeration, least rotations (`words.hpp`);
   - suffix arrays by SA-IS with references, LCP, factor complexity, sets of words (`suffix_array.hpp`);
   - the Burrows–Wheeler transform and its inverse, the extended and bijective transforms, Gessel–Reutenauer, runs and move-to-front (`bwt.hpp`);
   - rank and select, wavelet matrices (`succinct.hpp`); the FM-index over any alphabet (`fm_index.hpp`); the Aho–Corasick automaton (`aho_corasick.hpp`);
   - semirings as a concept, the exact ones, and their homomorphisms (`semiring.hpp`);
   - the edit graph of two words and the sum over its paths in any semiring, four modes, linear and affine schemes, traceback in selective semirings, enumeration (`alignment.hpp`);
   - the polytope algebra of lattice polygons, normal fans, and polytope propagation for parametric alignment (`polytope.hpp`).
2. **Words are integer sequences.** A word is a `std::vector<std::uint32_t>` (`Word`), and every algorithm that indexes or aligns words accepts any contiguous range of unsigned integers (`LetterRange`), so a domain passes its own byte codes without copying.
3. **Exact by default** (invariant 3): integer weights in the tropical semirings refuse to reach their infinities rather than wrap; counts are `Integer`; polygons have integer vertices and exact cross products. Floating point appears only where a domain instantiates a probability semiring (the genomics domain's pair HMM).
4. **Eulerian circuits and determinants.** The assembly of genomes counts the texts that share a k-mer composition, which is the number of Eulerian trails of a multigraph: `combinatorics/eulerian.hpp` (Hierholzer's circuits and trails, the matrix-tree theorem, the BEST theorem) and the exact determinant it needs, `algebra/determinant.hpp` (Bareiss). `combinatorics` therefore depends on `algebra` and `arithmetic`; the order stays acyclic, since `algebra` does not depend on `combinatorics`, and `topology` already stood on both.
5. **Every check program of the module is registered by a glob** (`tests/check_<topic>.cpp` → `check_morphe_formal_<topic>`), and one of them includes every header in one translation unit, so that two parts written apart cannot collide (a concept and a struct of one name once did).

## 3. Alternatives Considered

- **The words in the genomics domain.** Rejected: linguistics needs them too, and the user's rule forbids a domain to hold what a foundation should.
- **The alignment of words in the domain**, with only its fast schedules. Rejected: the definition is a semiring sum, exact and general; the domain keeps the schedules, the biological scorings and the probability semirings, and holds every schedule to this definition.
- **calx for the indexes**, as algorithms on arrays. Rejected: calx is numerical; these are exact discrete structures with statements of termination and complexity, which is morphe's mode.

## 4. Consequences

- morphe's README, CLAUDE.md §4 and the architecture list `formal` as written in part.
- The genomics domain's `index` builds its FMD-index on SA-IS, its `align` holds Gotoh, Myers–Miller, banded, WFA and the striped kernels to `alignment.hpp`, and its parametric alignment is `polytope.hpp`'s propagation; its `annotation` finds repeats by the suffix array and LCP.
- **Open:** three-dimensional polytopes (affine or local parametric alignment); weighted automata, whose shortest-distance problems are the same semiring sums (Mohri 2002); grammars and their parsing, for linguistics; external-memory construction of indexes for whole human genomes.
