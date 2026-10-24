// Checks of morphe/formal/alignment.hpp: on every pair of short random words, in the four modes,
// with linear and with affine gaps, the dynamic programme against the enumeration of every
// path: the optimal score (max, +), the least cost (min, +), the number of paths (Counting) and
// the number of optimal paths (the tropical semiring with multiplicities), whether a path exists
// (Boolean); every traceback a path the mode allows, of the optimal weight; forgetting the
// multiplicities of the sum is the sum of the forgotten weights (a homomorphism commutes with the
// dynamic programme); affine gaps with a zero open are linear gaps; the number of global
// alignments the Delannoy number, against the central values of OEIS A001850 and, beyond 64
// bits, the closed form; Durbin et al.'s HEAGAWGHEE and PAWHEAE with BLOSUM50 and gaps of −8.

#include <morphe/check.hpp>
#include <morphe/formal/alignment.hpp>
#include <morphe/formal/semiring.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace mf = morphe::formal;
using morphe::arithmetic::Integer;

namespace {

// splitmix64 (Steele, Lea and Flood 2014), for reproducible samples without a dependency.
struct Random {
  std::uint64_t state;
  std::uint64_t next() {
    std::uint64_t z = (state += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(next() % n); }
};

std::vector<std::uint8_t> word(Random& r, std::size_t n, std::size_t sigma) {
  std::vector<std::uint8_t> w(n);
  for (auto& c : w) c = static_cast<std::uint8_t>(r.below(sigma));
  return w;
}

bool allowed(const mf::EditPath& p, std::size_t m, std::size_t n, mf::Mode mode) {
  return mf::is_source(mode, p.x_begin, p.y_begin) && p.x_end() <= m && p.y_end() <= n && mf::is_sink(mode, p.x_end(), p.y_end(), m, n);
}

const char* name(mf::Mode m) {
  switch (m) {
    case mf::Mode::global: return "global";
    case mf::Mode::local: return "local";
    case mf::Mode::semiglobal: return "semiglobal";
    case mf::Mode::overlap: return "overlap";
  }
  return "?";
}

}  // namespace

int main() {
  morphe::Suite s;
  Random r{20261024};
  const mf::Mode modes[] = {mf::Mode::global, mf::Mode::local, mf::Mode::semiglobal, mf::Mode::overlap};

  // Scores: match 2, mismatch −3 (letters differ), gap letters −2, open −3 (affine); costs for
  // (min, +): mismatch 1, gap 1 (edit distance), open 2 (affine).
  const auto pair = [](std::uint8_t a, std::uint8_t b) { return std::int64_t(a == b ? 2 : -3); };
  const auto gap = [](std::uint8_t) { return std::int64_t(-2); };
  const auto cost = [](std::uint8_t a, std::uint8_t b) { return std::int64_t(a == b ? 0 : 1); };
  const auto unit = [](std::uint8_t) { return std::int64_t(1); };
  const auto one_count = [](auto...) { return Integer(1); };
  const auto counted_pair = [&](std::uint8_t a, std::uint8_t b) { return mf::MaxPlusCount::single(pair(a, b)); };
  const auto counted_gap = [&](std::uint8_t a) { return mf::MaxPlusCount::single(gap(a)); };
  const auto truth = [](auto...) { return true; };

  for (int affine = 0; affine < 2; affine++) {
    const auto lin_max = mf::linear<mf::MaxPlus<>>(pair, gap, gap);
    const auto aff_max = mf::affine<mf::MaxPlus<>>(pair, std::int64_t(-3), gap, gap);
    const auto lin_min = mf::linear<mf::MinPlus<>>(cost, unit, unit);
    const auto aff_min = mf::affine<mf::MinPlus<>>(cost, std::int64_t(2), unit, unit);
    const auto lin_count = mf::linear<mf::Counting>(one_count, one_count, one_count);
    const auto aff_count = mf::affine<mf::Counting>(one_count, Integer(1), one_count, one_count);
    const auto lin_mult = mf::linear<mf::MaxPlusCount>(counted_pair, counted_gap, counted_gap);
    const auto aff_mult = mf::affine<mf::MaxPlusCount>(counted_pair, mf::MaxPlusCount::single(-3), counted_gap, counted_gap);
    const auto lin_bool = mf::linear<mf::Boolean>(truth, truth, truth);
    const auto aff_bool = mf::affine<mf::Boolean>(truth, true, truth, truth);
    for (mf::Mode mode : modes) {
      bool score = true, least = true, count = true, mult = true, exists = true, trace = true, forget = true;
      int pairs = 0;
      std::size_t paths_seen = 0;
      for (std::size_t m = 0; m <= 5; m++)
        for (std::size_t n = 0; n <= 5; n++)
          for (int rep = 0; rep < 3; rep++) {
            const auto x = word(r, m, 3), y = word(r, n, 3);
            std::int64_t best = mf::MaxPlus<>::zero(), cheapest = mf::MinPlus<>::zero();
            Integer total = 0, at_best = 0;
            mf::enumerate_paths(m, n, mode, [&](const mf::EditPath& p) {
              const std::int64_t v = affine ? mf::path_weight(p, x, y, aff_max) : mf::path_weight(p, x, y, lin_max);
              const std::int64_t c = affine ? mf::path_weight(p, x, y, aff_min) : mf::path_weight(p, x, y, lin_min);
              if (v > best) { best = v; at_best = 0; }
              if (v == best) at_best += 1;
              cheapest = std::min(cheapest, c);
              total += 1;
            });
            paths_seen += static_cast<std::size_t>(total.to_int64());
            pairs++;
            const auto dp_max = affine ? mf::path_sum(x, y, aff_max, mode) : mf::path_sum(x, y, lin_max, mode);
            const auto dp_min = affine ? mf::path_sum(x, y, aff_min, mode) : mf::path_sum(x, y, lin_min, mode);
            const auto dp_count = affine ? mf::path_sum(x, y, aff_count, mode) : mf::path_sum(x, y, lin_count, mode);
            const auto dp_mult = affine ? mf::path_sum(x, y, aff_mult, mode) : mf::path_sum(x, y, lin_mult, mode);
            const bool dp_bool = affine ? mf::path_sum(x, y, aff_bool, mode) : mf::path_sum(x, y, lin_bool, mode);
            score = score && dp_max == best;
            least = least && dp_min == cheapest;
            count = count && dp_count == total;
            mult = mult && dp_mult == mf::ScoreCount{best, at_best};
            exists = exists && dp_bool;
            forget = forget && mf::ForgetCount{}(dp_mult) == dp_max;
            const auto opt = affine ? mf::optimal_path(x, y, aff_max, mode) : mf::optimal_path(x, y, lin_max, mode);
            const auto low = affine ? mf::optimal_path(x, y, aff_min, mode) : mf::optimal_path(x, y, lin_min, mode);
            const std::int64_t w = affine ? mf::path_weight(opt.path, x, y, aff_max) : mf::path_weight(opt.path, x, y, lin_max);
            const std::int64_t c = affine ? mf::path_weight(low.path, x, y, aff_min) : mf::path_weight(low.path, x, y, lin_min);
            trace = trace && opt.exists && low.exists && allowed(opt.path, m, n, mode) && allowed(low.path, m, n, mode) && w == best && opt.value == best && c == cheapest;
          }
      const std::string tag = std::string(affine ? "affine" : "linear") + ", " + name(mode) + ": ";
      const std::string detail = std::to_string(pairs) + " pairs of words of lengths 0..5 over 3 letters, " + std::to_string(paths_seen) + " paths";
      s.require(tag + "the (max, +) sum is the best score over every path", score, detail);
      s.require(tag + "the (min, +) sum is the least cost (edit distance) over every path", least);
      s.require(tag + "the Counting sum is the number of paths", count);
      s.require(tag + "the sum with multiplicities is the best score and the number of paths attaining it", mult);
      s.require(tag + "the Boolean sum says a path exists", exists);
      s.require(tag + "each traceback is a path the mode allows, of the optimal weight, in (max, +) and (min, +)", trace);
      s.require(tag + "forgetting multiplicities commutes with the sum over paths (a homomorphism)", forget);
    }
  }

  // ---- affine with a zero open is linear --------------------------------------------------------------
  {
    bool same = true;
    for (int rep = 0; rep < 200; rep++) {
      const auto x = word(r, r.below(9), 4), y = word(r, r.below(9), 4);
      for (mf::Mode mode : modes)
        same = same && mf::path_sum(x, y, mf::affine<mf::MaxPlus<>>(pair, std::int64_t(0), gap, gap), mode) == mf::path_sum(x, y, mf::linear<mf::MaxPlus<>>(pair, gap, gap), mode);
    }
    s.require("affine gaps with open = 1 (the semiring's one) are linear gaps, in every mode, on 200 pairs of length ≤ 8", same);
  }

  // ---- Delannoy numbers ---------------------------------------------------------------------------------
  {
    // The central Delannoy numbers D(n, n), n = 0 .. 10 (OEIS A001850).
    const long long central[] = {1, 3, 13, 63, 321, 1683, 8989, 48639, 265729, 1462563, 8097453};
    bool ok = true;
    for (std::size_t n = 0; n <= 10; n++) {
      const std::vector<std::uint8_t> x(n, 0), y(n, 1);
      ok = ok && mf::delannoy(n, n) == Integer(central[n]) && mf::path_sum(x, y, mf::linear<mf::Counting>(one_count, one_count, one_count), mf::Mode::global) == Integer(central[n]);
    }
    s.require("Delannoy: the global alignments of two words of length n number D(n, n) = 1, 3, 13, 63, …, 8097453 (OEIS A001850)", ok);
    bool rect = true;
    for (std::size_t m = 0; m <= 9; m++)
      for (std::size_t n = 0; n <= 9; n++) {
        const std::vector<std::uint8_t> x(m, 0), y(n, 0);
        rect = rect && mf::path_sum(x, y, mf::affine<mf::Counting>(one_count, Integer(1), one_count, one_count), mf::Mode::global) == mf::delannoy(m, n);
      }
    s.require("Delannoy: the three-state scheme counts D(m, n) for all m, n ≤ 9 (one path per word over {M, I, D})", rect);
    const std::vector<std::uint8_t> x(60, 0), y(60, 1);
    const Integer d = mf::path_sum(x, y, mf::linear<mf::Counting>(one_count, one_count, one_count), mf::Mode::global);
    s.require("Delannoy: D(60, 60), beyond 64 bits, is the closed form Σ C(60, k)² 2^k exactly", d == mf::delannoy(60, 60) && !d.fits_int64(), d.str());
  }

  // ---- a published example --------------------------------------------------------------------------------
  {
    // Durbin, Eddy, Krogh and Mitchison 1998, Biological Sequence Analysis, §2.3, Figs. 2.5, 2.6
    // and 2.8: HEAGAWGHEE and PAWHEAE, BLOSUM50, a linear gap of −8; global score 1
    // (HEAGAWGHE-E / --P-AW-HEAE), local 28 (AWGHE / AW-HE), overlap 25 (GAWGHEE / PAW-HEA).
    // BLOSUM50 entries used (NCBI's file): their letters only.
    const std::string L = "HEAGWP";
    const int B[6][6] = {// H   E   A   G   W   P
                         {10, 0, -2, -2, -3, -2},  // H
                         {0, 6, -1, -3, -3, -1},   // E
                         {-2, -1, 5, 0, -3, -1},   // A
                         {-2, -3, 0, 8, -3, -2},   // G
                         {-3, -3, -3, -3, 15, -4}, // W
                         {-2, -1, -1, -2, -4, 10}};// P
    const auto code = [&](const std::string& w) {
      std::vector<std::uint8_t> c;
      for (char ch : w) c.push_back(static_cast<std::uint8_t>(L.find(ch)));
      return c;
    };
    const auto x = code("HEAGAWGHEE"), y = code("PAWHEAE");
    const auto w = mf::linear<mf::MaxPlus<>>([&](std::uint8_t a, std::uint8_t b) { return std::int64_t(B[a][b]); }, [](std::uint8_t) { return std::int64_t(-8); },
                                             [](std::uint8_t) { return std::int64_t(-8); });
    const auto g = mf::optimal_path(x, y, w, mf::Mode::global), l = mf::optimal_path(x, y, w, mf::Mode::local), o = mf::optimal_path(x, y, w, mf::Mode::overlap);
    s.require("Durbin et al. Fig. 2.5: the global score of HEAGAWGHEE and PAWHEAE is 1", g.value == 1, std::to_string(g.value) + ", path " + g.path.ops);
    s.require("Durbin et al. Fig. 2.6: the local score is 28, AWGHE against AW-HE", l.value == 28 && l.path.ops == "MMIMM" && l.path.x_begin == 4 && l.path.y_begin == 1,
              std::to_string(l.value) + ", path " + l.path.ops);
    s.require("Durbin et al. Fig. 2.8: the overlap score is 25, GAWGHEE against PAW-HEA", o.value == 25 && o.path.ops == "MMMIMMM" && o.path.x_begin == 3 && o.path.y_begin == 0,
              std::to_string(o.value) + ", path " + o.path.ops);
  }
  return s.done();
}
