// Checks of morphe/combinatorics/eulerian.hpp: the matrix-tree theorem against Cayley's n^{n−2}
// and the complete bipartite m^{n−1} n^{m−1}, and against arborescences enumerated one by one
// on random multigraphs with loops and parallel edges; Hierholzer's circuits and trails, each
// edge once along a walk, on random Eulerian graphs and graphs with a trail, and none where
// the degrees forbid it; the BEST theorem against the enumeration of every circuit on random
// Eulerian multigraphs; the de Bruijn graphs B(σ, n − 1), whose circuits are the de Bruijn
// sequences, (σ!)^{σ^{n−1}} / σ^n of them (van Aardenne-Ehrenfest and de Bruijn 1951), the
// binary ones 2^{2^{n−1} − n} (de Bruijn 1946); long cycles through one vertex, counted after
// their transit vertices are suppressed; and the trails counted as edge sequences.

#include <morphe/check.hpp>
#include <morphe/combinatorics/eulerian.hpp>

#include <random>
#include <string>
#include <vector>

namespace mc = morphe::combinatorics;
using morphe::arithmetic::Integer;

namespace {

// Arborescences toward root by enumeration: one out-edge for each other vertex, every vertex
// reaching the root along the chosen edges.
long long brute_arborescences(const mc::Multigraph& G, std::size_t root) {
  std::vector<std::vector<std::size_t>> out(G.vertices);
  for (std::size_t e = 0; e < G.edges.size(); e++)
    if (G.edges[e][0] != G.edges[e][1]) out[G.edges[e][0]].push_back(e);
  std::vector<std::size_t> pick(G.vertices, 0);
  for (std::size_t v = 0; v < G.vertices; v++)
    if (v != root && out[v].empty()) return 0;
  long long count = 0;
  for (;;) {
    bool tree = true;
    for (std::size_t v = 0; v < G.vertices && tree; v++) {
      std::size_t x = v, steps = 0;
      while (x != root && steps <= G.vertices) { x = G.edges[out[x][pick[x]]][1]; steps++; }
      tree = x == root;
    }
    count += tree;
    std::size_t v = 0;
    for (; v < G.vertices; v++) {
      if (v == root) continue;
      if (++pick[v] < out[v].size()) break;
      pick[v] = 0;
    }
    if (v == G.vertices) return count;
  }
}

// A walk that uses every edge once, from edge sequence c; closed if asked.
bool is_eulerian_walk(const mc::Multigraph& G, const std::vector<std::size_t>& c, bool closed) {
  if (c.size() != G.edges.size()) return false;
  std::vector<char> seen(G.edges.size(), 0);
  for (std::size_t i = 0; i < c.size(); i++) {
    if (c[i] >= G.edges.size() || seen[c[i]]) return false;
    seen[c[i]] = 1;
    if (i + 1 < c.size() && G.edges[c[i]][1] != G.edges[c[i + 1]][0]) return false;
  }
  return !closed || c.empty() || G.edges[c.back()][1] == G.edges[c.front()][0];
}

// A random connected Eulerian multigraph: the union of random closed walks through vertex 0.
mc::Multigraph random_eulerian(std::mt19937_64& g, std::size_t n, std::size_t walks, std::size_t length) {
  mc::Multigraph G;
  G.vertices = n;
  for (std::size_t w = 0; w < walks; w++) {
    std::size_t v = 0;
    for (std::size_t i = 0; i + 1 < length; i++) {
      const std::size_t u = g() % n;
      G.edges.push_back({v, u});
      v = u;
    }
    G.edges.push_back({v, 0});
  }
  // Shuffle the edge indices, so that edge 0 is not special.
  std::shuffle(G.edges.begin(), G.edges.end(), g);
  return G;
}

// The de Bruijn graph B(σ, n): vertices the words of length n over σ letters (base σ), an
// edge x → y for each word of length n + 1 whose prefix is x and suffix y.
mc::Multigraph de_bruijn(std::size_t sigma, std::size_t n) {
  std::size_t V = 1;
  for (std::size_t i = 0; i < n; i++) V *= sigma;
  mc::Multigraph G;
  G.vertices = V;
  for (std::size_t x = 0; x < V; x++)
    for (std::size_t a = 0; a < sigma; a++) G.edges.push_back({x, (x * sigma + a) % V});
  return G;
}

}  // namespace

int main() {
  morphe::Suite s;
  std::mt19937_64 g(1951);

  // Cayley: the complete symmetric digraph on n vertices has n^{n−2} arborescences toward a root.
  {
    bool ok = true;
    for (std::size_t n = 1; n <= 12; n++) {
      mc::Multigraph K;
      K.vertices = n;
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t v = 0; v < n; v++)
          if (u != v) K.edges.push_back({u, v});
      ok = ok && mc::count_arborescences(K, n - 1) == (n == 1 ? Integer(1) : morphe::arithmetic::pow(Integer(static_cast<long long>(n)), n - 2));
    }
    s.require("matrix-tree: the complete digraph K_n has n^(n−2) arborescences toward each root, n = 1 to 12 (Cayley 1889)", ok);
    ok = true;
    for (std::size_t a = 1; a <= 5; a++)
      for (std::size_t b = 1; b <= 5; b++) {
        mc::Multigraph K;
        K.vertices = a + b;
        for (std::size_t u = 0; u < a; u++)
          for (std::size_t v = a; v < a + b; v++) {
            K.edges.push_back({u, v});
            K.edges.push_back({v, u});
          }
        const Integer want = morphe::arithmetic::pow(Integer(static_cast<long long>(a)), b - 1) * morphe::arithmetic::pow(Integer(static_cast<long long>(b)), a - 1);
        ok = ok && mc::count_arborescences(K, 0) == want;
      }
    s.require("matrix-tree: the complete bipartite K_{a,b} has a^(b−1) b^(a−1) spanning trees, a, b = 1 to 5 (Scoins 1962)", ok);
  }

  // Arborescences one by one, on random directed multigraphs with loops and parallel edges.
  {
    int bad = 0;
    for (int t = 0; t < 300; t++) {
      mc::Multigraph G;
      G.vertices = 1 + g() % 6;
      const std::size_t m = g() % 16;
      for (std::size_t e = 0; e < m; e++) G.edges.push_back({g() % G.vertices, g() % G.vertices});
      const std::size_t root = g() % G.vertices;
      if (mc::count_arborescences(G, root) != Integer(brute_arborescences(G, root))) bad++;
    }
    s.equal("matrix-tree: det L_w equals the arborescences enumerated on 300 random multigraphs with loops (failures)", bad, 0);
  }

  // Hierholzer on random Eulerian graphs and on graphs with a trail.
  {
    int bad = 0;
    for (int t = 0; t < 300; t++) {
      mc::Multigraph G = random_eulerian(g, 1 + g() % 8, 1 + g() % 4, 1 + g() % 12);
      const std::size_t first = g() % G.edges.size();
      const auto c = mc::eulerian_circuit(G, first);
      if (!c || !is_eulerian_walk(G, *c, true) || c->front() != first) bad++;
      // Remove one edge u → v: a trail from v to u remains (or a circuit, for a loop).
      const auto removed = G.edges.back();
      G.edges.pop_back();
      const auto st = mc::eulerian_status(G);
      const auto tr = mc::eulerian_trail(G);
      if (!tr || !is_eulerian_walk(G, *tr, false)) bad++;
      if (removed[0] != removed[1] && !G.edges.empty() &&
          (st.kind != mc::EulerianKind::trail || st.start != removed[1] || st.end != removed[0] || G.edges[tr->front()][0] != removed[1])) bad++;
    }
    s.equal("Hierholzer: circuits from a given edge and trails from the source, every edge once along a walk, on 300 random graphs (failures)", bad, 0);
    mc::Multigraph P{3, {{0, 2}, {1, 2}, {2, 0}, {2, 1}, {0, 2}, {1, 2}}};
    s.require("Hierholzer: no circuit or trail when two vertices have out − in = 1", !mc::eulerian_circuit(P) && !mc::eulerian_trail(P) && mc::count_eulerian_circuits(P).is_zero());
    mc::Multigraph D{4, {{0, 1}, {1, 0}, {2, 3}, {3, 2}}};
    s.require("Hierholzer: no circuit when the edges lie in two components", !mc::eulerian_circuit(D) && mc::count_eulerian_trails(D).is_zero());
  }

  // BEST against the enumeration of every circuit.
  {
    int bad = 0, nontrivial = 0, tried = 0;
    long long largest = 0;
    while (tried < 200) {
      const mc::Multigraph G = random_eulerian(g, 1 + g() % 5, 1 + g() % 3, 1 + g() % 5);
      // Enumeration is exponential: graphs with more than 10^5 circuits are passed over.
      if (mc::count_eulerian_circuits(G) > Integer(100000)) continue;
      tried++;
      std::size_t lexicographic = 1;
      std::vector<std::size_t> previous;
      const std::size_t n = mc::for_each_eulerian_circuit(G, [&](const std::vector<std::size_t>& c) {
        if (!is_eulerian_walk(G, c, true) || c.front() != 0 || (!previous.empty() && !(previous < c))) lexicographic = 0;
        previous = c;
        return true;
      });
      if (mc::count_eulerian_circuits(G) != Integer(static_cast<unsigned long long>(n)) || !lexicographic) bad++;
      if (n > 1) nontrivial++;
      largest = std::max(largest, static_cast<long long>(n));
    }
    s.require("BEST: t_w ∏ (out(v) − 1)! equals the circuits enumerated on 200 random Eulerian multigraphs",
              bad == 0, std::to_string(nontrivial) + " with several circuits, at most " + std::to_string(largest));
    // The complete symmetric digraph K_n*: ec = n^{n−2} ((n − 2)!)^n, enumerated for n ≤ 4.
    bool ok = true;
    for (std::size_t n = 2; n <= 4; n++) {
      mc::Multigraph K;
      K.vertices = n;
      for (std::size_t u = 0; u < n; u++)
        for (std::size_t v = 0; v < n; v++)
          if (u != v) K.edges.push_back({u, v});
      const Integer want = morphe::arithmetic::pow(Integer(static_cast<long long>(n)), n - 2) *
                           morphe::arithmetic::pow(morphe::arithmetic::factorial(n - 2), n);
      const std::size_t e = mc::for_each_eulerian_circuit(K, [](const std::vector<std::size_t>&) { return true; });
      ok = ok && mc::count_eulerian_circuits(K) == want && Integer(static_cast<unsigned long long>(e)) == want;
    }
    s.require("BEST: the complete digraph K_n* has n^(n−2) ((n−2)!)^n circuits, enumerated for n = 2 to 4", ok);
  }

  // Large graphs: three cycles of 50000, 30000 and 20000 vertices through vertex 0 have
  // (3 − 1)! = 2 circuits; one cycle of 100000 vertices has one. The transit vertices are
  // suppressed, so the determinant is of order 0.
  {
    mc::Multigraph F;
    F.vertices = 1;
    for (std::size_t len : {50000u, 30000u, 20000u}) {
      std::size_t prev = 0;
      for (std::size_t i = 1; i < len; i++) {
        F.edges.push_back({prev, F.vertices});
        prev = F.vertices++;
      }
      F.edges.push_back({prev, 0});
    }
    mc::Multigraph C;
    C.vertices = 100000;
    for (std::size_t v = 0; v < C.vertices; v++) C.edges.push_back({v, (v + 1) % C.vertices});
    const auto H = mc::suppress_transit_vertices(F);
    s.require("BEST: three long cycles through one vertex have 2 circuits, one long cycle 1, after suppressing transit vertices",
              mc::count_eulerian_circuits(F) == Integer(2) && mc::count_eulerian_circuits(C) == Integer(1) && H.vertices == 1 && H.edges.size() == 3,
              std::to_string(F.vertices) + " vertices reduced to " + std::to_string(H.vertices));
  }

  // De Bruijn sequences: the circuits of B(σ, n − 1).
  {
    bool ok = true;
    std::string detail;
    for (std::size_t n = 1; n <= 8; n++) {
      // de Bruijn 1946: 2^{2^{n−1} − n} binary sequences of order n.
      const Integer want = morphe::arithmetic::pow(Integer(2), (1ull << (n - 1)) - n);
      ok = ok && mc::count_eulerian_circuits(de_bruijn(2, n - 1)) == want;
      if (n == 8) detail = "B(2, 7): " + want.str();
    }
    s.require("BEST: B(2, n − 1) has 2^(2^(n−1) − n) Eulerian circuits, the binary de Bruijn sequences of order n = 1 to 8", ok, detail);
    ok = true;
    for (std::size_t sigma = 2; sigma <= 4; sigma++)
      for (std::size_t n = 1; n <= (sigma == 4 ? 4u : 5u); n++) {
        std::size_t power = 1;  // σ^{n−1}
        for (std::size_t i = 1; i < n; i++) power *= sigma;
        const Integer want = morphe::arithmetic::exact_div(morphe::arithmetic::pow(morphe::arithmetic::factorial(sigma), power),
                                                           morphe::arithmetic::pow(Integer(static_cast<long long>(sigma)), n));
        if (mc::count_eulerian_circuits(de_bruijn(sigma, n - 1)) != want) ok = false;
        if (sigma == 4 && n == 4) detail = "σ = 4, n = 4: " + want.str();
      }
    s.require("BEST: B(σ, n − 1) has (σ!)^(σ^(n−1)) / σ^n circuits, σ = 2, 3 (n ≤ 5) and 4 (n ≤ 4)", ok, detail);
    // Enumerated: 1, 1, 2, 16, 2048 binary sequences of orders 1 to 5.
    std::vector<std::size_t> counted;
    for (std::size_t n = 1; n <= 5; n++) counted.push_back(mc::for_each_eulerian_circuit(de_bruijn(2, n - 1), [](const std::vector<std::size_t>&) { return true; }));
    s.equal("BEST: the binary de Bruijn sequences of orders 1 to 5, enumerated", counted, std::vector<std::size_t>{1, 1, 2, 16, 2048});
    counted.clear();
    counted.push_back(mc::for_each_eulerian_circuit(de_bruijn(3, 1), [](const std::vector<std::size_t>&) { return true; }));
    counted.push_back(mc::for_each_eulerian_circuit(de_bruijn(4, 1), [](const std::vector<std::size_t>&) { return true; }));
    s.equal("BEST: the de Bruijn sequences of order 2 over 3 and 4 letters, enumerated: 24 and 20736", counted, std::vector<std::size_t>{24, 20736});
  }

  // Trails as edge sequences: m ec(G) on an Eulerian graph; on a graph with a trail, the trails
  // enumerated as the circuits of G + (t → s) cut at the added edge.
  {
    int bad = 0;
    for (int t = 0; t < 100; t++) {
      mc::Multigraph G = random_eulerian(g, 1 + g() % 4, 1 + g() % 3, 1 + g() % 4);
      if (mc::count_eulerian_trails(G) != mc::count_eulerian_circuits(G) * Integer(static_cast<unsigned long long>(G.edges.size()))) bad++;
      // Make the trail graph: drop the last edge u → v (u ≠ v).
      const auto e = G.edges.back();
      if (e[0] == e[1]) continue;
      G.edges.pop_back();
      // Count trails by brute force: circuits of G + (u → v), each from the added edge.
      mc::Multigraph H = G;
      H.edges.insert(H.edges.begin(), e);
      if (mc::count_eulerian_circuits(H) > Integer(100000)) continue;
      std::size_t trails = mc::for_each_eulerian_circuit(H, [](const std::vector<std::size_t>&) { return true; });
      if (mc::count_eulerian_trails(G) != Integer(static_cast<unsigned long long>(trails))) bad++;
    }
    s.equal("BEST: Eulerian trails counted as edge sequences against enumeration on 100 graphs (failures)", bad, 0);
  }

  return s.done();
}
