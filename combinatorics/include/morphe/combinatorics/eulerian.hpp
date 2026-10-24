// Eulerian circuits and trails of directed multigraphs: their existence, Hierholzer's
// algorithm, the enumeration of all of them on small graphs, the number of spanning
// arborescences by the matrix-tree theorem and the number of Eulerian circuits by the BEST
// theorem, exactly.
//
// A Multigraph (graph.hpp) is read here as directed: edge e runs from edges[e][0] (its tail)
// to edges[e][1] (its head). Loops and parallel edges are allowed, and parallel edges are
// distinct edges. An Eulerian circuit is a closed walk that uses every edge exactly once,
// taken as a cyclic sequence of edges: its m rotations are one circuit. An Eulerian trail is a
// walk, with a first edge, that uses every edge exactly once.
//
// - Existence (Euler 1736; Hierholzer 1873, "Über die Möglichkeit, einen Linienzug ohne
//   Wiederholung und ohne Unterbrechung zu umfahren", Math. Ann. 6, 30–32; for directed graphs
//   Good 1946, J. London Math. Soc. 21, 167–169): a directed multigraph has an Eulerian circuit
//   if and only if in(v) = out(v) at every vertex and its edges lie in one weakly connected
//   component; it has an Eulerian trail from s ≠ t if and only if out(s) = in(s) + 1,
//   in(t) = out(t) + 1, every other vertex is balanced, and its edges are connected.
// - Hierholzer's algorithm, in its iterative form: walk unused edges from the start until
//   stuck, then back up along the walk, emitting edges, and walk again from the first vertex
//   that still has unused edges. O(n + m). At every vertex the unused edge of least index is
//   taken first, so that the result is a function of the graph.
// - The matrix-tree theorem for directed graphs (Tutte 1948, "The dissection of equilateral
//   triangles into equilateral triangles", Proc. Cambridge Philos. Soc. 44, 463–482;
//   Borchardt 1860 and Kirchhoff 1847 for undirected graphs; Stanley 1999, Enumerative
//   Combinatorics 2, §5.6): the number of spanning arborescences oriented toward a root w,
//   in which every vertex has exactly one directed path to w, is det L_w, the determinant of
//   the out-degree Laplacian L = D_out − A (A_uv the number of edges u → v; loops cancel)
//   with row and column w deleted. Computed exactly by Bareiss's elimination
//   (morphe/algebra/determinant.hpp), O(n³) operations on integers of O(n log n) bits.
// - The BEST theorem (van Aardenne-Ehrenfest and de Bruijn 1951, "Circuits and trees in
//   oriented linear graphs", Simon Stevin 28, 203–217; Smith and Tutte 1941, "On unicursal
//   paths in a network of degree 4", Amer. Math. Monthly 48, 233–237; Stanley 1999, §5.6):
//   for a connected Eulerian directed multigraph, the number of Eulerian circuits is
//
//     ec(G) = t_w(G) ∏_v (out(v) − 1)!,
//
//   the same for every vertex w, the product over the vertices that have edges. A vertex with
//   one edge in and one out contributes 0! = 1 and no choice, so it is suppressed first (its
//   two edges joined into one), and the determinant is taken over the vertices where a circuit
//   chooses: a graph of a million vertices with a hundred branchings costs a determinant of
//   order a hundred.
// - Eulerian trails, as edge sequences: m · ec(G) for an Eulerian G (each circuit cut at each
//   of its m edges), and ec(G + (t → s)) for a trail from s to t (each circuit of G with the
//   edge t → s added passes it once, and cut there is a trail).
// - Enumeration: depth-first search over the unused out-edges, each circuit written from edge
//   0 so that it appears once. It visits ec(G) circuits with O(m) work for each extension,
//   exponential in general; for small graphs, and for the checks of the counts.
//
// A graph with no edges has one (empty) circuit and one (empty) trail. Vertices with no edges
// take no part in circuits; they are left out of the arborescences of the BEST theorem.

#pragma once

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

#include <morphe/algebra/determinant.hpp>
#include <morphe/arithmetic/integer.hpp>

#include "graph.hpp"

namespace morphe::combinatorics {

using arithmetic::Integer;

struct Degrees {
  std::vector<std::size_t> in, out;
};

inline Degrees degrees(const Multigraph& G) {
  Degrees d{std::vector<std::size_t>(G.vertices, 0), std::vector<std::size_t>(G.vertices, 0)};
  for (const auto& [tail, head] : G.edges) {
    if (tail >= G.vertices || head >= G.vertices) throw std::out_of_range("morphe::combinatorics::degrees: an edge leaves the vertices");
    d.out[tail]++;
    d.in[head]++;
  }
  return d;
}

// Whether the edges lie in one weakly connected component (true when there are none).
inline bool edges_connected(const Multigraph& G) {
  if (G.edges.empty()) return true;
  const SpanningForest F = spanning_forest(G, G.edges[0][0]);
  const std::size_t c = F.component[G.edges[0][0]];
  for (const auto& e : G.edges)
    if (F.component[e[0]] != c) return false;
  return true;
}

enum class EulerianKind { none, circuit, trail };

struct EulerianStatus {
  EulerianKind kind = EulerianKind::none;
  std::size_t start = 0, end = 0;  // of a trail; for a circuit, the tail of edge 0, twice
};

inline EulerianStatus eulerian_status(const Multigraph& G) {
  EulerianStatus st;
  if (G.edges.empty()) { st.kind = EulerianKind::circuit; return st; }
  if (!edges_connected(G)) return st;
  const Degrees d = degrees(G);
  std::size_t sources = 0, sinks = 0;
  for (std::size_t v = 0; v < G.vertices; v++) {
    if (d.out[v] == d.in[v]) continue;
    if (d.out[v] == d.in[v] + 1) { sources++; st.start = v; continue; }
    if (d.in[v] == d.out[v] + 1) { sinks++; st.end = v; continue; }
    return st;
  }
  if (sources == 0 && sinks == 0) {
    st.kind = EulerianKind::circuit;
    st.start = st.end = G.edges[0][0];
  } else if (sources == 1 && sinks == 1) {
    st.kind = EulerianKind::trail;
  }
  return st;
}

namespace detail {

// Out-edges of each vertex, in increasing edge index.
inline std::vector<std::vector<std::size_t>> out_edges(const Multigraph& G) {
  std::vector<std::vector<std::size_t>> out(G.vertices);
  for (std::size_t e = 0; e < G.edges.size(); e++) out[G.edges[e][0]].push_back(e);
  return out;
}

// Hierholzer from vertex `start`: the edges in walk order; the caller has checked existence.
inline std::vector<std::size_t> hierholzer(const Multigraph& G, std::size_t start) {
  const auto out = out_edges(G);
  std::vector<std::size_t> next(G.vertices, 0);  // the next unused out-edge of each vertex
  std::vector<std::size_t> walk;                  // edges, from the start
  std::vector<std::size_t> stack;                 // edges of the current walk
  std::vector<std::size_t> at;                    // the vertex reached by each edge of stack
  walk.reserve(G.edges.size());
  std::size_t v = start;
  for (;;) {
    if (next[v] < out[v].size()) {
      const std::size_t e = out[v][next[v]++];
      stack.push_back(e);
      v = G.edges[e][1];
      at.push_back(v);
      continue;
    }
    if (stack.empty()) break;
    walk.push_back(stack.back());
    stack.pop_back();
    at.pop_back();
    v = at.empty() ? start : at.back();
  }
  std::reverse(walk.begin(), walk.end());
  return walk;
}

}  // namespace detail

// An Eulerian circuit as a sequence of edges beginning with edge `first`; none if G has no
// Eulerian circuit. Empty for a graph without edges.
inline std::optional<std::vector<std::size_t>> eulerian_circuit(const Multigraph& G, std::size_t first = 0) {
  const EulerianStatus st = eulerian_status(G);
  if (st.kind != EulerianKind::circuit) return std::nullopt;
  if (G.edges.empty()) return std::vector<std::size_t>{};
  if (first >= G.edges.size()) throw std::out_of_range("morphe::combinatorics::eulerian_circuit: no such edge");
  std::vector<std::size_t> c = detail::hierholzer(G, G.edges[first][0]);
  std::rotate(c.begin(), std::find(c.begin(), c.end(), first), c.end());
  return c;
}

// An Eulerian trail: from the unique source to the unique sink when G has one, else an
// Eulerian circuit from edge 0; none if G has neither.
inline std::optional<std::vector<std::size_t>> eulerian_trail(const Multigraph& G) {
  const EulerianStatus st = eulerian_status(G);
  if (st.kind == EulerianKind::none) return std::nullopt;
  if (G.edges.empty()) return std::vector<std::size_t>{};
  return detail::hierholzer(G, st.start);
}

// Calls f(circuit) for every Eulerian circuit, written from edge 0, in lexicographic order of
// its edge sequence; f returns false to stop. Returns the number of circuits visited.
template <class F>
std::size_t for_each_eulerian_circuit(const Multigraph& G, F&& f) {
  if (eulerian_status(G).kind != EulerianKind::circuit) return 0;
  const std::size_t m = G.edges.size();
  std::vector<std::size_t> circuit;
  if (m == 0) { f(circuit); return 1; }
  const auto out = detail::out_edges(G);
  std::vector<char> used(m, 0);
  std::vector<std::size_t> choice;  // per depth, the index into out[] of the edge taken
  circuit.push_back(0);
  used[0] = 1;
  choice.push_back(0);
  std::size_t visited = 0;
  // Depth-first: at depth d the circuit has d + 1 edges; extend from the head of the last.
  std::size_t from = 0;  // the next index into out[v] to try at this depth
  while (!circuit.empty()) {
    if (circuit.size() == m) {
      visited++;
      if (!f(static_cast<const std::vector<std::size_t>&>(circuit))) return visited;
    } else {
      const std::size_t v = G.edges[circuit.back()][1];
      std::size_t i = from;
      while (i < out[v].size() && used[out[v][i]]) i++;
      if (i < out[v].size()) {
        const std::size_t e = out[v][i];
        used[e] = 1;
        circuit.push_back(e);
        choice.push_back(i);
        from = 0;
        continue;
      }
    }
    // Back up: undo the last edge and try the next choice at its depth.
    if (circuit.size() == 1) break;
    used[circuit.back()] = 0;
    circuit.pop_back();
    from = choice.back() + 1;
    choice.pop_back();
  }
  return visited;
}

// The out-degree Laplacian L = D_out − A, n × n.
inline algebra::Matrix laplacian(const Multigraph& G) {
  algebra::Matrix L(G.vertices, G.vertices);
  for (const auto& [tail, head] : G.edges) {
    if (tail >= G.vertices || head >= G.vertices) throw std::out_of_range("morphe::combinatorics::laplacian: an edge leaves the vertices");
    if (tail == head) continue;
    L(tail, tail) += 1;
    L(tail, head) -= 1;
  }
  return L;
}

// The number of spanning arborescences oriented toward `root` (every vertex has one directed
// path to it), by the matrix-tree theorem.
inline Integer count_arborescences(const Multigraph& G, std::size_t root) {
  if (root >= G.vertices) throw std::out_of_range("morphe::combinatorics::count_arborescences: no such root");
  return algebra::minor(laplacian(G), root, root);
}

// The graph with its transit vertices suppressed: a vertex with one edge in and one edge out
// (not a loop) is passed straight through by every walk, so each maximal chain of edges through
// transit vertices becomes one edge, and vertices without edges are dropped. On a connected
// Eulerian graph the circuits of G and of the result correspond one to one, as cyclic
// sequences, and every vertex kept has the out-degree it had; a graph that is one cycle becomes
// one vertex with one loop. O(n + m).
inline Multigraph suppress_transit_vertices(const Multigraph& G) {
  const Degrees d = degrees(G);
  std::vector<std::size_t> out_edge(G.vertices, SIZE_MAX), id(G.vertices, SIZE_MAX);
  for (std::size_t e = 0; e < G.edges.size(); e++) out_edge[G.edges[e][0]] = e;
  const auto transit = [&](std::size_t v) { return d.in[v] == 1 && d.out[v] == 1 && G.edges[out_edge[v]][1] != v; };
  Multigraph H;
  for (std::size_t v = 0; v < G.vertices; v++)
    if ((d.in[v] > 0 || d.out[v] > 0) && !transit(v)) id[v] = H.vertices++;
  if (H.vertices == 0) {
    if (!G.edges.empty()) H = Multigraph{1, {{0, 0}}};
    return H;
  }
  for (const auto& [tail, head] : G.edges) {
    if (id[tail] == SIZE_MAX) continue;
    std::size_t v = head;
    for (std::size_t steps = 0; id[v] == SIZE_MAX && steps <= G.edges.size(); steps++) v = G.edges[out_edge[v]][1];
    H.edges.push_back({id[tail], id[v]});
  }
  return H;
}

// The number of Eulerian circuits (cyclic sequences of edges) by the BEST theorem; 0 when G
// has none. The transit vertices are suppressed first, so that the determinant is over the
// vertices where a circuit has a choice: O(n + m + b³) for b of them.
inline Integer count_eulerian_circuits(const Multigraph& G) {
  const EulerianStatus st = eulerian_status(G);
  if (st.kind != EulerianKind::circuit) return 0;
  if (G.edges.empty()) return 1;
  const Multigraph H = suppress_transit_vertices(G);
  const Degrees d = degrees(H);
  Integer ec = count_arborescences(H, 0);
  for (std::size_t v = 0; v < H.vertices; v++)
    if (d.out[v] > 1) ec *= arithmetic::factorial(d.out[v] - 1);
  return ec;
}

// The number of Eulerian trails, as edge sequences: m ec(G) when G is Eulerian, ec(G + (t → s))
// when G has a trail from s to t, and 0 otherwise.
inline Integer count_eulerian_trails(const Multigraph& G) {
  const EulerianStatus st = eulerian_status(G);
  if (st.kind == EulerianKind::none) return 0;
  if (G.edges.empty()) return 1;
  if (st.kind == EulerianKind::circuit) return count_eulerian_circuits(G) * Integer(static_cast<unsigned long long>(G.edges.size()));
  Multigraph H = G;
  H.edges.push_back({st.end, st.start});
  return count_eulerian_circuits(H);
}

}  // namespace morphe::combinatorics
