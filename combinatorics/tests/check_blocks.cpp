// Checks of blocks, cut vertices and bridges (blocks.hpp).
//
// - Small graphs whose answers are known: a path, a triangle, two triangles sharing a vertex, a bridge between two cycles, parallel
//   edges (never a bridge), a loop (in no block), an isolated vertex.
// - Random multigraphs: a vertex is a cut vertex exactly when removing it leaves more components among the others; an edge is a
//   bridge exactly when removing it leaves more components; every non-loop edge lies in one block, and two blocks share at most one
//   vertex, a cut vertex.
// - A path of 200,000 vertices: every inner vertex is a cut vertex and every edge a bridge, the search being iterative.

#include <morphe/check.hpp>
#include <morphe/combinatorics/blocks.hpp>

#include <random>
#include <set>
#include <string>
#include <vector>

namespace mc = morphe::combinatorics;

namespace {

// Components of the graph with some vertices or edges removed, counted over the vertices that remain.
std::size_t components(const mc::Multigraph& G, std::size_t without_vertex, std::size_t without_edge) {
  mc::Multigraph H;
  std::vector<std::size_t> map(G.vertices, SIZE_MAX);
  for (std::size_t v = 0; v < G.vertices; v++)
    if (v != without_vertex) map[v] = H.vertices++;
  for (std::size_t e = 0; e < G.edges.size(); e++) {
    if (e == without_edge) continue;
    const auto [a, b] = G.edges[e];
    if (a == without_vertex || b == without_vertex) continue;
    H.edges.push_back({map[a], map[b]});
  }
  return mc::spanning_forest(H).components();
}

mc::Multigraph graph(std::size_t n, std::vector<std::array<std::size_t, 2>> edges) { return {n, std::move(edges)}; }

}  // namespace

int main() {
  morphe::Suite s;
  {
    const auto b = mc::blocks(graph(3, {{0, 1}, {1, 2}}));
    s.require("a path: its middle vertex is a cut vertex, both edges are bridges, each its own block",
              b.cut == std::vector<char>{0, 1, 0} && b.bridges() == 2 && b.blocks.size() == 2);
  }
  {
    const auto b = mc::blocks(graph(3, {{0, 1}, {1, 2}, {2, 0}}));
    s.require("a triangle: no cut vertex, no bridge, one block", b.cut_vertices() == 0 && b.bridges() == 0 && b.blocks.size() == 1 && b.blocks[0].size() == 3);
  }
  {
    const auto b = mc::blocks(graph(5, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 2}}));
    s.require("two triangles sharing a vertex: it is the one cut vertex, between two blocks", b.cut == std::vector<char>{0, 0, 1, 0, 0} && b.blocks.size() == 2 && b.bridges() == 0);
  }
  {
    const auto b = mc::blocks(graph(6, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 5}, {5, 3}}));
    s.require("two cycles joined by an edge: the edge is the bridge, its ends the cut vertices, three blocks",
              b.bridges() == 1 && b.bridge[3] && b.cut == std::vector<char>{0, 0, 1, 1, 0, 0} && b.blocks.size() == 3);
  }
  {
    const auto b = mc::blocks(graph(3, {{0, 1}, {1, 0}, {1, 2}, {2, 2}}));
    s.require("parallel edges are never a bridge; the single edge is; a loop lies in no block",
              !b.bridge[0] && !b.bridge[1] && b.bridge[2] && b.cut == std::vector<char>{0, 1, 0} && b.block_of[3] == mc::Blocks::none && b.block_of[0] == b.block_of[1]);
  }
  {
    const auto b = mc::blocks(graph(4, {{0, 1}}));
    s.require("isolated vertices are no cut vertices", b.cut_vertices() == 0 && b.bridges() == 1);
  }
  {
    std::mt19937_64 rng(20261119);
    std::size_t graphs = 0, bad = 0;
    std::string why;
    for (int k = 0; k < 400; k++) {
      std::uniform_int_distribution<std::size_t> nv(1, 12);
      const std::size_t n = nv(rng);
      std::uniform_int_distribution<std::size_t> pick(0, n - 1), many(0, 2 * n);
      mc::Multigraph G{n, {}};
      const std::size_t m = many(rng);
      for (std::size_t e = 0; e < m; e++) G.edges.push_back({pick(rng), pick(rng)});
      const auto b = mc::blocks(G);
      const std::size_t base = mc::spanning_forest(G).components();
      for (std::size_t v = 0; v < n; v++) {
        // removing an isolated vertex removes a component; any other vertex is a cut vertex when the count rises above the graph's
        const bool isolated = [&] {
          for (const auto& e : G.edges)
            if ((e[0] == v || e[1] == v) && e[0] != e[1]) return false;
          return true;
        }();
        const std::size_t after = components(G, v, SIZE_MAX);
        const bool is_cut = !isolated && after > base;
        if (static_cast<bool>(b.cut[v]) != is_cut) bad++, why = "vertex " + std::to_string(v) + " of a graph with " + std::to_string(n) + " vertices";
      }
      for (std::size_t e = 0; e < m; e++) {
        const bool is_bridge = G.edges[e][0] != G.edges[e][1] && components(G, SIZE_MAX, e) > base;
        if (static_cast<bool>(b.bridge[e]) != is_bridge) bad++, why = "edge " + std::to_string(e);
        const bool loop = G.edges[e][0] == G.edges[e][1];
        if (loop != (b.block_of[e] == mc::Blocks::none)) bad++, why = "edge " + std::to_string(e) + " in no block or a loop in one";
      }
      // two blocks share at most one vertex, and it is a cut vertex
      for (std::size_t i = 0; i < b.blocks.size(); i++)
        for (std::size_t j = i + 1; j < b.blocks.size(); j++) {
          std::set<std::size_t> vi, vj, shared;
          for (auto e : b.blocks[i]) vi.insert(G.edges[e][0]), vi.insert(G.edges[e][1]);
          for (auto e : b.blocks[j]) vj.insert(G.edges[e][0]), vj.insert(G.edges[e][1]);
          for (auto v : vi)
            if (vj.count(v)) shared.insert(v);
          if (shared.size() > 1 || (shared.size() == 1 && !b.cut[*shared.begin()])) bad++, why = "two blocks share more than a cut vertex";
        }
      graphs++;
    }
    s.require("on 400 random multigraphs, cut vertices and bridges are exactly what removal shows, and blocks meet only at cut vertices", bad == 0,
              std::to_string(graphs) + " graphs" + (bad ? "; " + why : ""));
  }
  {
    const std::size_t n = 200000;
    mc::Multigraph P{n, {}};
    for (std::size_t v = 0; v + 1 < n; v++) P.edges.push_back({v, v + 1});
    const auto b = mc::blocks(P);
    s.require("a path of 200,000 vertices: every inner vertex is a cut vertex and every edge a bridge, without exhausting the stack",
              b.cut_vertices() == n - 2 && b.bridges() == n - 1 && !b.cut[0] && !b.cut[n - 1]);
  }
  return s.done();
}
