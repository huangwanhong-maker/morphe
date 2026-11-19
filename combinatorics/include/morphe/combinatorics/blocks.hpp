// Blocks, cut vertices and bridges of a multigraph (decision 0020).
//
// Read as undirected, a connected multigraph falls apart at a cut vertex (removing it leaves more components among the vertices
// that remain) or at a bridge (an edge whose removal leaves more components). Its blocks are its maximal subgraphs with no cut
// vertex of their own; every edge lies in exactly one block, two blocks share at most one vertex, which is a cut vertex, and the
// blocks and cut vertices form a forest (the block–cut tree). A bridge is a block of one edge whose ends are joined by nothing else.
//
// All three are found in one depth-first search, O(n + m), by Hopcroft and Tarjan's lowpoints (Hopcroft and Tarjan 1973,
// "Algorithm 447: efficient algorithms for graph manipulation", Communications of the ACM 16(6), 372–378; Tarjan 1972, "Depth-first
// search and linear graph algorithms", SIAM Journal on Computing 1(2), 146–160): low(v) is the earliest discovery time reachable
// from v's subtree by tree edges down and one other edge up. A child w of v with low(w) ≥ disc(v) cannot reach above v without it,
// so v is a cut vertex (a root only when it has two children or more) and the edges stacked since w form a block; low(w) > disc(v)
// makes the tree edge a bridge.
//
// The search is iterative, so that a long path cannot exhaust the stack. A loop lies in no block that disconnects anything and is
// left out of the blocks (block_of says so); parallel edges are told apart by their identities, so that two edges joining the same
// pair are never a bridge.

#pragma once

#include <morphe/combinatorics/graph.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace morphe::combinatorics {

struct Blocks {
  static constexpr std::size_t none = SIZE_MAX;
  std::vector<char> cut;                          // per vertex: a cut vertex
  std::vector<char> bridge;                       // per edge: a bridge
  std::vector<std::vector<std::size_t>> blocks;   // each block's edges
  std::vector<std::size_t> block_of;              // per edge: its block, none for a loop
  std::size_t cut_vertices() const { return static_cast<std::size_t>(std::count(cut.begin(), cut.end(), 1)); }
  std::size_t bridges() const { return static_cast<std::size_t>(std::count(bridge.begin(), bridge.end(), 1)); }
};

inline Blocks blocks(const Multigraph& G) {
  const std::size_t n = G.vertices, m = G.edges.size();
  const auto inc = G.incidence();
  Blocks out;
  out.cut.assign(n, 0);
  out.bridge.assign(m, 0);
  out.block_of.assign(m, Blocks::none);
  constexpr std::size_t unseen = SIZE_MAX;
  std::vector<std::size_t> disc(n, unseen), low(n, 0);
  std::vector<char> edge_stacked(m, 0);
  std::vector<std::size_t> edge_stack;
  struct Frame {
    std::size_t v, parent_edge, next;  // the vertex, the tree edge into it (none at a root), the next incidence to look at
  };
  std::size_t time = 0;
  for (std::size_t root = 0; root < n; root++) {
    if (disc[root] != unseen) continue;
    std::vector<Frame> stack{{root, Blocks::none, 0}};
    disc[root] = low[root] = time++;
    std::size_t root_children = 0;
    while (!stack.empty()) {
      Frame& f = stack.back();
      if (f.next < inc[f.v].size()) {
        const auto [e, end] = inc[f.v][f.next++];
        const std::size_t w = G.edges[e][end > 0 ? 1 : 0];
        if (w == f.v) continue;               // a loop: in no block
        if (e == f.parent_edge) continue;     // the tree edge we came by (a parallel edge is another identity)
        if (disc[w] == unseen) {
          edge_stack.push_back(e), edge_stacked[e] = 1;
          if (f.v == root) root_children++;
          disc[w] = low[w] = time++;
          stack.push_back({w, e, 0});
        } else if (disc[w] < disc[f.v]) {     // an edge back to an ancestor, stacked once
          if (!edge_stacked[e]) edge_stack.push_back(e), edge_stacked[e] = 1;
          low[f.v] = std::min(low[f.v], disc[w]);
        }
        continue;
      }
      // f.v is finished: fold it into its parent.
      const Frame done = f;
      stack.pop_back();
      if (stack.empty()) break;
      const std::size_t v = stack.back().v, w = done.v;
      low[v] = std::min(low[v], low[w]);
      if (low[w] >= disc[v]) {
        if (v != root) out.cut[v] = 1;
        std::vector<std::size_t> block;
        while (!edge_stack.empty()) {
          const std::size_t e = edge_stack.back();
          edge_stack.pop_back();
          block.push_back(e);
          out.block_of[e] = out.blocks.size();
          if (e == done.parent_edge) break;
        }
        std::sort(block.begin(), block.end());
        out.blocks.push_back(std::move(block));
      }
      if (low[w] > disc[v]) out.bridge[done.parent_edge] = 1;
    }
    if (root_children >= 2) out.cut[root] = 1;
  }
  return out;
}

}  // namespace morphe::combinatorics
