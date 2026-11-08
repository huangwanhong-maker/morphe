"""Reference centralities for check_centrality.cpp, computed by networkx.

Writes centrality.txt: for each graph, a line "graph <name> <n> <m>", m lines "u v w" (undirected edges and their weights), then
n lines "v betweenness total_cost harmonic reached" with networkx's values:
- betweenness: betweenness_centrality(G, normalized=False, weight="weight"), each unordered pair once;
- total_cost, harmonic, reached: from single_source_dijkstra_path_length, the source excluded.

Run with networkx 3.4.2 and numpy (seeded): python make_centrality_fixtures.py
"""

import os

import networkx as nx
import numpy as np


def graphs():
    rng = np.random.default_rng(20261108)
    # A random graph with real weights: least paths are unique, almost surely.
    g = nx.gnp_random_graph(40, 0.12, seed=7)
    for u, v in g.edges:
        g[u][v]["weight"] = float(np.round(rng.uniform(1, 10), 6))
    g = g.subgraph(max(nx.connected_components(g), key=len)).copy()
    g = nx.convert_node_labels_to_integers(g)
    yield "random-real", g
    # A grid with unit weights: many least paths of equal length.
    h = nx.convert_node_labels_to_integers(nx.grid_2d_graph(6, 5))
    for u, v in h.edges:
        h[u][v]["weight"] = 1.0
    yield "grid-6x5", h
    # A random geometric graph with small integer weights: ties among paths of different shapes.
    k = nx.random_geometric_graph(30, 0.32, seed=11)
    for u, v in k.edges:
        k[u][v]["weight"] = float(rng.integers(1, 4))
    k = k.subgraph(max(nx.connected_components(k), key=len)).copy()
    k = nx.convert_node_labels_to_integers(k)
    yield "geometric-integer", k


def main():
    out = []
    for name, g in graphs():
        n, edges = g.number_of_nodes(), list(g.edges(data="weight"))
        out.append(f"graph {name} {n} {len(edges)}")
        out += [f"{u} {v} {w!r}" for u, v, w in edges]
        btw = nx.betweenness_centrality(g, normalized=False, weight="weight")
        for v in range(n):
            d = nx.single_source_dijkstra_path_length(g, v, weight="weight")
            others = [x for t, x in d.items() if t != v]
            out.append(f"{v} {btw[v]!r} {sum(others)!r} {sum(1 / x for x in others if x > 0)!r} {len(others)}")
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "centrality.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# networkx {nx.__version__}, numpy {np.__version__}; written by make_centrality_fixtures.py\n")
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
