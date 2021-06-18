#! /usr/bin/env python
# -*- coding: utf-8 -*-


import logging
from dataclasses import dataclass
from typing import List

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from ..core.util import BB_DIST
from .structure import Structure

""" fre3DNA connectivity graph class.
"""


@dataclass
class Graph(object):
    edges: List[tuple]
    struct: Structure

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)

        self.G = nx.DiGraph()
        self.G.add_nodes_from(self.struct.scaffold_routing)
        self.G.add_weighted_edges_from(self.edges)

        self.G_scaffold = nx.DiGraph()
        self.G_scaffold.add_nodes_from(self.struct.scaffold_routing)
        scaffold_edges = [(u, v, w) for (u, v, w) in self.edges if w == 0.]
        self.G_scaffold.add_weighted_edges_from(scaffold_edges)

        self.G_staple = nx.DiGraph()
        self.G_staple.add_nodes_from(self.struct.scaffold_routing)
        staple_edges = [(u, v, w) for (u, v, w) in self.edges if w != 0.]
        self.G_staple.add_weighted_edges_from(staple_edges)

    def recompose_G(self):
        self.G = nx.compose(self.G_scaffold, self.G_staple)

    def draw_graph(self, arrows=True):
        self.recompose_G()
        plt.figure(figsize=(8.0, 10.0))
        plt.axis('off')

        pos = nx.shell_layout(self.G, rotate=0)
        scaffold = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if d["weight"] == 0.]
        short = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 0. < d["weight"]/BB_DIST < 0.8]
        direct = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 0.8 < d["weight"]/BB_DIST < 1.2]
        middle = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 1.2 < d["weight"]/BB_DIST < 1.7]
        spaceT = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 1.7 < d["weight"]/BB_DIST < 2.3]
        long = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 2.3 < d["weight"]/BB_DIST]
        nx.draw_networkx_edges(self.G, pos, edgelist=scaffold,
                               edge_color="b", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=short,
                               edge_color="k", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=direct,
                               edge_color="g", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=middle,
                               edge_color="y", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=spaceT,
                               edge_color="r", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=long,
                               edge_color="k", alpha=0.3, arrows=arrows, connectionstyle="arc3,rad=0.2")

        nx.draw_networkx_nodes(self.G, pos, node_size=10, node_color="k")

        self.logger.info("Plotting graph. Close to continue")
        plt.title("network_graph (b = sc, g < 1.5, r > 1.5)")
        plt.show()

    def network_stats(self):
        """ 5'3' stats
        """
        self.logger.info(f"staple_edges: {len(self.G_staple.edges)}")
        self.logger.info(f"scaffold_edges: {len(self.G_scaffold.edges)}")

        bb_multiples = [d["weight"]/BB_DIST
                        for (_, _, d) in self.G_staple.edges(data=True)]
        return bb_multiples

    def draw_network_stats(self):
        bb_multiples = self.network_stats()
        f, axs = plt.subplots(1, 2, figsize=(
            8, 10), sharey=False, tight_layout=True)
        hist, bins = np.histogram(bb_multiples, range=(
            0, 3.5), bins=60, density=False)
        cum = np.cumsum(hist)
        bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]

        axs[0].plot(bincentres, hist)
        axs[1].plot(bincentres, cum)

        f.suptitle("53_pairs_distr")
        plt.show()

    def reduce_graph(self):
        """ reduce graph for staple routing
            conditions:
                * every node can only have one in- one out-edge (except scaffold)
                # * nodes must form triplets if possible

            simple:
                * pick the best outgoing per node
                * pick best incomming per node
        """
        # simple reduction
        out_edges = list()
        # single exit
        for node_ids in self.G_staple.nodes:
            all_edges = [(node_ids, v, d["weight"])
                         for v, d in self.G_staple[node_ids].items()]
            if all_edges:
                best_edge = min(all_edges, key=lambda e: e[2])
                out_edges.append(best_edge)

        # single enter
        edges = list()
        clear_ids = set()
        for edge in sorted(out_edges, key=lambda e: (e[1], e[2])):
            enter_id = edge[1]
            if enter_id not in clear_ids:
                clear_ids.add(enter_id)
                edges.append(edge)

        self.logger.info(f"all: {len(self.G_staple.edges)}")
        self.logger.info(f"simple: {len(edges)}")

        self.G_staple.clear_edges()
        self.G_staple.add_weighted_edges_from(edges)
