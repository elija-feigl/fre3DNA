#! /usr/bin/env python
# -*- coding: utf-8 -*-


import logging
from dataclasses import dataclass
from typing import List

import networkx as nx
import matplotlib.pyplot as plt

from .structure import Structure
from ..core.util import BB_DIST

""" fre3DNA connectivity graph class.
"""


@dataclass
class Graph(object):
    edges: List[tuple]
    struct: Structure

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)
        self.G = nx.MultiDiGraph()
        self.G.add_nodes_from(self.struct.scaffold_routing)
        self.G.add_weighted_edges_from(self.edges)

    def draw_graph(self, arrows=False):
        plt.figure(figsize=(8.0, 10.0))
        plt.axis('off')

        pos = nx.shell_layout(self.G, rotate=0)
        scaffold = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if d["weight"] == 0.]
        direct = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 0.5 < d["weight"]/BB_DIST < 1.5]
        spaceT = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 1.5 < d["weight"]/BB_DIST < 2.5]
        nx.draw_networkx_edges(self.G, pos, edgelist=scaffold,
                               edge_color="b", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=direct,
                               edge_color="g", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=spaceT,
                               edge_color="r", arrows=arrows, connectionstyle="arc3,rad=0.2")

        nx.draw_networkx_nodes(self.G, pos, node_size=10, node_color="k")

        self.logger.info("Plotting graph. Close to continue")
        plt.title("network_graph (b = sc, g < 1.5, r > 1.5)")
        plt.show()

    def network_stats(self):
        """ 5'3' stats
        """
        return False

    def draw_network_stats(self):
        stats = self.network_stats()
        plt.show()

    def reduce_graph(self):
        """ reduce graph for staple routing
            conditions:
                * every node can only have one in- one out-edge (except scaffold)
                * nodes must form triplets if possible
        """
        return
