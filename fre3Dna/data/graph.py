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
        self._expand_graph_data()

    def _expand_graph_data(self):

        for edge in self.G.edges(data=True):
            w = edge[2]["weight"]
            if w is None:
                edge[2]["is_scaffold"] = True
                edge[2]["is_nick"] = False
                edge[2]["is_53"] = False
                edge[2]["distance"] = 0.
                edge[2]["bb_multi"] = 0.
            else:
                edge[2]["is_scaffold"] = False
                edge[2]["is_nick"] = True if w < 0. else False
                edge[2]["is_53"] = True if w > 0. else False
                edge[2]["distance"] = abs(w)
                edge[2]["bb_multi"] = abs(w) / BB_DIST

    def draw_graph(self, arrows=True):
        plt.figure(figsize=(8.0, 10.0))
        plt.axis('off')

        pos = nx.shell_layout(self.G, rotate=0)
        scaffold = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if d["is_scaffold"]]
        nick = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if d["is_nick"]]
        short = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 0.01 < d["bb_multi"] < 0.8 and d["is_53"]]
        direct = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 0.8 < d["bb_multi"] < 1.2 and d["is_53"]]
        middle = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 1.2 < d["bb_multi"] < 1.7 and d["is_53"]]
        spaceT = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 1.7 < d["bb_multi"] < 2.3 and d["is_53"]]
        long = [(u, v) for (u, v, d) in self.G.edges(
            data=True) if 2.3 < d["bb_multi"] and d["is_53"]]
        nx.draw_networkx_edges(self.G, pos, edgelist=scaffold,
                               edge_color="b", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=nick,
                               edge_color="r", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=short,
                               edge_color="k", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=direct,
                               edge_color="g", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=middle,
                               edge_color="g", alpha=0.8, arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=spaceT,
                               edge_color="y", arrows=arrows, connectionstyle="arc3,rad=0.2")
        nx.draw_networkx_edges(self.G, pos, edgelist=long,
                               edge_color="k", alpha=0.3, arrows=arrows, connectionstyle="arc3,rad=0.2")

        nx.draw_networkx_nodes(self.G, pos, node_size=10, node_color="k")

        self.logger.info("Plotting graph. Close to continue")
        plt.title("network_graph (b = sc, g < 1.5, r > 1.5)")
        plt.show()

    def get_edges(self, typ=""):
        if typ == "53":
            return [(u, v, d) for (u, v, d) in self.G.edges(data=True) if d["is_53"]]
        elif typ == "nick":
            return [(u, v, d) for (u, v, d) in self.G.edges(data=True) if d["is_nick"]]
        elif typ == "scaffold":
            return [(u, v, d) for (u, v, d) in self.G.edges(data=True) if d["is_scaffold"]]
        elif typ == "":
            return [(u, v, d) for (u, v, d) in self.G.edges(data=True)]
        else:
            raise KeyError

    def network_stats(self):
        """ 5'3' stats
        """
        n_staple_edges = len(self.get_edges("53"))
        self.logger.info(f"staple_edges: {n_staple_edges}")
        n_scaffold_edges = len(self.get_edges("scaffold"))
        self.logger.info(f"scaffold_edges: {n_scaffold_edges}")
        n_nick_edges = len(self.get_edges("nick"))
        self.logger.info(f"nick_edges: {n_nick_edges}")

        bb_multiples = [d["bb_multi"] for (_, _, d) in self.get_edges("53")]
        return bb_multiples

    def draw_network_stats(self):
        bb_multiples = self.network_stats()
        f, axs = plt.subplots(1, 2, figsize=(
            8, 10), sharey=False, tight_layout=True)
        hist, bins = np.histogram(bb_multiples, range=(
            0, 3.5), bins=35, density=False)
        cum = np.cumsum(hist)
        bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]

        axs[0].plot(bincentres, hist)
        axs[1].plot(bincentres, cum)

        f.suptitle("53_pairs_distr")
        plt.show()

    def reduce_graph_simple(self):
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
        for node_ids in self.G.nodes:
            all_edges = [(node_ids, v, d["weight"])
                         for v, d in self.G[node_ids].items() if d["is_53"]]
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

        n_53_edges = len(self.get_edges("53"))
        self.logger.info(f"all: {n_53_edges}")
        self.logger.info(f"simple: {len(edges)}")

        all_nicks = [(u, v, d["weight"])
                     for (u, v, d) in self.get_edges("nick")]
        nicks = list()
        for nick in all_nicks:
            is_53start = nick[0] in [u for u, _, _ in edges]
            is_53end = nick[1] in [v for _, v, _ in edges]
            if not is_53start and not is_53end:
                nicks.append(nick)
        edges += nicks

        edges += [(u, v, d["weight"])
                  for (u, v, d) in self.get_edges("scaffold")]

        self.G.clear_edges()
        self.G.add_weighted_edges_from(edges)
        self._expand_graph_data()

    def get_routing(self, max_bb_multi=2.3):
        pairs = [(u, v, d["distance"]) for (u, v, d) in self.get_edges(
            "53") if d["bb_multi"] < max_bb_multi]
        # NOTE: add penalty to force filler base
        nicks = [(u, v, (d["distance"] + BB_DIST)) for (u, v, d) in self.get_edges(
            "nick") if d["bb_multi"] < max_bb_multi]
        return nicks + pairs
