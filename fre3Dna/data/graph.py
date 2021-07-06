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
        def scaffold_distance(node1, node2, relative=False):
            routing = self.struct.scaffold_routing
            idx1 = routing.index(node1)
            idx2 = routing.index(node2)
            dist = min((idx1 - idx2) % len(routing),
                       (idx2 - idx1) % len(routing)) - 1
            if relative:
                return dist / len(routing)
            else:
                return dist

        for edge in self.G.edges(data=True):
            w = edge[2]["weight"]
            if w is None:
                edge[2]["is_scaffold"] = True
                edge[2]["is_nick"] = False
                edge[2]["is_53"] = False
                edge[2]["distance"] = 0.
                edge[2]["bb_multi"] = 0.
                edge[2]["penalty"] = 0.
            else:
                edge[2]["is_scaffold"] = False
                edge[2]["is_nick"] = True if w < 0. else False
                edge[2]["is_53"] = True if w > 0. else False
                edge[2]["distance"] = abs(w)
                edge[2]["bb_multi"] = abs(w) / BB_DIST

                edge[2]["sc_distance"] = scaffold_distance(edge[0], edge[1])
                edge[2]["penalty"] = (
                    abs(w) + BB_DIST * (w < 0.) + scaffold_distance(edge[0], edge[1], True) * BB_DIST)

    def draw_graph(self, arrows=True):
        plt.figure(figsize=(8.0, 10.0))
        plt.axis('off')

        pos = nx.shell_layout(self.G, rotate=0)
        _ = nx.draw_networkx_labels(self.G, pos=pos)
        # scaffold = [(u, v) for (u, v, d) in self.G.edges(
        #    data=True) if d["is_scaffold"]]
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
        # nx.draw_networkx_edges(self.G, pos, edgelist=scaffold,
        #                      edge_color="b", arrows=arrows, connectionstyle="arc3,rad=0.2")
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
        if typ in ["53", "nick", "scaffold"]:
            return [(u, v, d) for (u, v, d) in self.G.edges(data=True) if d[f"is_{typ}"]]
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
        out_edges = list()
        # single exit
        for node_ids in self.G.nodes:
            all_edges = [(node_ids, v, d["weight"], d["penalty"])
                         for v, d in self.G[node_ids].items() if d["is_53"]]
            if all_edges:
                best_edge = min(all_edges, key=lambda e: e[3])
                out_edges.append(best_edge[:-1])

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

    def reduce_graph_isofree(self):
        """ reduce graph for staple routing, improved version
            conditions:
                * every node can only have one in- one out-edge (except scaffold)
                * all staples between 21 and 84 bases
            approach:
                * perform simple routing
                * iterate over isolated nodes and try to add node
                * iterate over double nodes and try to add node
                * cut long staples
        """
        def try_replace_out(_g, edge) -> bool:
            target = edge[1]
            sources_old = [u for (u, v) in _g.edges() if v == target]
            if sources_old:
                source_old = sources_old[0]
            else:  # NOTE: if old source has been remodef by isolator fix
                _g.add_edge(edge[0], edge[1], weight=edge[2]["weight"])
                return True

            if [u for (u, v) in _g.edges() if v == source_old]:
                _g.remove_edge(source_old, target)
                _g.add_edge(edge[0], edge[1], weight=edge[2]["weight"])
                return True
            return False

        def try_replace_in(_g, edge) -> bool:
            source = edge[0]
            targets_old = [v for (u, v) in _g.edges() if u == source]
            if targets_old:
                target_old = targets_old[0]
            else:  # NOTE: if old source has been remodef by isolator fix
                _g.add_edge(edge[0], edge[1], weight=edge[2]["weight"])
                return True

            if [v for (u, v) in _g.edges() if u == target_old]:
                _g.remove_edge(source, target_old)
                _g.add_edge(edge[0], edge[1], weight=edge[2]["weight"])
                return True
            return False

        old_edges = self.get_edges(typ="53") + self.get_edges(typ="nick")
        self.reduce_graph_simple()

        _g = self.G.copy()
        scaffold_edges = [(u, v, d["weight"])
                          for (u, v, d) in self.get_edges("scaffold")]
        for u, v, _ in scaffold_edges:
            _g.remove_edge(u, v)

        # single element
        # outgoing
        for isolate in list(nx.isolates(_g)):
            candidates_out = sorted([(u, v, d) for (
                u, v, d) in old_edges if isolate == u and abs(u-v) > 1], key=lambda e: e[2]["penalty"])
            while candidates_out:
                candidate = candidates_out.pop()
                if try_replace_out(_g, candidate):
                    break

        # incomming
        for isolate in list(nx.isolates(_g)):
            candidates_in = sorted([(u, v, d) for (
                u, v, d) in old_edges if isolate == v and abs(u-v) > 1], key=lambda e: e[2]["penalty"])
            while candidates_in:
                candidate = candidates_in.pop()
                if try_replace_in(_g, candidate):
                    break

        # cycles
        for cycle in nx.simple_cycles(_g):
            edges = list()
            for i, u in enumerate(cycle):
                v = cycle[(i+1) % len(cycle)]
                penalty = _g[u][v]["penalty"]
                edges.append((u, v, penalty))
            max_edge = max(edges, key=lambda e: e[2])
            _g.remove_edge(max_edge[0], max_edge[1])

        # long
        # TODO

        self.G = _g
        self.G.add_weighted_edges_from(scaffold_edges)

        self._expand_graph_data()

    def reduce_graph_reverse_simple(self):
        """ reduce graph for staple routing, improved version
            conditions:
                * every node can only have one in- one out-edge (except scaffold)
                * all staples between 21 and 84 bases
            approach:
                iteratively add edges to empty graph
        """
        def edge_step(_g: nx.DiGraph, edge):
            u, v, weight, penalty = edge

            u_is_free = not bool(_g[u])
            v_is_free = not bool([w for _, w in _g.edges if w == v])

            if u_is_free and v_is_free:
                _g.add_edge(u, v, weight=weight)
                _g[u][v]["penalty"] = penalty

        _53_edges = [(u, v, d["weight"], d["penalty"])
                     for (u, v, d) in self.get_edges(typ="53")]
        nick_edges = [(u, v, d["weight"], d["penalty"])
                      for (u, v, d) in self.get_edges(typ="nick")]
        scaffold_edges = [(u, v, d["weight"])
                          for (u, v, d) in self.get_edges("scaffold")]

        _g = nx.DiGraph()
        _g.add_nodes_from(self.struct.scaffold_routing)

        for edge in sorted(_53_edges, key=lambda e: e[3]) + sorted(nick_edges, key=lambda e: e[3]):
            edge_step(_g, edge)

        for cycle in nx.simple_cycles(_g):
            edges = list()
            for i, u in enumerate(cycle):
                v = cycle[(i+1) % len(cycle)]
                penalty = _g[u][v]["penalty"]
                edges.append((u, v, penalty))
            max_edge = max(edges, key=lambda e: e[2])
            _g.remove_edge(max_edge[0], max_edge[1])

        self.G = _g
        self.G.add_weighted_edges_from(scaffold_edges)
        self._expand_graph_data()

    def get_routing(self, max_bb_multi=2.3):
        pairs = [(u, v, d["distance"]) for (u, v, d) in self.get_edges(
            "53") if d["bb_multi"] < max_bb_multi]
        # NOTE: add penalty to force filler base
        nicks = [(u, v, (d["distance"] + BB_DIST)) for (u, v, d) in self.get_edges(
            "nick") if d["bb_multi"] < max_bb_multi]
        return nicks + pairs
