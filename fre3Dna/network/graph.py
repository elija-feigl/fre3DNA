#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html

import logging
from dataclasses import dataclass
from typing import List
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from ..core.util import BB_DIST
from ..data.structure import Structure

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

    def _expand_graph_data(self, _g=None):
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

        if _g is None:
            _g = self.G

        for edge in _g.edges(data=True):
            w = edge[2]["weight"]
            edge[2]["is_nick"] = True if w < 0. else False
            edge[2]["is_53"] = True if w > 0. else False
            edge[2]["distance"] = abs(w)
            edge[2]["bb_multi"] = abs(w) / BB_DIST

            edge[2]["sc_distance"] = scaffold_distance(edge[0], edge[1])
            edge[2]["penalty"] = (abs(w) / BB_DIST + BB_DIST * (w < 0.) +
                                  scaffold_distance(edge[0], edge[1], True) * BB_DIST)

    def draw_graph(self, arrows=True):
        plt.figure(figsize=(8.0, 10.0))
        plt.axis('off')

        pos = nx.shell_layout(self.G, rotate=0)
        _ = nx.draw_networkx_labels(self.G, pos=pos)
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
        if typ in ["53", "nick"]:
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
        n_nick_edges = len(self.get_edges("nick"))
        self.logger.info(f"nick_edges: {n_nick_edges}")
        self.logger.info(f"isolates: {len(list(nx.isolates(self.G)))}")

        edges = self.get_edges(typ="")
        avg_penalty = sum(d["penalty"] for _, _, d in edges) / len(edges)
        self.logger.info(f"avg_penalty: {avg_penalty}")
        avg_weight = sum(d["weight"] for _, _, d in edges) / len(edges)
        self.logger.info(f"avg_weight: {avg_weight}")

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

    def split_long_paths(self):
        # TODO: weak_components might be faster
        longest_path = nx.dag_longest_path(self.G, weight="n")
        while len(longest_path) > 11:  # TODO: replace with cumulative n < 80
            edges = list()
            for i, u in enumerate(longest_path):
                if (1 < i < len(longest_path)-3):  # NOTE: avoid short strands
                    v = longest_path[(i+1) % len(longest_path)]
                    penalty = self.G[u][v]["penalty"]
                    edges.append((u, v, penalty))
            max_edge = max(edges, key=lambda e: e[2])
            self.G.remove_edge(max_edge[0], max_edge[1])
            longest_path = nx.dag_longest_path(self.G, weight="n")

    def split_cylces(self):
        for cycle in nx.simple_cycles(self.G):
            edges = list()
            for i, u in enumerate(cycle):
                v = cycle[(i+1) % len(cycle)]
                penalty = self.G[u][v]["penalty"]
                edges.append((u, v, penalty))
            max_edge = max(edges, key=lambda e: e[2])
            self.G.remove_edge(max_edge[0], max_edge[1])

    def reduce_isolates(self, possible_edges):
        """ reduce number of isolates, improved version
                * iterate over isolated nodes and try to add node
                * iterate over double nodes and try to add node
                * cut long staples
        """
        def try_replace_out(_g: nx.DiGraph, edge) -> bool:
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

        def try_replace_in(_g: nx.DiGraph, edge) -> bool:
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

        _g = self.G.copy()

        # single element
        # outgoing
        for isolate in list(nx.isolates(_g)):
            candidates_out = sorted([(u, v, d) for (
                u, v, d) in possible_edges if isolate == u and abs(u-v) > 1], key=lambda e: e[2]["penalty"])
            while candidates_out:
                candidate = candidates_out.pop()
                if try_replace_out(_g, candidate):
                    break

        # incomming
        for isolate in list(nx.isolates(_g)):
            candidates_in = sorted([(u, v, d) for (
                u, v, d) in possible_edges if isolate == v and abs(u-v) > 1], key=lambda e: e[2]["penalty"])
            while candidates_in:
                candidate = candidates_in.pop()
                if try_replace_in(_g, candidate):
                    break

        self.G = _g
        self._expand_graph_data()

    def optimize_routing_metropolis(self, possible_edges, steps=1e3):
        """ metropolis like algorithm that reduces the average edge penalty.
            possible steps:
            * remove an edge
            * add an edge
            * replace an edge
            rules:
            * must not close a ring
        """

        def get_score(_g: nx.DiGraph, measure="penalty") -> float:
            edges = [(u, v, d) for (u, v, d) in _g.edges(data=True)]
            avg_penalty = sum(d[measure]
                              for _, _, d in edges) / len(edges) - 1.

            isolate_penalty = sum(
                3. for c in nx.weakly_connected_components(_g) if len(c) == 1)
            disolate_penalty = sum(
                1. for c in nx.weakly_connected_components(_g) if len(c) == 2)
            try:
                next(nx.simple_cycles(_g))
                cycle_penalty = 100.
            except StopIteration:
                cycle_penalty = 0.
            # length penalty
            score = avg_penalty + isolate_penalty + disolate_penalty + cycle_penalty
            return score

        def get_node_weights(_g: nx.DiGraph):
            node_weights = dict()
            for component in nx.weakly_connected_components(_g):
                weight = 1./len(component)
                for node in component:
                    node_weights[node] = weight
            return node_weights

        # TODO: don't work with full tree to save time
        _g: nx.DiGraph = self.G.copy()
        score = get_score(_g)
        node_weights = get_node_weights(_g)

        for _ in range(steps):
            # select node using weighted random choice -> get_node_weights()
            node = random.choices(list(node_weights.keys()),
                                  weights=list(node_weights.values())).pop()
            # select move type:
            move_type = random.choice(["in", "out"])
            # print(list(nx.isolates(_g)))

            # prepare move
            _h: nx.DiGraph = _g.copy()
            if move_type == "in":
                old_edge = next(((u, v)
                                for (u, v) in _h.edges() if v == node), None)
                candidates = [(u, v, d) for (u, v, d)
                              in possible_edges if ((v == node) and (abs(u-v) > 1) and ((u, v) != old_edge))]
                if not candidates:
                    # self.logger.info(
                    #    f"Rejected move {move_type} for node {node}. impossible")
                    continue
                # TODO: maybe add penalty as a weight for the choice?
                new_edge = random.choice(candidates)
                old_edge_comp = next(((u, v)
                                      for (u, v) in _h.edges() if u == new_edge[0]), None)

            elif move_type == "out":
                old_edge = next(((u, v)
                                for (u, v) in _h.edges() if u == node), None)
                candidates = [(u, v, d) for (u, v, d)
                              in possible_edges if ((u == node) and (abs(u-v) > 1) and ((u, v) != old_edge))]
                if not candidates:
                    # self.logger.info(
                    #    f"Rejected move {move_type} for node {node}. impossible")
                    continue
                # TODO: maybe add penalty as a weight for the choice?
                new_edge = random.choice(candidates)
                old_edge_comp = next(((u, v)
                                      for (u, v) in _h.edges() if v == new_edge[1]), None)

            # update move
            if old_edge is not None:
                _h.remove_edge(old_edge[0], old_edge[1])
            if old_edge_comp is not None:
                _h.remove_edge(old_edge_comp[0], old_edge_comp[1])
            _h.add_edge(new_edge[0], new_edge[1], weight=new_edge[2]["weight"])
            self._expand_graph_data(_h)

            # metropolis decision
            new_score = get_score(_h)
            alpha = min(1., np.exp(-new_score)/np.exp(-score))

            if alpha >= random.random():
                _g = _h
                node_weights = get_node_weights(_g)
                score = new_score
                # self.logger.info(
                #    f"Accepted move: {old_edge}->{new_edge[:2]}; {move_type} with probability {alpha}")
            else:
                _h = _g
                # self.logger.info(
                #    f"Rejected move: {old_edge}->{new_edge[:2]}; {move_type} with probability {alpha}")

        self.G = _g

    def reduce_graph_direct(self, reduce_iso=True):
        """ reduce graph for staple routing
            conditions:
                * every node can only have one in- one out-edge (except scaffold)
                # * nodes must form triplets if possible

            simple:
                * pick the best outgoing per node
                * pick best incomming per node
        """
        possible_edges = self.get_edges(typ="53") + self.get_edges(typ="nick")

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

        # nicks
        all_nicks = [(u, v, d["weight"])
                     for (u, v, d) in self.get_edges("nick")]
        nicks = list()
        for nick in all_nicks:
            is_53start = nick[0] in [u for u, _, _ in edges]
            is_53end = nick[1] in [v for _, v, _ in edges]
            if not is_53start and not is_53end:
                nicks.append(nick)
        edges += nicks

        self.G.clear_edges()
        self.G.add_weighted_edges_from(edges)
        self._expand_graph_data()

        self.split_cylces()
        self.split_long_paths()

        if reduce_iso:
            self.reduce_isolates(possible_edges=possible_edges)

    def reduce_graph_reverse(self, is_simple=True, reduce_iso=True):
        """ reduce graph for staple routing, improved version
            conditions:
                * every node can only have one in- one out-edge
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
                _g[u][v]["n"] = 7  # TODO: node size?

        def edge_step_recursive(_g: nx.DiGraph, edge):
            raise NotImplementedError

        _53_edges = [(u, v, d["weight"], d["penalty"])
                     for (u, v, d) in self.get_edges(typ="53")]
        nick_edges = [(u, v, d["weight"], d["penalty"])
                      for (u, v, d) in self.get_edges(typ="nick")]
        possible_edges = self.get_edges(typ="53") + self.get_edges(typ="nick")

        _g = nx.DiGraph()
        _g.add_nodes_from(self.struct.scaffold_routing)

        for edge in sorted(_53_edges, key=lambda e: e[3]) + sorted(nick_edges, key=lambda e: e[3]):
            if is_simple:
                edge_step(_g, edge)
            else:
                edge_step_recursive(_g, edge)

        self.G = _g
        self._expand_graph_data()

        # if reduce_iso:
        self.reduce_isolates(possible_edges=possible_edges)
        self.split_cylces()
        self.split_long_paths()

        # self.optimize_routing_metropolis(
        #    possible_edges=possible_edges, steps=int(1e6))

    def get_routing(self, max_bb_multi=2.3):
        pairs = [(u, v, d["distance"]) for (u, v, d) in self.get_edges(
            "53") if d["bb_multi"] < max_bb_multi]
        # NOTE: add penalty to force filler base
        nicks = [(u, v, (d["distance"] + BB_DIST)) for (u, v, d) in self.get_edges(
            "nick") if d["bb_multi"] < max_bb_multi]
        return nicks + pairs
