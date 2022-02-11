#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import itertools

import numpy as np

# TODO: get correct value
DISTANCE_BB = 0.64


def get_scaffold_id(top):
    # TODO: properly handle multiscaffold or single oligo designs
    try:
        return int(top["strand"].mode())
    except TypeError:
        print("NOTE: design has no single longest strand, assuming scaffold_id=1")
        return 1


def get_nicks(top, conf, scaffold_id):
    """select only staple ends that form an actual nick"""
    staple_ends = top[((top["5p"] == -1) | (top["3p"] == -1)) & (top["strand"] != scaffold_id)]
    nicks = list()
    for pair in itertools.combinations(staple_ends.iterrows(), r=2):
        coords = [conf.iloc[p[0]] for p in pair]

        # NOTE: checking only for pairdistance works solely for crystal
        # TODO: alternative: check if bound to same scaffold position
        if distance(coords) < DISTANCE_BB:
            nicks.append([p[0] for p in pair])

    str_nicks = ""
    for nick in nicks:
        x, y = nick
        str_nicks += f"{x}-{y},"
    return str_nicks


def distance(coords):
    pos = [np.array([coor["x"], coor["y"], coor["z"]]) for coor in coords]
    return np.linalg.norm(pos[1] - pos[0])
