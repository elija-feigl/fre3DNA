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
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

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
    """ select only staple ends that form an actual nick
    """
    staple_ends = top[
        ((top["5p"] == -1) | (top["3p"] == -1))
        & (top["strand"] != scaffold_id)
    ]
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
        str_nicks += "{}-{},".format(x, y)
    return str_nicks


def distance(coords):
    pos = [
        np.array([coor["x"], coor["y"], coor["z"]])
        for coor in coords
    ]
    return np.linalg.norm(pos[1] - pos[0])
