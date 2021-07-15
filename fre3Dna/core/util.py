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

import numpy as np

BB_DIST = 0.7564  # oxdna1 0.7525
OXDNA_L = 0.8518  # nm


def basepair_trap(particle: int, ref_particle: int) -> str:
    trap = ["{",
            "type = mutual_trap",
            f"particle = {particle}",
            "stiff = 0.9",
            "r0 = 1.2",
            f"ref_particle = {ref_particle}",
            "PBC=1",
            "}"]
    return "\n".join(trap)


def bb_position_2_base_position(bb_position, normversor, bb2base_versor) -> np.ndarray:
    """ base position calculated from backbone position for oxDNA2 (w. grooving)
    """
    thrid_versor = np.cross(normversor, bb2base_versor)
    return bb_position + 0.34 * bb2base_versor - 0.3408 * thrid_versor
