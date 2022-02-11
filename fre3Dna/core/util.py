#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import numpy as np

BB_DIST = 0.7564  # oxdna1 0.7525
OXDNA_L = 0.8518  # nm


def basepair_trap(particle: int, ref_particle: int) -> str:
    trap = [
        "{",
        "type = mutual_trap",
        f"particle = {particle}",
        "stiff = 0.9",
        "r0 = 1.2",
        f"ref_particle = {ref_particle}",
        "PBC=1",
        "}",
    ]
    return "\n".join(trap)


def bb_position_2_base_position(bb_position, normversor, bb2base_versor):
    """base position calculated from backbone position for oxDNA2 (w. grooving)"""
    thrid_versor = np.cross(normversor, bb2base_versor)
    return bb_position + 0.34 * bb2base_versor - 0.3408 * thrid_versor
