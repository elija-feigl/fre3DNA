#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
