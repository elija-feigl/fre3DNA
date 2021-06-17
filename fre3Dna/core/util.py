#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
