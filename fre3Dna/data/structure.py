#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

import numpy as np

from .base import Base
from .strand import Strand

""" OxDNA structure class. A structure consists of strands and bases.
    It not only contains information from topology and configuration file,
    but also additional information: intended basepairs, nick positions, etc
"""


@dataclass
class Structure(object):
    bases: Dict[int, Base] = field(default_factory=dict)
    strands: Dict[int, Strand] = field(default_factory=dict)

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)

    def generate_from_oxDNA(self, top: Path, conf: Path):
        if len(self.bases) + len(self.strands):
            self.logger.error(
                "mistakes have been made. struc already generated")
            sys.exit(1)

        with top.open() as t:
            topology = t.readlines()
            n_bases, n_strands = (int(n) for n in topology.pop(0).split())
        with conf.open() as c:
            configuration = c.readlines()
            configuration = configuration[3:]  # discard temp, box, energy

        strand_ids = set()
        for i, line_top in enumerate(topology):
            line_conf = configuration[i]

            base_id = i
            strand_id_str, seq, p3_str, p5_str = line_top.split()
            strand_id = int(strand_id_str)
            strand_ids.add(strand_id)
            p3 = int(p3_str) if p3_str != "-1" else None
            p5 = int(p5_str) if p5_str != "-1" else None
            position = np.array(line_conf.split()[0:3])

            self.bases[base_id] = Base(
                id=base_id,
                seq=seq,
                strand_id=strand_id,
                struct=self,
                position=position,
                _p5_id=p5,
                _p3_id=p3,
            )

        for base in self.bases.values():
            base.p3 = self.bases[base._p3_id] if base._p3_id is not None else None
            base.p5 = self.bases[base._p5_id] if base._p5_id is not None else None

        for id in strand_ids:
            tour = list()
            tour_unsorted = [
                base for base in self.bases.values() if base.strand_id == id]
            start_base = next((b for b in tour_unsorted if b.p5 is None), None)
            if start_base is None:
                start_base = tour_unsorted[0]
                self.logger.warn("Circular strand encountered")
            tour.append(start_base)
            next_base = start_base.p3
            while next_base is not None:
                tour.append(next_base)
                next_base = next_base.p3
                if next_base is start_base:
                    break  # circular strand condition

            self.strands[id] = Strand(
                id=id,
                struct=self,
                tour=tour,
            )

    def assign_basepairs(self, basepairs):
        # assign base.across according to some external input file
        return
