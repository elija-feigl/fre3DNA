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

BB_DIST = 0.63


@dataclass
class Structure(object):
    bases: Dict[int, Base] = field(default_factory=dict)
    strands: Dict[int, Strand] = field(default_factory=dict)
    basepairs: Dict[int, int] = field(default_factory=dict)
    nicks: Dict[int, int] = field(default_factory=dict)

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
            position = np.array([float(x) for x in line_conf.split()[0:3]])

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
            is_scaffold = False
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
                    is_scaffold = True
                    break  # circular strand condition

            self.strands[id] = Strand(
                id=id,
                struct=self,
                tour=tour,
                is_scaffold=is_scaffold,
            )

        if (len(self.bases) != n_bases) or (len(self.strands) != n_strands):
            self.logger.error(
                "mistakes have been made. number of strands and bases inconsistent")
            sys.exit(1)

    def assign_basepairs(self, forces: Path):
        with forces.open() as f:
            lines = f.readlines()

        is_pair = False
        for line in lines:
            if line.startswith("type"):
                if "mutual_trap" in line:
                    is_pair = True
                else:
                    is_pair = False

            if is_pair and line.startswith("particle"):
                particle = int(line.split()[-1])
            if is_pair and line.startswith("ref_particle"):
                ref_particle = int(line.split()[-1])
                self.basepairs[particle] = ref_particle
                is_pair = False

        for i, j in self.basepairs.items():
            self.bases[i].across = self.bases[j]

        self.logger.debug(f"Assigned {len(self.basepairs)/2} basepairs.")
        try:
            for i, j in self.basepairs.items():
                i_ = self.basepairs[j]
            if i_ != i:
                raise KeyError
        except KeyError:
            self.logger.error(
                "mistakes have been made. basepairs inconsistent")
            sys.exit(1)

    def categorise_structure(self):
        p5s = (s.tour[0] for s in self.strands.values() if not s.is_scaffold)
        p3s = [s.tour[-1] for s in self.strands.values() if not s.is_scaffold]

        for p5 in p5s:
            p3 = p5.across.p3.across
            if p3 in p3s:
                self.nicks[p5.id] = p3.id
                self.nicks[p3.id] = p5.id

    def generate_connectivity_graph(self, cutoff=2.5) -> list:
        weighted_edges = list()

        p5s = [s.tour[0] for s in self.strands.values() if not s.is_scaffold]
        p3s = [s.tour[-1] for s in self.strands.values() if not s.is_scaffold]

        for i, j in self.nicks.items():
            base5 = self.bases[i]
            if base5 in p5s:
                base3 = self.bases[j]
                edge = tuple([base5.strand_id, base3.strand_id, 0.])
                weighted_edges.append(edge)

        for p5 in p5s:
            for p3 in p3s:
                distance = np.linalg.norm(p5.position-p3.position)
                is_close = (distance <= cutoff * BB_DIST)
                not_self = (p5.strand_id != p3.strand_id)
                not_nick = (self.bases[self.nicks[p5.id]
                                       ].strand_id != p3.strand_id)

                if is_close and not_self and not_nick:
                    edge = tuple([p5.strand_id, p3.strand_id, distance])
                    weighted_edges.append(edge)

        return weighted_edges
