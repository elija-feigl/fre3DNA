#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys
from dataclasses import dataclass, field
from operator import attrgetter
from pathlib import Path
from typing import Dict, List

import numpy as np

from ..core.util import basepair_trap
from .base import Base
from .strand import Strand

""" OxDNA structure class. A structure consists of strands and bases.
    It not only contains information from topology and configuration file,
    but also additional information: intended basepairs, nick positions, etc
"""

BB_DIST = 0.7564  # oxdna1 0.7525
OXDNA_L = 0.8518  # nm


@dataclass
class Structure(object):
    bases: Dict[int, Base] = field(default_factory=dict)
    strands: Dict[int, Strand] = field(default_factory=dict)
    basepairs: Dict[int, int] = field(default_factory=dict)
    nicks: Dict[int, int] = field(default_factory=dict)
    scaffold_routing: List[int] = field(default_factory=list)

    header: List[str] = field(default_factory=list)

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
            lines = c.readlines()
            self.header = [line[:-1] for line in lines[:3]]
            configuration = lines[3:]

        strand_ids = set()
        for i, line_top in enumerate(topology):
            line_conf = configuration[i]

            base_id = i
            strand_id_str, seq, p3_str, p5_str = line_top.split()
            strand_id = int(strand_id_str)
            strand_ids.add(strand_id)
            p3 = int(p3_str) if p3_str != "-1" else None
            p5 = int(p5_str) if p5_str != "-1" else None

            base_data = [float(x) for x in line_conf.split()]
            position = np.array(base_data[0:3])
            bb2base_versor = np.array(base_data[3:6])
            normversor = np.array(base_data[6:9])
            velocity = np.array(base_data[9:12])
            angular_velocity = np.array(base_data[12:15])

            self.bases[base_id] = Base(
                id=base_id,
                seq=seq,
                strand_id=strand_id,
                struct=self,
                position=position,
                bb2base_versor=bb2base_versor,
                normversor=normversor,
                velocity=velocity,
                angular_velocity=angular_velocity,
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

        self.logger.debug(f"Assigned {len(self.basepairs)//2} basepairs.")
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
        nicks_staple = dict()

        for p5 in p5s:
            p3 = p5.across.p3.across
            if p3 in p3s:
                self.nicks[p5.id] = p3.id
                self.nicks[p3.id] = p5.id
                nicks_staple[p5.strand_id] = p3.strand_id

        first_staple_id = min(nicks_staple.keys())
        self.scaffold_routing.append(first_staple_id)
        next_staple_id = nicks_staple[first_staple_id]
        while next_staple_id != first_staple_id:
            self.scaffold_routing.append(next_staple_id)
            next_staple_id = nicks_staple[next_staple_id]

    def generate_connectivity_graph(self, cutoff=2.5) -> list:
        """ generate a list of weighted edges for graph construction.
            an edge runs from 5' to 3' at a potential junction (scaffold-CO or staple-staple)
            its weight is the distance between the end or 0 for scaffold-CO
        """
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
                distance = np.linalg.norm(p5.bb_position()-p3.bb_position())
                is_close = (distance <= cutoff * BB_DIST)
                not_self = (p5.strand_id != p3.strand_id)
                not_nick = (self.bases[self.nicks[p5.id]
                                       ].strand_id != p3.strand_id)

                if is_close and not_self and not_nick:
                    edge = tuple([p5.strand_id, p3.strand_id, distance])
                    weighted_edges.append(edge)

        return weighted_edges

    def write_structure(self, prefix: str):
        top_file = Path(f"{prefix}.top")
        conf_file = Path(f"{prefix}.oxdna")
        trap_file = Path(f"{prefix}_bp-trap.txt")
        n_bakk = 1
        while top_file.exists() or conf_file.exists():
            top_file = Path(f"{prefix}_{n_bakk}.top")
            conf_file = Path(f"{prefix}_{n_bakk}.oxdna")
            trap_file = Path(f"{prefix}_bp-trap_{n_bakk}.txt")
            n_bakk += 1
        self.logger.info(f"Writing files with prefix {prefix}_{n_bakk}")

        topology = [f"{len(self.bases)} {len(self.strands)}"]
        configuration = self.header.copy()
        for _, strand in sorted(self.strands.items()):
            # NOTE: for now assume that base.id is correct (renumberd on insert etc.)
            for base in sorted(strand.tour, key=attrgetter("id")):
                topology.append(base.top())
                configuration.append(base.conf())

        with top_file.open(mode="w") as top:
            top.write("\n".join(topology))

        with conf_file.open(mode="w") as conf:
            conf.write("\n".join(configuration))

        with trap_file.open(mode="w") as trap:
            forces = [basepair_trap(i, j) for i, j in self.basepairs.items()]
            trap.write("\n".join(forces))
