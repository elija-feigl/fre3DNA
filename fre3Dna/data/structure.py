#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys
from dataclasses import dataclass, field
from operator import attrgetter
from pathlib import Path
from typing import Dict, List

import numpy as np

from ..core.util import basepair_trap, bb_position_2_base_position, BB_DIST
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
                nicks_staple[p5.strand_id] = p3.strand_id

        first_staple_id = min(nicks_staple.keys())
        self.scaffold_routing.append(first_staple_id)
        next_staple_id = nicks_staple[first_staple_id]
        while next_staple_id != first_staple_id:
            self.scaffold_routing.append(next_staple_id)
            next_staple_id = nicks_staple[next_staple_id]

    def generate_connectivity(self, cutoff=2.5) -> list:
        """ generate a list of weighted edges for graph construction.
            an edge runs from 5' to 3' at a potential junction (scaffold-CO or staple-staple)
            its weight is the distance between the end or 0 for scaffold-CO
        """
        weighted_edges = list()
        p5_bases = [s.tour[0]
                    for s in self.strands.values() if not s.is_scaffold]
        p3_bases = [s.tour[-1]
                    for s in self.strands.values() if not s.is_scaffold]

        for p3, p5 in self.nicks.items():
            base5 = self.bases[p5]
            base3 = self.bases[p3]
            distance = np.linalg.norm(
                base5.bb_position()-base3.bb_position())
            edge_scaffold = tuple([base3.strand_id, base5.strand_id, 0.])
            edge_staple = tuple([base5.strand_id, base3.strand_id, -distance])

            weighted_edges.append(edge_scaffold)
            weighted_edges.append(edge_staple)

        for base5 in p5_bases:
            for base3 in p3_bases:
                distance = np.linalg.norm(
                    base5.bb_position()-base3.bb_position())
                is_close = (distance <= cutoff * BB_DIST)
                not_self = (base5.strand_id != base3.strand_id)
                not_nick = (self.bases[self.nicks[base5.id]
                                       ].strand_id != base3.strand_id)

                if is_close and not_self and not_nick:
                    edge = tuple([base5.strand_id, base3.strand_id, distance])
                    weighted_edges.append(edge)
        return weighted_edges

    def _reorder_structure(self):
        """ many oxDNA tools require specific base ordering to work properly:
            * strictly increasing staple number
            * strictly increasing base number per staple
        """
        self.logger.debug("reordering structure")

        tmp_strands = self.strands.copy()
        self.strands.clear()
        self.bases.clear()

        tmp_strand_idx = 1
        tmp_base_idx = 0

        for strand in tmp_strands.values():
            for base in strand.tour:
                base.id = tmp_base_idx
                base.strand_id = tmp_strand_idx
                self.bases[tmp_base_idx] = base
                tmp_base_idx += 1

            strand.id = tmp_strand_idx
            self.strands[tmp_strand_idx] = strand
            tmp_strand_idx += 1

    def write_structure(self, prefix: str, is_reorder=True):
        if is_reorder:
            self._reorder_structure()

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
            # NOTE: might produce unusable topology on unordered strucutre
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

    def link_strands(self, strand_id_p5: int, strand_id_p3: int, n_insert=0):
        """ connect two strands
            this operation requires full renumbering of the structure before writing
                to guarantee compatibility with oxDNA tools
            if strand ends are far apart they can be connected by adding n_insert T bases
        """
        def unit(v: np.ndarray) -> np.ndarray:
            return v / np.linalg.norm(v)

        strand1 = self.strands[strand_id_p3]
        strand2 = self.strands[strand_id_p5]

        for base in strand2.tour:
            base.strand_id = strand1.id

        last_strand1 = strand1.tour[-1]
        first_strand2 = strand2.tour[0]
        if n_insert:
            last2first = first_strand2.bb_position() - last_strand1.bb_position()
            # NOTE: flipping the T base, as 53 pairs usualy point away from each other
            bb2base_versor = -1.0 * unit(
                0.5 * (last_strand1.bb2base_versor + first_strand2.bb2base_versor))
            normversor = unit(
                0.5 * (last_strand1.normversor + first_strand2.normversor))
            for _ in range(n_insert):
                spacing = 1. / (n_insert + 1)
                base_id = max(self.bases.keys()) + 1

                bb_position = last_strand1.bb_position() + spacing * last2first
                position = bb_position_2_base_position(
                    bb_position, normversor=normversor, bb2base_versor=bb2base_versor)

                new_base = Base(
                    id=base_id,
                    seq="T",
                    strand_id=strand1.id,
                    struct=self,
                    position=position,
                    bb2base_versor=bb2base_versor,
                    normversor=normversor,
                    velocity=last_strand1.velocity,
                    angular_velocity=last_strand1.angular_velocity,
                )

                new_base.p3 = first_strand2
                first_strand2.p5 = new_base
                new_base.p5 = last_strand1
                last_strand1.p3 = new_base

                self.bases[base_id] = new_base
                strand1.tour.append(new_base)
                last_strand1 = new_base

        last_strand1.p3 = first_strand2
        first_strand2.p5 = last_strand1
        strand1.tour += strand2.tour

        self.strands.pop(strand2.id)

        return strand1.id, strand2.id

    def staple(self, connections: list):
        staple_id_chain: Dict[int, int] = dict()
        for u, v, w in connections:
            # TODO: check if safe
            while u in staple_id_chain:
                u = staple_id_chain[u]
            while v in staple_id_chain:
                v = staple_id_chain[v]

            if u == v:
                self.logger.debug("skipping circularizing connection")
                continue
            n_insert = int((w/BB_DIST) // 1.2)
            new_id, old_id = self.link_strands(u, v, n_insert)
            staple_id_chain[old_id] = new_id
