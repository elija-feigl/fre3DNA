#! /usr/bin/env python
# -*- coding: utf-8 -*-


from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

""" OxDNA base class.
"""


@dataclass
class Base(object):
    id: int
    seq: str
    strand_id: Optional[int]
    struct: Any
    position: np.ndarray
    bb2base_versor: np.ndarray
    normversor: np.ndarray
    velocity: np.ndarray
    angular_velocity: np.ndarray

    _p5_id: Optional[int]
    _p3_id: Optional[int]

    p5: Optional["Base"] = None
    p3: Optional["Base"] = None
    across: Optional["Base"] = None

    def bb_position(self) -> np.ndarray:
        """ backbone position for oxDNA2 (w. grooving)
            oxDNA1: self.position - 0.40 * self.bb2base_versor
        """
        thrid_versor = np.cross(self.normversor, self.bb2base_versor)
        return self.position - 0.34 * self.bb2base_versor + 0.3408 * thrid_versor

    def base_position(self) -> np.ndarray:
        return self.position + 0.40 * self.bb2base_versor

    def stack_position(self) -> np.ndarray:
        return self.position + 0.34 * self.bb2base_versor

    def __repr__(self) -> str:
        return f"Base {self.id}, strand {self.strand_id}, seq {self.seq}"
