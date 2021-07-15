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

    _p5_id: Optional[int] = None
    _p3_id: Optional[int] = None

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

    def conf(self) -> str:
        pos = ' '.join(map(str, self.position))
        v_1 = ' '.join(map(str, self.bb2base_versor))
        v_3 = ' '.join(map(str, self.normversor))
        vel = ' '.join(map(str, self.velocity))
        a_vel = ' '.join(map(str, self.angular_velocity))
        return f"{pos} {v_1} {v_3} {vel} {a_vel}"

    def top(self) -> str:
        p3_id = "-1" if self.p3 is None else self.p3.id
        p5_id = "-1" if self.p5 is None else self.p5.id
        return f"{self.strand_id} {self.seq} {p3_id} {p5_id}"

    def __repr__(self) -> str:
        return f"{self.id} {self.top()}"

    def __str__(self) -> str:
        return f"{self.id}: {self.top()}: {self.conf()}"
