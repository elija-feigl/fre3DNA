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

    def __repr__(self) -> str:
        return f"Base {self.id}, strand {self.strand_id}, seq {self.seq}"
