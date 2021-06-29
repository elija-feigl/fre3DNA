#! /usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from .base import Base

from typing import Any, List

""" OxDNA strand class.
    direction from 5' to 3' end
    strands must not be circular
"""


@dataclass
class Strand(object):
    id: int
    struct: Any
    tour: List[Base]
    is_scaffold: bool

    def __repr__(self) -> str:
        return f"Strand {self.id}, length {len(self)}"

    def p5(self) -> Base:
        return self.tour[0]

    def p3(self) -> Base:
        return self.tour[-1]

    def __len__(self) -> int:
        return len(self.tour)
