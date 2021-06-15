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

    def __repr__(self) -> str:
        return f"Strand {self.id}, length {self.length()}"

    def p5(self) -> Base:
        return self.tour[0]

    def p3(self) -> Base:
        return self.tour[-1]

    def length(self) -> int:
        return len(self.tour)
