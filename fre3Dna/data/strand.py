#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
from dataclasses import dataclass
from typing import Any
from typing import List

from .base import Base

""" OxDNA strand class.
    direction from 5' to 3' end
    strands must not be circular
"""


@dataclass
class Strand:
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
