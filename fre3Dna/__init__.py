#!/usr/bin/env python
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
import logging
import os
from pathlib import Path


def get_version() -> str:
    return __version__


def get_resource(resources: str) -> Path:
    return Path(__file__).parent / "resources" / resources


def _init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s | [%(module)s]\t%(levelname)s\t- %(message)s", "%Y.%m.%d %H:%M"
    )
    handler.setFormatter(formatter)

    if os.name != "nt":  # no ANSI escape support on Windows
        logging.addLevelName(
            logging.DEBUG, "\033[1;34m%s\033[1;0m" % logging.getLevelName(logging.DEBUG)
        )
        logging.addLevelName(
            logging.WARNING, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.WARNING)
        )
        logging.addLevelName(
            logging.ERROR, "\033[1;33m%s\033[1;0m" % logging.getLevelName(logging.ERROR)
        )
        logging.addLevelName(
            logging.CRITICAL, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.CRITICAL)
        )
    logger.addHandler(handler)


_init_logging()


version_info = [0, 7, 0, "dev0"]

__version__ = ".".join([str(sub) for sub in version_info])
__all__ = ["__version__"]
