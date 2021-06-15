#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os

from pathlib import Path
from shutil import copyfile, copytree

import click
from fre3Dna.version import get_version


""" cascaded mrDNA-driven MD flexible fitting:
"""
logger = logging.getLogger(__name__)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def cli():
    pass


@cli.command()
def main(cadnano, mrc, sequence, gpu, prefix, multidomain):
    """
    """
    return
