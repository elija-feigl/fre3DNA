#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging

import click
import pandas as pd
from pathlib import Path

from fre3Dna.version import get_version
from fre3Dna.core.attract_prep import get_scaffold_id, get_nicks
from fre3Dna.data.structure import Structure

""" free form DNA Origami
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
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
def prep(top, conf):
    """ prep run with modified potential by creating required nick_list and scaffold number

        TOP is the name of the design file [.top]\n
        CONF is the scaffold strand sequence file [.oxdna, .dat]\n
    """
    df_top = pd.read_csv(
        top, delim_whitespace=True,
        skiprows=1, usecols=[0, 2, 3], names=["strand", "5p", "3p"]
    )
    df_conf = pd.read_csv(
        conf, delim_whitespace=True,
        skiprows=[0, 1, 2], usecols=[0, 1, 2], names=["x", "y", "z"]
    )
    scaffold_id = get_scaffold_id(df_top)
    logger.info(f"53_scaffold = {scaffold_id}")
    print(f"53_nicks = {get_nicks(df_top, df_conf, scaffold_id)}")


@cli.command()
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
def staple(top, conf):
    """ close staples

        TOP is the name of the design file [.top]\n
        CONF is the scaffold strand sequence file [.oxdna, .dat]\n
    """
    struct = Structure()
    struct.generate_from_oxDNA(top=Path(top), conf=Path(conf))
