#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging

import click
import pandas as pd
from pathlib import Path


from fre3Dna.version import get_version
from fre3Dna.core.attract_prep import get_scaffold_id, get_nicks
from fre3Dna.data.structure import Structure
from fre3Dna.data.graph import Graph

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
@click.argument('forces', type=click.Path(exists=True))
@click.argument('cutoff', type=float)
def staple(top, conf, forces, cutoff):
    """ close staples

        TOP is the name of the design file [.top]\n
        CONF is the scaffold strand sequence file [.oxdna, .dat]\n
        FORCES is the oxDNA external force file containing the basepair information []\n
        CUTOFF multiple of Backbone-Backbone distance to be considered for staple routing\n
    """
    struct = Structure()
    struct.generate_from_oxDNA(top=Path(top), conf=Path(conf))
    struct.assign_basepairs(forces=Path(forces))
    struct.categorise_structure()
    weighted_edges = struct.generate_connectivity(cutoff=cutoff)

    graph = Graph(struct=struct, edges=weighted_edges)

    graph.draw_network_stats()
    graph.reduce_graph_simple()
    struct.staple(graph.get_routing(max_bb_multi=3.3))
    struct.structure_stats()

    struct.write_structure("test")

    graph.draw_graph()
    graph.draw_network_stats()
