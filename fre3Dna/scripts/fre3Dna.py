#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" free form DNA Origami
"""
import logging
from pathlib import Path

import click
import pandas as pd
from fre3Dna import get_version
from fre3Dna.core.attract_prep import get_nicks
from fre3Dna.core.attract_prep import get_scaffold_id
from fre3Dna.data.structure import Structure
from fre3Dna.network.graph import Graph


logger = logging.getLogger(__name__)


def print_version(ctx, _, value):
    """click print version."""
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-v",
    "--version",
    is_flag=True,
    help="Show __version__ and exit.",
    callback=print_version,
    expose_value=False,
    is_eager=True,
)
def cli():
    pass


@cli.command()
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
def prep(top, conf):
    """\b
    prep run with modified potential by creating required nick_list and scaffold number
    \b
    TOP is the name of the design file [.top]
    CONF is the scaffold strand sequence file [.oxdna, .dat]
    """
    df_top = pd.read_csv(
        top, delim_whitespace=True, skiprows=1, usecols=[0, 2, 3], names=["strand", "5p", "3p"]
    )
    df_conf = pd.read_csv(
        conf, delim_whitespace=True, skiprows=[0, 1, 2], usecols=[0, 1, 2], names=["x", "y", "z"]
    )
    scaffold_id = get_scaffold_id(df_top)
    logger.info("53_scaffold = %s", scaffold_id)
    print(f"53_nicks = {get_nicks(df_top, df_conf, scaffold_id)}")


@cli.command()
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("forces", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("cutoff", type=float)
@click.option("--o", "out", type=str, default=None, help="new file name, default to conf-name_out")
def pipeline(top, conf, forces, cutoff, out):
    """\b
    close staples perform analysis and print new files and sequences
    \b
    TOP is the name of the design file [.top]
    CONF is the scaffold strand sequence file [.oxdna, .dat]
    FORCES is the oxDNA external force file containing the basepair information []
    CUTOFF multiple of Backbone-Backbone distance to be considered for staple routing
    """
    if out is None:
        out = f"{conf.stem}_out"
    struct = Structure()
    struct.generate_from_oxDNA(top=top, conf=conf)
    struct.assign_basepairs(forces=forces)
    struct.categorise_structure()
    weighted_edges = struct.generate_connectivity(cutoff=cutoff)

    graph = Graph(struct=struct, edges=weighted_edges)
    graph.reduce_graph_reverse(reduce_iso=True, optimize_mc=False)
    staple_linkage = graph.get_routing(max_bb_multi=3.5)
    struct.staple(staple_linkage)

    struct.write_structure(out)
    struct.write_sequences(out)

    struct.structure_stats()
    graph.draw_graph()
    graph.draw_network_stats()


@cli.command()
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("forces", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("cutoff", type=float)
@click.option("--o", "out", type=str, default=None, help="new file name, default to conf-name_out")
def staple(top, conf, forces, cutoff, out):
    """\b
    close staples
    \b
    TOP is the name of the design file [.top]
    CONF is the scaffold strand sequence file [.oxdna, .dat]
    FORCES is the oxDNA external force file containing the basepair information []
    CUTOFF multiple of Backbone-Backbone distance to be considered for staple routing
    """
    if out is None:
        out = f"{conf.stem}_out"
    struct = Structure()
    struct.generate_from_oxDNA(top=top, conf=conf)
    struct.assign_basepairs(forces=forces)
    struct.categorise_structure()
    weighted_edges = struct.generate_connectivity(cutoff=cutoff)

    graph = Graph(struct=struct, edges=weighted_edges)
    graph.reduce_graph_reverse(reduce_iso=True, optimize_mc=False)
    staple_linkage = graph.get_routing(max_bb_multi=3.5)

    struct.staple(staple_linkage)
    struct.write_structure(out)


@cli.command()
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("forces", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("cutoff", type=float)
def analyse(top, conf, forces, cutoff):
    """\b
    analyse structure graph
    \b
    TOP is the name of the design file [.top]
    CONF is the scaffold strand sequence file [.oxdna, .dat]
    FORCES is the oxDNA external force file containing the basepair information []
    CUTOFF multiple of Backbone-Backbone distance to be considered for staple routing
    """
    struct = Structure()
    struct.generate_from_oxDNA(top=top, conf=conf)
    struct.assign_basepairs(forces=forces)
    struct.categorise_structure()
    struct.structure_stats()

    weighted_edges = struct.generate_connectivity(cutoff=cutoff)
    graph = Graph(struct=struct, edges=weighted_edges)
    graph.draw_graph()
    graph.draw_network_stats()


@cli.command()
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option("--o", "out", type=str, default=None, help="sequence file name, default to conf name")
def sequence(top, conf, out):
    """\b
    analyse structure graph
    \b
    TOP is the name of the design file [.top]
    CONF is the scaffold strand sequence file [.oxdna, .dat]
    """
    if out is None:
        out = conf.stem
    struct = Structure()
    struct.generate_from_oxDNA(top=top, conf=conf)
    struct.write_sequences(out)
