import networkx as nx
import pandas as pd
import typing

from .alife_dataframe_to_dict_of_lists import alife_dataframe_to_dict_of_lists


# see also https://github.com/alife-data-standards/alife-std-dev-python/blob/f21a63d70077345441b9b52c3f470f5dca1127c1/tests/test_phylogeny_loader.py
def alife_dataframe_to_networkx_digraph(
    df: pd.DataFrame,
    setup_edge_lengths: bool = False,
) -> nx.DiGraph:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as a networkx directed graph.

    Directed edges point from child to parent. Clades that do not share a
    common ancestor are supported. Call `.reverse()` to orient so directed
    edges point from parent to child.

    If enabled, branch lengths will be set up based on the origin_time
    attribute.

    The following column values will automatically be applied as node attributes, if available:
        * branch_length,
        * edge_length,
        * length,
        * weight,
        * label,
        * name,
        * origin_time, and
        * taxon_label.

    Will raise ValueError if nonequivalent branch_length and edge_length
    columns are provided.

    Parameters
    ----------
    df:
        Pandas DataFrame to convert.
    setup_edge_lengths: bool, optional
        Should we try to set up edge lengths using the origin_time column?
        Will not override if branch_length or edge_length is provided as a
        column of df.
    """

    if "branch_length" in df and "edge_length" in df and not (
        df["branch_length"].equals(df["edge_length"])
    ):
        raise ValueError

    g = nx.from_dict_of_lists(
        alife_dataframe_to_dict_of_lists(df),
        create_using=nx.DiGraph,
    )

    nx.set_node_attributes(
        g,
        df.set_index("id")[[
            attr
            for attr in (
                "branch_length",
                "edge_length",
                "length",
                "weight",
                "label",
                "name",
                "origin_time",
                "taxon_label",
            )
            if attr in df
        ]].to_dict(orient="index"),
    )

    if setup_edge_lengths and "length" in df:
        for from_, to in g.edges:
            g[from_][to]["length"] = g.nodes[from_]["length"]
    elif setup_edge_lengths and "edge_length" in df:
        for from_, to in g.edges:
            g[from_][to]["length"] = g.nodes[from_]["edge_length"]
    elif setup_edge_lengths and "branch_length" in df:
        for from_, to in g.edges:
            g[from_][to]["length"] = g.nodes[from_]["branch_length"]
    elif setup_edge_lengths and "origin_time" in df:
        for from_, to in g.edges:
            g[from_][to]["length"] = (
                g.nodes[from_]["origin_time"] - g.nodes[to]["origin_time"]
            )
    elif setup_edge_lengths:
        for from_, to in g.edges:
            g[from_][to]["length"] = 1

    return g
