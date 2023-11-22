from Bio import Phylo
from collections import defaultdict
from collections.abc import Mapping
import functools
from iterpop import iterpop as ip
import itertools as it
from nanto import isanan
import networkx as nx
import numbers
import numpy as np
import opytional as opyt
import pandas as pd
import typing

from ._impl import rgetattr as _rgetattr


def networkx_digraph_to_alife_dataframe(
    graph: nx.Graph,
) -> pd.DataFrame:
    """Convert a networkx digraph to a dataframe formatted to the artificial
    life community data format standards.

    The following node attributes will automatically be exported to
    dataframe columns, if available:
        * branch_length,
        * edge_length,
        * length,
        * weight,
        * label,
        * name,
        * origin_time, and
        * taxon_label.

    Parameters
    ----------
    graph:
        networkx graph to convert.
    """

    if not all(
        isinstance(x, numbers.Integral) and x >= 0 for x in graph.nodes
    ):
        assert not any("label" in data for __, data in graph.nodes(data=True))
        graph = nx.convert_node_labels_to_integers(
            graph,
            label_attribute="label",
        )

    if len(graph.nodes) == 0:
        return pd.DataFrame({"id": [], "ancestor_list": []})

    records = []
    visited_nodes = set(graph.nodes)
    for from_, to, edge_data in it.chain(
        graph.edges(data=True),
        (
            (ip.popsingleton(visited_nodes), None, {})
            for __ in range(1)
        ),
    ):
        attrs = graph.nodes[from_]
        def try_get_attr(attr):
            return opyt.or_else(
                attrs.get(attr, None),
                lambda: edge_data.get(attr, np.nan),
            )

        record = {
            "id": from_,
            "ancestor_list" : str([to]),
            "ancestor_id" : opyt.or_value(to, from_),
            "branch_length" : try_get_attr("branch_length"),
            "edge_length" : try_get_attr("edge_length"),
            "length" : try_get_attr("length"),
            "weight" : try_get_attr("weight"),
            "label" : try_get_attr("label"),
            "name" : try_get_attr("name"),
            "origin_time" : try_get_attr("origin_time"),
            "taxon_label" : try_get_attr("taxon_label"),
        }
        visited_nodes.remove(from_)
        records.append(record)

    return pd.DataFrame.from_records(records).dropna(axis=1, how='all')
