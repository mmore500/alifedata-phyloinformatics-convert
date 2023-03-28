#!/usr/bin/env python

'''
`dendropy_tree_to_alife_dataframe` tests for
`alifedata-phyloinformatics-convert` package.
'''

import dendropy
from dendropy.calculate import treecompare
from iterpop import iterpop as ip
import networkx as nx
import numpy as np
from os.path import dirname, realpath

import alifedata_phyloinformatics_convert as apc


def test_newick1():

    source_tree = dendropy.Tree.get(
        path=f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        schema='newick',
    )

    source_df = apc.dendropy_tree_to_alife_dataframe(source_tree).dropna(
        axis=1, how='all'
    )
    networkx_tree = apc.alife_dataframe_to_networkx_digraph(source_df)

    reconverted_df = apc.networkx_digraph_to_alife_dataframe(networkx_tree)

    assert reconverted_df.drop(
        "ancestor_id", axis=1
    ).sort_values(by="id").sort_index(axis=1).reset_index(drop=True).equals(
        source_df.fillna(value=np.nan).sort_values(by="id").sort_index(axis=1).reset_index(
            drop=True
        )
    )


def test_newick2():

    source_tree = dendropy.Tree.get(
        path=f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        schema='newick',
    )
    source_tree.resolve_polytomies()
    source_tree.suppress_unifurcations()
    source_tree.collapse_basal_bifurcation()

    source_df = apc.dendropy_tree_to_alife_dataframe(source_tree).dropna(
        axis=1, how='all'
    )
    networkx_tree = apc.alife_dataframe_to_networkx_digraph(source_df)

    reconverted_df = apc.networkx_digraph_to_alife_dataframe(networkx_tree)

    assert reconverted_df.drop(
        "ancestor_id", axis=1
    ).sort_values(by="id").sort_index(axis=1).reset_index(drop=True).equals(
        source_df.fillna(value=np.nan).sort_values(
            by="id"
        ).sort_index(axis=1).reset_index(
            drop=True
        )
    )


def test_nexml():

    script_directory = dirname(realpath(__file__))
    source_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    source_df = apc.dendropy_tree_to_alife_dataframe(source_tree).dropna(
        axis=1, how='all'
    )

    networkx_tree = apc.alife_dataframe_to_networkx_digraph(source_df)
    reconverted_df = apc.networkx_digraph_to_alife_dataframe(networkx_tree)

    assert reconverted_df.drop(
        "ancestor_id", axis=1
    ).sort_values(by="id").sort_index(axis=1).reset_index(drop=True).equals(
        source_df.fillna(value=np.nan).sort_values(
            by="id"
        ).sort_index(axis=1).reset_index(
            drop=True
        )
    )


def test_relabel_ints():

    script_directory = dirname(realpath(__file__))
    source_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    source_df = apc.dendropy_tree_to_alife_dataframe(source_tree).dropna(
        axis=1, how='all'
    )

    networkx_tree = apc.alife_dataframe_to_networkx_digraph(source_df)
    nx.relabel_nodes(
        networkx_tree, {n : f"{n}!" for n in networkx_tree.nodes}, copy=False
    )
    reconverted_df = apc.networkx_digraph_to_alife_dataframe(networkx_tree)

    source_df["label"] = source_df["id"].apply(lambda x: f"{x}!")

    assert reconverted_df.drop(
        "ancestor_id", axis=1
    ).sort_values(by="id").sort_index(axis=1).reset_index(drop=True).equals(
        source_df.fillna(value=np.nan).sort_values(
            by="id"
        ).sort_index(axis=1).reset_index(
            drop=True
        )
    )
