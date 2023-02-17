#!/usr/bin/env python

'''
`alife_dataframe_to_dict_of_lists` tests for
`alifedata-phyloinformatics-convert` package.
'''

import networkx as nx
import pandas as pd

import alifedata_phyloinformatics_convert as apc

def test_alife_dataframe_to_network_digraph_empty():
    df = pd.DataFrame({'id': [], 'ancestor_list': []})
    expected_output = nx.DiGraph()
    assert nx.is_isomorphic(
        apc.alife_dataframe_to_networkx_digraph(df),
        expected_output,
        node_match=True,
        edge_match=True,
    )

def test_alife_dataframe_to_network_digraph_asexual():
    df = pd.DataFrame(
        {
            'id': [0, 1, 7, 9],
            'ancestor_list': ['[None]', '[0]', '[1]', '[1]'],
            'taxon_label': ['a', 'b', 'c', 'd'],
        },
    )
    expected_output = nx.DiGraph(
        {0: [], 1: [0], 7: [1], 9: [1]}
    )

    for df_ in df, df.sample(frac=1):
        result = apc.alife_dataframe_to_networkx_digraph(df_)
        assert nx.to_dict_of_lists(expected_output) == nx.to_dict_of_lists(
            result
        )
        assert result.nodes[0]["taxon_label"] == "a"
        assert result.nodes[1]["taxon_label"] == "b"
        assert result.nodes[7]["taxon_label"] == "c"
        assert result.nodes[9]["taxon_label"] == "d"


def test_alife_dataframe_to_network_digraph_asexual_origin_time():
    df = pd.DataFrame(
        {
            'id': [0, 1, 7, 9],
            'ancestor_list': ['[None]', '[0]', '[1]', '[1]'],
            'origin_time': [0, 1, 12, 12],
        },
    )
    for df_ in df, df.sample(frac=1):
        result = apc.alife_dataframe_to_networkx_digraph(
            df_, setup_edge_lengths=True
        )
        assert result[1][0]["length"] == 1
        assert result[7][1]["length"] == 11
        assert result[9][1]["length"] == 11

        result = apc.alife_dataframe_to_networkx_digraph(
            df_, setup_edge_lengths=False
        )
        assert "length" not in result[1][0]
        assert "length" not in result[7][1]
        assert "length" not in result[9][1]


def test_alife_dataframe_to_network_digraph_asexual_edge_length():
    df = pd.DataFrame(
        {
            'id': [0, 1, 7, 9],
            'ancestor_list': ['[None]', '[0]', '[1]', '[1]'],
            'edge_length': [0, 1, 12, 12],
        },
    )
    for df_ in df, df.sample(frac=1):
        result = apc.alife_dataframe_to_networkx_digraph(
            df_, setup_edge_lengths=True
        )
        assert result[1][0]["length"] == 1
        assert result[7][1]["length"] == 12
        assert result[9][1]["length"] == 12

        result = apc.alife_dataframe_to_networkx_digraph(
            df_, setup_edge_lengths=False
        )
        assert "length" not in result[1][0]
        assert "length" not in result[7][1]
        assert "length" not in result[9][1]


def test_alife_dataframe_to_network_digraph_sexual():
    df = pd.DataFrame(
        {
            'id': [0, 1, 7, 9],
            'ancestor_list': ['[None]', '[None]', '[1,0]', '[1,7]'],
            'taxon_label': ['a', 'b', 'c', 'd'],
        }
    )

    expected_output = nx.DiGraph(
        {0: [], 1: [], 7: [1, 0], 9: [1, 7]}
    )

    for df_ in df, df.sample(frac=1):
        result = apc.alife_dataframe_to_networkx_digraph(df_)
        assert nx.to_dict_of_lists(expected_output) == nx.to_dict_of_lists(
            result
        )
        assert result.nodes[0]["taxon_label"] == "a"
        assert result.nodes[1]["taxon_label"] == "b"
        assert result.nodes[7]["taxon_label"] == "c"
        assert result.nodes[9]["taxon_label"] == "d"
