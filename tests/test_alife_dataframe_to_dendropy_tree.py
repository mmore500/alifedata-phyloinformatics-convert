#!/usr/bin/env python

'''
`alife_dataframe_to_dendropy_tree` tests for
`alifedata-phyloinformatics-convert` package.
'''

from os.path import dirname, realpath
import pandas as pd

import alifedata_phyloinformatics_convert as apc


def _setify(df):
    return set(
        (
            int(row['id']),
            row['ancestor_list'].upper(),
            int(row['origin_time']),
        )
        for __, row in df.iterrows()
    )


def _setify_no_origin_time(df):
    return set(
        (
            int(row['id']),
            row['ancestor_list'].upper(),
        )
        for __, row in df.iterrows()
    )


def test():
    original_df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata.csv',
    )
    converted_tree = apc.alife_dataframe_to_dendropy_tree(original_df)
    reconverted_df = apc.dendropy_tree_to_alife_dataframe(converted_tree)

    assert _setify(original_df) == _setify(reconverted_df)


def test_minimal():
    original_df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata_minimal.csv'
    )
    converted_tree = apc.alife_dataframe_to_dendropy_tree(original_df)
    reconverted_df = apc.dendropy_tree_to_alife_dataframe(converted_tree)

    assert _setify_no_origin_time(original_df) \
        == _setify_no_origin_time(reconverted_df)


def test_empty():
    df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata_empty.csv'
    )
    assert apc.alife_dataframe_to_dendropy_tree(df) is None


def test_minimal_empty():
    df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata_minimal_empty.csv'
    )
    assert apc.alife_dataframe_to_dendropy_tree(df) is None


def test_setattrs_iterable():
    original_df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata.csv',
    )
    converted_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df,
        setattrs=['num_orgs'],
    )
    for node in converted_tree:
        assert hasattr(node, 'num_orgs')
        assert isinstance(node.num_orgs, int)
        assert node.num_orgs >= 0


def test_setattrs_mapping():
    original_df = pd.read_csv(
        f'{dirname(realpath(__file__))}/assets/alifedata.csv',
    )
    converted_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df,
        setattrs={
            'num_orgs': 'num_orgs_',
            'tot_orgs': 'tot_orgs_',
        },
    )
    for node in converted_tree:
        assert hasattr(node, 'num_orgs_')
        assert not hasattr(node, 'num_orgs')
        assert hasattr(node, 'tot_orgs_')
        assert not hasattr(node, 'tot_orgs')
        assert isinstance(node.num_orgs_, int)
        assert isinstance(node.tot_orgs_, int)
        assert node.num_orgs_ >= 0
        assert node.tot_orgs_ >= 0
