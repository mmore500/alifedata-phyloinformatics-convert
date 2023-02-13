#!/usr/bin/env python

'''
`alife_dataframe_to_dict_of_lists` tests for
`alifedata-phyloinformatics-convert` package.
'''

import pandas as pd

import alifedata_phyloinformatics_convert as apc

def test_alife_dataframe_to_dict_of_lists_empty():
    df = pd.DataFrame({'id': [], 'ancestor_list': []})
    expected_output = {}
    assert apc.alife_dataframe_to_dict_of_lists(df) == expected_output

def test_alife_dataframe_to_dict_of_lists_asexual():
    df = pd.DataFrame(
        {'id': [0, 1, 7, 9], 'ancestor_list': ['[None]', '[0]', '[1]', '[1]']}
    )
    expected_output = {0: [], 1: [0], 7: [1], 9: [1]}
    assert apc.alife_dataframe_to_dict_of_lists(df) == expected_output
    assert apc.alife_dataframe_to_dict_of_lists(
        df.sample(frac=1)
    ) == expected_output

def test_alife_dataframe_to_dict_of_lists_sexual():
    df = pd.DataFrame(
        {
            'id': [0, 1, 7, 9],
            'ancestor_list': ['[None]', '[None]', '[1,0]', '[1,7]']
        }
    )

    expected_output = {0: [], 1: [], 7: [1, 0], 9: [1, 7]}
    assert apc.alife_dataframe_to_dict_of_lists(df) == expected_output
    assert apc.alife_dataframe_to_dict_of_lists(
        df.sample(frac=1)
    ) == expected_output
