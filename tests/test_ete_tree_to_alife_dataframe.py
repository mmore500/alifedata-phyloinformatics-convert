#!/usr/bin/env python

'''
`ete_tree_to_alife_dataframe` tests for
`alifedata-phyloinformatics-convert` package.
'''

import ete3
from iterpop import iterpop as ip
from os.path import dirname, realpath

import alifedata_phyloinformatics_convert as apc


def test_newick1():

    original_tree = ete3.Tree(
        newick=f'{dirname(realpath(__file__))}/assets/pythonidae.newick',
    )

    converted_df = apc.ete_tree_to_alife_dataframe(original_tree)
    reconverted_tree = apc.alife_dataframe_to_ete_tree(converted_df)
    reconverted_df = apc.ete_tree_to_alife_dataframe(reconverted_tree)

    assert len(converted_df.compare(reconverted_df)) == 0
    assert len(converted_df) == len(reconverted_df)
    assert len(converted_df) == len([*original_tree.traverse()])
    assert len(reconverted_df) == len([*original_tree.traverse()])
    assert len(reconverted_df) == len([*reconverted_tree.traverse()])

    assert len(reconverted_tree.get_children()) == \
        len(original_tree.get_children())
    assert len(reconverted_tree.get_leaves()) \
        == len(original_tree.get_leaves())
    assert original_tree.get_ascii() == reconverted_tree.get_ascii()
    assert str(original_tree) == str(reconverted_tree)


def test_newick2():

    original_tree = ete3.Tree(
        newick=f'{dirname(realpath(__file__))}/assets/pythonidae.newick',
    )
    original_tree.standardize()

    converted_df = apc.ete_tree_to_alife_dataframe(original_tree)
    reconverted_tree = apc.alife_dataframe_to_ete_tree(converted_df)
    reconverted_df = apc.ete_tree_to_alife_dataframe(reconverted_tree)

    assert len(converted_df.compare(reconverted_df)) == 0
    assert len(converted_df) == len(reconverted_df)
    assert len(converted_df) == len([*original_tree.traverse()])
    assert len(reconverted_df) == len([*original_tree.traverse()])
    assert len(reconverted_df) == len([*reconverted_tree.traverse()])

    assert original_tree.get_ascii() == reconverted_tree.get_ascii()
    assert original_tree.robinson_foulds(reconverted_tree)[0] == 0


def test_exportattrs_iterable():

    script_directory = dirname(realpath(__file__))
    original_tree = ete3.Tree(
        newick=f'{dirname(realpath(__file__))}/assets/pythonidae.newick',
    )

    for node in original_tree.get_leaves():
        node.add_features(fish = 'Salmon')

    for node in original_tree.traverse():
        if not node.is_leaf():
            node.add_features(fish = 'Tilapia')

    converted_df = apc.ete_tree_to_alife_dataframe(
        original_tree,
        exportattrs=['fish'],
    )

    assert 'fish' in converted_df
    assert sum(converted_df['fish'] == 'Salmon') \
        == len(original_tree.get_leaves())
    assert sum(converted_df['fish'] == 'Tilapia')


def test_exportattrs_mapping():

    script_directory = dirname(realpath(__file__))
    original_tree = ete3.Tree(
        newick=f'{dirname(realpath(__file__))}/assets/pythonidae.newick',
    )
    for i, node in enumerate(original_tree.traverse()):
        if node.name is None:
            node.name = str(i)

    for node in original_tree.get_leaves():
        node.add_features(fish = 'Salmon', soup = None)

    for node in original_tree.traverse():
        if not node.is_leaf():
            node.add_features(fish = 'Tilapia', soup = None)

    converted_df = apc.ete_tree_to_alife_dataframe(
        original_tree,
        exportattrs={'fish': 'The Fish', 'soup': 'soup', 'name': 'label'},
    )

    assert 'The Fish' in converted_df
    assert sum(converted_df['The Fish'] == 'Salmon') \
        == len(original_tree.get_leaves())
    assert sum(converted_df['The Fish'] == 'Tilapia')

    assert 'soup' in converted_df
    assert ip.popsingleton(converted_df['soup'].unique()) is None

    assert 'label' in converted_df
    assert {*converted_df['label']} == {
        node.name for node in original_tree.traverse()
    }
    assert len(converted_df['label'].unique()) > 1
