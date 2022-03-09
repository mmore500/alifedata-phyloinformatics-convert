#!/usr/bin/env python

'''
`dendropy_tree_to_alife_dataframe` tests for
`alifedata-phyloinformatics-convert` package.
'''

from Bio import Phylo
import dendropy
from iterpop import iterpop as ip
from os.path import dirname, realpath

import alifedata_phyloinformatics_convert as apc


def test_newick():

    dendropy_tree = dendropy.Tree.get(
        path=f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        schema='newick',
    )
    dendropy_df = apc.dendropy_tree_to_alife_dataframe(dendropy_tree)

    biopython_tree = Phylo.read(
        f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        'newick',
    )
    biopython_df = apc.biopython_tree_to_alife_dataframe(biopython_tree)

    assert len(dendropy_df) == len(biopython_df)

    assert len(dendropy_df['id'].compare(biopython_df['id'])) == 0
    assert len(dendropy_df['ancestor_list'].compare(
        biopython_df['ancestor_list']),
    ) == 0
    assert len(dendropy_df['origin_time'].compare(
        biopython_df['origin_time']),
    ) == 0
    assert len(dendropy_df['edge_length'].compare(
        biopython_df['branch_length']),
    ) == 0


def test_nexml():

    script_directory = dirname(realpath(__file__))
    dendropy_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    dendropy_df = apc.dendropy_tree_to_alife_dataframe(dendropy_tree)

    biopython_tree = Phylo.read(
        f'{script_directory}/assets/pythonidae.annotated.nexml',
        'nexml',
    )
    biopython_df = apc.biopython_tree_to_alife_dataframe(biopython_tree)
    assert len(dendropy_df) == len(biopython_df)

    assert len(dendropy_df) == len(biopython_df)


def test_exportattrs_iterable():

    original_tree = Phylo.read(
        f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        'newick',
    )

    for node in original_tree.get_terminals():
        node.fish = 'Salmon'

    for node in original_tree.get_nonterminals():
        node.fish = 'Tilapia'

    converted_df = apc.biopython_tree_to_alife_dataframe(
        original_tree,
        exportattrs=['fish'],
    )

    assert 'fish' in converted_df
    assert sum(converted_df['fish'] == 'Salmon') \
        == len(original_tree.get_terminals())
    assert sum(converted_df['fish'] == 'Tilapia') \
        == len(original_tree.get_nonterminals())


def test_exportattrs_mapping():

    original_tree = Phylo.read(
        f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        'newick',
    )

    for clade in original_tree.get_terminals():
        clade.fish = 'Salmon'
        clade.soup = None

    for clade in original_tree.get_nonterminals():
        clade.fish = 'Tilapia'
        clade.soup = None

    converted_df = apc.biopython_tree_to_alife_dataframe(
        original_tree,
        exportattrs={'fish': 'The Fish', 'soup': 'soup'},
    )

    assert 'The Fish' in converted_df
    assert sum(converted_df['The Fish'] == 'Salmon') \
        == len(original_tree.get_terminals())
    assert sum(converted_df['The Fish'] == 'Tilapia') \
        == len(original_tree.get_nonterminals())

    assert 'soup' in converted_df
    assert ip.popsingleton(converted_df['soup'].unique()) is None
