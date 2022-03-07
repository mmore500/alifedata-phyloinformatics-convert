#!/usr/bin/env python

'''
`dendropy_tree_to_alife_dataframe` tests for
`alifedata-phyloinformatics-convert` package.
'''

from Bio import Phylo
import dendropy
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
    # print(dendropy_tree)
    dendropy_df = apc.dendropy_tree_to_alife_dataframe(dendropy_tree)

    biopython_tree = Phylo.read(
        f'{script_directory}/assets/pythonidae.annotated.nexml',
        'nexml',
    )
    # print(biopython_tree)
    biopython_df = apc.biopython_tree_to_alife_dataframe(biopython_tree)
    assert len(dendropy_df) == len(biopython_df)

    assert len(dendropy_df) == len(biopython_df)
