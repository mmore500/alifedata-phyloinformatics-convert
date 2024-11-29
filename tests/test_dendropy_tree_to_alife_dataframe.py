#!/usr/bin/env python

'''
`dendropy_tree_to_alife_dataframe` tests for
`alifedata-phyloinformatics-convert` package.
'''

import dendropy
from dendropy.calculate import treecompare
from iterpop import iterpop as ip
from os.path import dirname, realpath

import alifedata_phyloinformatics_convert as apc


def test_newick1():

    original_tree = dendropy.Tree.get(
        path=f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        schema='newick',
    )
    original_tree.is_rooted = True

    converted_df = apc.dendropy_tree_to_alife_dataframe(original_tree)
    reconverted_tree = apc.alife_dataframe_to_dendropy_tree(converted_df)
    reconverted_df = apc.dendropy_tree_to_alife_dataframe(reconverted_tree)

    assert len(converted_df.compare(reconverted_df)) == 0
    assert len(converted_df) == len(reconverted_df)
    assert len(converted_df) == len(original_tree.nodes())
    assert len(reconverted_df) == len(original_tree.nodes())
    assert len(reconverted_df) == len(reconverted_tree.nodes())

    assert len(reconverted_tree.nodes()) == len(original_tree.nodes())
    assert len(reconverted_tree.leaf_nodes()) \
        == len(original_tree.leaf_nodes())
    assert original_tree.as_ascii_plot() == reconverted_tree.as_ascii_plot()
    assert str(original_tree).strip() == reconverted_tree.as_string('newick').strip()


def test_newick2():

    original_tree = dendropy.Tree.get(
        path=f'{dirname(realpath(__file__))}/assets/APG_Angiosperms.newick',
        schema='newick',
    )
    original_tree.resolve_polytomies()
    original_tree.suppress_unifurcations()
    original_tree.collapse_basal_bifurcation()
    original_tree.is_rooted = True

    converted_df = apc.dendropy_tree_to_alife_dataframe(original_tree)
    reconverted_tree = apc.alife_dataframe_to_dendropy_tree(converted_df)
    reconverted_df = apc.dendropy_tree_to_alife_dataframe(reconverted_tree)

    assert len(converted_df.compare(reconverted_df)) == 0
    assert len(converted_df) == len(reconverted_df)
    assert len(converted_df) == len(original_tree.nodes())
    assert len(reconverted_df) == len(original_tree.nodes())
    assert len(reconverted_df) == len(reconverted_tree.nodes())

    reconverted_tree.migrate_taxon_namespace(original_tree.taxon_namespace)

    assert original_tree.as_ascii_plot() == reconverted_tree.as_ascii_plot()
    assert str(original_tree).strip() == str(reconverted_tree).strip()
    assert treecompare.symmetric_difference(
        original_tree,
        reconverted_tree,
    ) == 0


def test_nexml():

    script_directory = dirname(realpath(__file__))
    original_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    original_tree.is_rooted = True
    converted_df = apc.dendropy_tree_to_alife_dataframe(original_tree)
    reconverted_tree = apc.alife_dataframe_to_dendropy_tree(
        converted_df,
    )

    # ensure taxons can be compared between trees
    reconverted_tree.migrate_taxon_namespace(original_tree.taxon_namespace)

    assert treecompare.symmetric_difference(
        original_tree,
        reconverted_tree,
    ) == 0


def test_exportattrs_iterable():

    script_directory = dirname(realpath(__file__))
    original_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    original_tree.is_rooted = True

    for node in original_tree.leaf_node_iter():
        node.fish = 'Salmon'

    for node in original_tree.preorder_internal_node_iter():
        node.fish = 'Tilapia'

    converted_df = apc.dendropy_tree_to_alife_dataframe(
        original_tree,
        exportattrs=['fish'],
    )

    assert 'fish' in converted_df
    assert sum(converted_df['fish'] == 'Salmon') \
        == len(original_tree.leaf_nodes())
    assert sum(converted_df['fish'] == 'Tilapia') \
        == len(original_tree.internal_nodes())


def test_exportattrs_mapping():

    script_directory = dirname(realpath(__file__))
    original_tree = dendropy.Tree.get(
        path=f'{script_directory}/assets/pythonidae.annotated.nexml',
        schema='nexml',
    )
    original_tree.is_rooted = True
    for i, node in enumerate(original_tree):
        if node.taxon is None:
            node.taxon = original_tree.taxon_namespace.new_taxon(str(i))

    for node in original_tree.leaf_node_iter():
        node.fish = 'Salmon'
        node.soup = None

    for node in original_tree.preorder_internal_node_iter():
        node.fish = 'Tilapia'
        node.soup = None

    converted_df = apc.dendropy_tree_to_alife_dataframe(
        original_tree,
        exportattrs={'fish': 'The Fish', 'soup': 'soup', 'taxon.label': 'name'},
    )

    assert 'The Fish' in converted_df
    assert sum(converted_df['The Fish'] == 'Salmon') \
        == len(original_tree.leaf_nodes())
    assert sum(converted_df['The Fish'] == 'Tilapia') \
        == len(original_tree.internal_nodes())

    assert 'soup' in converted_df
    assert ip.popsingleton(converted_df['soup'].unique()) is None

    assert 'name' in converted_df
    assert {*converted_df['name']} ==  {node.taxon.label for node in original_tree}
    assert len(converted_df['name'].unique()) > 1
