#!/usr/bin/env python

'''
`scipy_linkage_matrix_to_dendropy_tree` tests for
`alifedata-phyloinformatics-convert` package.
'''

import dendropy
from dendropy.calculate import treecompare
import numpy as np

import alifedata_phyloinformatics_convert as apc


def test():

    # from https://stackoverflow.com/a/40983611/17332200
    linkage_matrix = np.array([
        [7.0, 9.0, 0.3, 2.0],
        [4.0, 6.0, 0.5, 2.0],
        [5.0, 12.0, 0.5, 3.0],
        [2.0, 13.0, 0.53851648, 4.0],
        [3.0, 14.0, 0.58309519, 5.0],
        [1.0, 15.0, 0.64031242, 6.0],
        [10.0, 11.0, 0.72801099, 3.0],
        [8.0, 17.0, 1.2083046, 4.0],
        [0.0, 16.0, 1.5132746, 7.0],
        [18.0, 19.0, 1.92353841, 11.0],
    ])

    tree = apc.scipy_linkage_matrix_to_dendropy_tree(linkage_matrix)

    assert sum(1 for __ in tree.leaf_node_iter()) == 11

    def node(cluster_id: float) -> dendropy.Node:
        return next(
            node for node in tree if node.cluster_id == cluster_id
        )

    assert {n.cluster_id for n in node(11).child_nodes()} == {7., 9.}
    assert node(7).edge_length == 0.3 / 2
    assert node(9).edge_length == 0.3 / 2

    assert {n.cluster_id for n in node(12).child_nodes()} == {4., 6.}
    assert node(4).edge_length == 0.5 / 2
    assert node(6).edge_length == 0.5 / 2

    assert {n.cluster_id for n in node(17).child_nodes()} == {10., 11.}
    assert node(10).edge_length == 0.72801099 / 2
    assert node(11).edge_length == 0.72801099 / 2

    assert {n.cluster_id for n in node(20).child_nodes()} == {18., 19.}
    assert node(18).edge_length == 1.92353841 / 2
    assert node(19).edge_length == 1.92353841 / 2
    assert node(20).parent_node is None


def test_leaf_taxon_labels():

    # from https://stackoverflow.com/a/40983611/17332200
    linkage_matrix = np.array([
        [7.0, 9.0, 0.3, 2.0],
        [4.0, 6.0, 0.5, 2.0],
        [5.0, 12.0, 0.5, 3.0],
        [2.0, 13.0, 0.53851648, 4.0],
        [3.0, 14.0, 0.58309519, 5.0],
        [1.0, 15.0, 0.64031242, 6.0],
        [10.0, 11.0, 0.72801099, 3.0],
        [8.0, 17.0, 1.2083046, 4.0],
        [0.0, 16.0, 1.5132746, 7.0],
        [18.0, 19.0, 1.92353841, 11.0],
    ])

    tree1 = apc.scipy_linkage_matrix_to_dendropy_tree(
        linkage_matrix,
        [*range(11)],
    )

    for node in tree1.leaf_node_iter():
        # check that all nodes have a corresponding taxon
        assert node.taxon is not None
        # check that the taxons have the right labels
        assert node.taxon.label is not None
        assert node.cluster_id == node.taxon.label
        # check that the nodes themselves have the right labels
        assert node.label is not None
        assert node.cluster_id == node.label

    for node in tree1.preorder_internal_node_iter():
        assert node.taxon is None

    tree2 = apc.scipy_linkage_matrix_to_dendropy_tree(
        linkage_matrix,
        [*range(11)],
    )

    tree2.migrate_taxon_namespace(tree1.taxon_namespace)

    assert treecompare.symmetric_difference(tree1, tree2) == 0

    tree3 = apc.scipy_linkage_matrix_to_dendropy_tree(
        linkage_matrix,
        leaf_labels=[*range(11)],
        label_leaf_nodes=False,
    )

    for node in tree3.leaf_node_iter():
        # check that all nodes have a corresponding taxon
        assert node.taxon is not None
        # check that the taxons have the right labels
        assert node.taxon.label is not None
        assert node.cluster_id == node.taxon.label
        # check that the nodes themselves have the right labels
        assert node.label is None

    tree4 = apc.scipy_linkage_matrix_to_dendropy_tree(
        linkage_matrix,
        leaf_labels=[*range(11)],
        label_leaf_taxa=False,
    )

    for node in tree4.leaf_node_iter():
        # check that all nodes don't have corresponding taxon labels
        assert not hasattr(node.taxon, "label")
        # check that the nodes themselves have the right labels
        assert node.label is not None
        assert node.cluster_id == node.label
