#!/usr/bin/env python

'''
`dendropy_tree_to_scipy_linkage_matrix` tests for
`alifedata-phyloinformatics-convert` package.
'''

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
    reconstructed_linkage_matrix \
        = apc.dendropy_tree_to_scipy_linkage_matrix(tree)

    assert set(map(tuple, linkage_matrix)) \
        == set(map(tuple, reconstructed_linkage_matrix))


def test_no_cluster_ids():
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
    for node in tree:
        del node.cluster_id
    reconstructed_linkage_matrix \
        = apc.dendropy_tree_to_scipy_linkage_matrix(tree)

    assert [*linkage_matrix[:, 2]] == [*reconstructed_linkage_matrix[:, 2]]
    assert [*linkage_matrix[:, 3]] == [*reconstructed_linkage_matrix[:, 3]]
