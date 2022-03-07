from collections import defaultdict
import dendropy
from nanto import nantonone
import numpy as np
import typing


def scipy_linkage_matrix_to_dendropy_tree(
    matrix: np.array,
    leaf_taxon_labels: typing.Optional[typing.List] = None,
) -> dendropy.Tree:
    """
    Parameters
    ----------
    leaf_taxon_labels: list, optional
        Taxon labels for leaf nodes of tree. If provided, must exactly equal
        the number of leaf nodes (i.e., num matrix rows + 1).
    """

    # scipy linkage format
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
    num_rows = len(matrix)

    # cluster id -> node
    nodes = defaultdict(dendropy.Node)

    if leaf_taxon_labels is not None:
        assert len(leaf_taxon_labels) == num_rows + 1
        for label, cluster_id in enumerate(leaf_taxon_labels):
            nodes[cluster_id].taxon = dendropy.Taxon(label=label)

    for row_idx, row in enumerate(matrix):
        parent_cluster = row_idx + num_rows + 1
        joined_cluster1, joined_cluster2, cluster_distance, cluster_size = row

        nodes[parent_cluster].cluster_id = parent_cluster
        for child_cluster in joined_cluster1, joined_cluster2:
            nodes[parent_cluster].add_child(nodes[child_cluster])
            nodes[child_cluster].cluster_id = child_cluster
            nodes[child_cluster].edge_length = nantonone(cluster_distance / 2)

    return dendropy.Tree(seed_node=nodes[parent_cluster])
