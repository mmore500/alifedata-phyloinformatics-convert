from collections import defaultdict
import dendropy
from nanto import nantonone
import numpy as np
import typing


def scipy_linkage_matrix_to_dendropy_tree(
    matrix: np.array,
    leaf_labels: typing.Optional[typing.List] = None,
    label_leaf_nodes: bool = True,
    label_leaf_taxa: bool = True,
) -> dendropy.Tree:
    """
    Parameters
    ----------
    leaf_labels: list, optional
        Labels for leaves of tree. If provided, must exactly equal the number
        of leaf nodes (i.e., num matrix rows + 1).
    label_leaf_nodes: bool, default True
        Apply leaf labels directly to leaf nodes.
    label_leaf_taxa: bool, default True
        Apply leaf labels to leaf node taxa.
    """

    # scipy linkage format
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
    num_rows = len(matrix)

    # cluster id -> node
    nodes = defaultdict(dendropy.Node)

    if leaf_labels is not None:
        assert len(leaf_labels) == num_rows + 1
        for cluster_id, label in enumerate(leaf_labels):
            if label_leaf_taxa:
                nodes[cluster_id].taxon = dendropy.Taxon(label=label)
            if label_leaf_nodes:
                nodes[cluster_id].label = label

    for row_idx, row in enumerate(matrix):
        parent_cluster = row_idx + num_rows + 1
        joined_cluster1, joined_cluster2, cluster_distance, cluster_size = row

        nodes[parent_cluster].cluster_id = parent_cluster
        for child_cluster in joined_cluster1, joined_cluster2:
            nodes[parent_cluster].add_child(nodes[child_cluster])
            nodes[child_cluster].cluster_id = child_cluster
            nodes[child_cluster].edge_length = nantonone(cluster_distance / 2)

    return dendropy.Tree(seed_node=nodes[parent_cluster])
