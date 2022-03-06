from collections import defaultdict
import dendropy
from nanto import nantonone
import numpy as np


def scipy_linkage_matrix_to_dendropy_tree(matrix: np.array) -> dendropy.Tree:
    # scipy linkage format
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    # cluster id -> node
    nodes = defaultdict(dendropy.Node)
    for row_idx, row in enumerate(matrix):
        num_rows = len(matrix)
        parent_cluster = row_idx + num_rows + 1
        joined_cluster1, joined_cluster2, cluster_distance, cluster_size = row

        nodes[parent_cluster].cluster_id = parent_cluster
        for child_cluster in joined_cluster1, joined_cluster2:
            nodes[parent_cluster].add_child(nodes[child_cluster])
            nodes[child_cluster].cluster_id = child_cluster
            nodes[child_cluster].edge_length = nantonone(cluster_distance / 2)

    return dendropy.Tree(seed_node=nodes[parent_cluster])
