from copy import deepcopy
import dendropy
from iterpop import iterpop as ip
import itertools as it
import numpy as np
from sortedcontainers import SortedSet


def dendropy_tree_to_scipy_linkage_matrix(tree: dendropy.Tree) -> np.array:
    # scipy linkage format
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    # simplify tree
    tree = deepcopy(tree)
    tree.resolve_polytomies()
    tree.suppress_unifurcations()
    assert all(len(node.child_nodes()) in [0, 2] for node in tree)

    for node in tree.postorder_node_iter():
        node.num_leaf_descendants = max(
            sum(chld.num_leaf_descendants for chld in node.child_node_iter()),
            1,
        )

    cluster_id_generator = it.count()
    for leaf in tree.leaf_node_iter():
        if not hasattr(leaf, 'cluster_id'):
            leaf.cluster_id = next(cluster_id_generator)

    one_leaf_parents = {
        leaf.parent_node
        for leaf in tree.leaf_node_iter()
        if leaf.parent_node is not None
        and not ip.popsingleton(leaf.sibling_nodes()).is_leaf()
    }
    two_leaf_parents = SortedSet(
        (
            leaf.parent_node
            for leaf in tree.leaf_node_iter()
            if leaf.parent_node is not None
            and ip.popsingleton(leaf.sibling_nodes()).is_leaf()
        ),
        key=lambda node: sum(n.edge_length for n in node.child_node_iter()),
    )

    res = []
    while len(two_leaf_parents):
        two_leaf_parent = two_leaf_parents.pop(0)
        if not hasattr(two_leaf_parent, 'cluster_id'):
            two_leaf_parent.cluster_id = next(cluster_id_generator)
        if two_leaf_parent.parent_node is not None:
            if two_leaf_parent.parent_node in one_leaf_parents:
                one_leaf_parents.remove(two_leaf_parent.parent_node)
                two_leaf_parents.add(two_leaf_parent.parent_node)
            else:
                one_leaf_parents.add(two_leaf_parent.parent_node)

        child1, child2 = two_leaf_parent.child_node_iter()
        assert child1 not in two_leaf_parents
        assert child2 not in two_leaf_parents

        # see https://stackoverflow.com/a/40983611/17332200
        # for explainer on scipy linkage format
        joined_cluster1 = child1.cluster_id
        joined_cluster2 = child2.cluster_id
        assert None not in (joined_cluster1, joined_cluster2)
        cluster_distance = child1.edge_length + child2.edge_length
        cluster_size = two_leaf_parent.num_leaf_descendants

        res.append([
            joined_cluster1,
            joined_cluster2,
            float(cluster_distance),
            cluster_size,
        ])

    return np.array(res)
