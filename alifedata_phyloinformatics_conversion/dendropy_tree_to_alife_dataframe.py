import dendropy
from nanto import isanan
import opytional as opyt
import pandas as pd


def dendropy_tree_to_alife_dataframe(tree: dendropy.Tree) -> pd.DataFrame:

    # set up node origin times if any edge lengths set
    if any(node.edge_length is not None for node in tree):
        if not hasattr(tree.seed_node, 'origin_time'):
            tree.seed_node.origin_time = opyt.or_value(
                tree.seed_node.edge_length,
                0
            )
        else:
            assert tree.seed_node.origin_time is not None
        for node in tree:
            parent = node.parent_node
            if parent is not None and not hasattr(node, 'origin_time'):
                if None not in (
                    getattr(parent, 'origin_time', None),
                    node.edge_length,
                ) and not isanan(parent.origin_time):
                    node.origin_time = parent.origin_time + node.edge_length

    # attach ids to nodes, if needed
    for fallback_id, node in enumerate(tree):
        if not hasattr(node, 'id'):
            node.id = fallback_id
        else:
            assert isinstance(node.id, int)

    return pd.DataFrame.from_records([
        {
            'id': node.id,
            'ancestor_list': str(
                [opyt.apply_if(node.parent_node, lambda x: x.id)]
            ),
            'origin_time': getattr(node, 'origin_time', None),
            'edge_length': node.edge_length,
            'label': node.label,
            'taxon_label': opyt.apply_if(node.taxon, lambda x: x.label),
        }
        for node in tree
    ])
