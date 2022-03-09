from collections.abc import Mapping
import dendropy
import functools
from nanto import isanan
import opytional as opyt
import pandas as pd
import typing


# adapted from https://stackoverflow.com/a/31174427/17332200
def _rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)
    return functools.reduce(_getattr, [obj] + attr.split('.'))


def dendropy_tree_to_alife_dataframe(
    tree: dendropy.Tree,
    exportattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
) -> pd.DataFrame:
    """Convert a dendropy phylogenetic tree to a dataframe formatted to the
    artificial life communit data format standards.

    The following Node object attributes will automatically be exported to
    dataframe columns, if available:
        * edge_length,
        * id,
        * label,
        * origin_time, and
        * taxon.label.

    Parameters
    ----------
    tree:
        dendropy tree to convert.
    exportattrs: optional
        Node attrs that should be copied as columns into the generated
        dataframe. If a map is provided, attr values in keys will be inserted
        into the dataframe with the corresponding value as the column name.
    """

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

    assert not any((
        attr in opyt.or_value(exportattrs, [])
        for attr
        in ('origin_time', 'id', 'edge_length', 'label', 'taxon.label')
    ))
    return pd.DataFrame.from_records([
        {
            **{
                'id': node.id,
                'ancestor_list': str(
                    [opyt.apply_if(node.parent_node, lambda x: x.id)]
                ),
                'origin_time': getattr(node, 'origin_time', None),
                'edge_length': node.edge_length,
                'label': node.label,
                'taxon_label': opyt.apply_if(node.taxon, lambda x: x.label),
            },
            **{
                exportattrs[attr_name]
                if isinstance(exportattrs, Mapping)
                else attr_name: _rgetattr(node, attr_name)
                for attr_name in opyt.or_value(exportattrs, [])
            },
        }
        for node in tree
    ])
