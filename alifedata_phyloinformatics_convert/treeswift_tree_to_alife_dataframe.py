from collections.abc import Mapping
from nanto import isanan
import opytional as opyt
import pandas as pd
import treeswift
import typing

from ._impl import rgetattr as _rgetattr


def treeswift_tree_to_alife_dataframe(
    tree: treeswift.Tree,
    exportattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
) -> pd.DataFrame:
    """Convert a treeswift phylogenetic tree to a dataframe formatted to the
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
        treeswift tree to convert.
    exportattrs: optional
        Node attrs that should be copied as columns into the generated
        dataframe. If a map is provided, attr values in keys will be inserted
        into the dataframe with the corresponding value as the column name.
    """

    # set up node origin times if any edge lengths set
    if any(node.edge_length is not None for node in tree.traverse_postorder()):
        if not hasattr(tree.root, 'origin_time'):
            tree.root.origin_time = opyt.or_value(
                tree.root.edge_length,
                0
            )
        else:
            assert tree.root.origin_time is not None
        for node in tree.traverse_postorder():
            parent = node.parent
            if parent is not None and not hasattr(node, 'origin_time'):
                if None not in (
                    getattr(parent, 'origin_time', None),
                    node.edge_length,
                ) and not isanan(parent.origin_time):
                    node.origin_time = parent.origin_time + node.edge_length

    # attach ids to nodes, if needed
    for fallback_id, node in enumerate(tree.traverse_postorder()):
        if not hasattr(node, 'id'):
            node.id = fallback_id
        else:
            assert isinstance(node.id, int)

    assert not any((
        attr in (
            exportattrs.values()
            if isinstance(exportattrs, Mapping)
            else opyt.or_value(exportattrs, [])
        )
        for attr in ('origin_time', 'id', 'edge_length', 'label', 'taxon.label')
    ))
    return pd.DataFrame.from_records([
        {
            **{
                'id': node.id,
                'ancestor_list': str(
                    [opyt.apply_if(node.parent, lambda x: x.id)]
                ),
                'origin_time': getattr(node, 'origin_time', None),
                'edge_length': node.edge_length,
                'label': node.label,
                'taxon_label': node.label,
            },
            **{
                exportattrs[attr_name]
                if isinstance(exportattrs, Mapping)
                else attr_name: _rgetattr(node, attr_name)
                for attr_name in opyt.or_value(exportattrs, [])
            },
        }
        for node in tree.traverse_postorder()
    ])
