from collections.abc import Mapping
from nanto import isanan
import opytional as opyt
import pandas as pd
import typing

from ._impl import ete3
from ._impl import rgetattr as _rgetattr


def ete_tree_to_alife_dataframe(
    tree: typing.Union[ete3.Tree, ete3.TreeNode],
    exportattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
) -> pd.DataFrame:
    """Convert a ete phylogenetic tree to a dataframe formatted to the
    artificial life communit data format standards.

    The following TreeNode object attributes will automatically be exported to
    dataframe columns, if available:
        * dist,
        * id,
        * origin_time, and
        * name.

    Parameters
    ----------
    tree:
        ete tree to convert.
    exportattrs: optional
        Node attrs that should be copied as columns into the generated
        dataframe. If a map is provided, attr values in keys will be inserted
        into the dataframe with the corresponding value as the column name.
    """

    # set up node origin times if any edge lengths set
    if any(node.dist != 1.0 for node in tree):
        if not hasattr(tree, 'origin_time'):
            tree.add_features(origin_time=tree.dist)
        else:
            assert tree.origin_time is not None
        for node in tree:
            parent = next(node.iter_ancestors(), None)
            if parent is not None and not hasattr(node, 'origin_time'):
                if (
                    getattr(parent, 'origin_time', None) is not None
                    and not isanan(parent.origin_time)
                ):
                    node.add_features(
                        origin_time=parent.origin_time + node.dist,
                    )

    # attach ids to nodes, if needed
    for fallback_id, node in enumerate(tree.traverse()):
        if not hasattr(node, 'id'):
            node.add_features(id=fallback_id)
        else:
            assert isinstance(node.id, int)

    assert not any((
        attr in (
            exportattrs.values()
            if isinstance(exportattrs, Mapping)
            else opyt.or_value(exportattrs, [])
        )
        for attr in ('origin_time', 'id', 'dist', 'name')
    ))
    return pd.DataFrame.from_records([
        {
            **{
                'id': node.id,
                'ancestor_list':
                    str([opyt.apply_if(node.up, lambda x: x.id)]),
                'origin_time': getattr(node, 'origin_time', None),
                'dist': node.dist,
                'name': node.name,
            },
            **{
                exportattrs[attr_name]
                if isinstance(exportattrs, Mapping)
                else attr_name: _rgetattr(node, attr_name)
                for attr_name in opyt.or_value(exportattrs, [])
            },
        }
        for node in tree.traverse()
    ])
