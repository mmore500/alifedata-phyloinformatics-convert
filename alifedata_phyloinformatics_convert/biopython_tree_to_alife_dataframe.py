from Bio import Phylo
from collections import defaultdict
from collections.abc import Mapping
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


def biopython_tree_to_alife_dataframe(
    tree: Phylo.BaseTree,
    exportattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
) -> pd.DataFrame:
    """Convert a biopython phylogenetic tree to a dataframe formatted to the
    artificial life communit data format standards.

    The following Clade object attributes will automatically be exported to
    dataframe columns, if available:
        * branch_length,
        * id,
        * name, and
        * origin_time.

    Parameters
    ----------
    tree:
        biopython tree to convert.
    exportattrs: optional
        Clade attrs that should be copied as columns into the generated
        dataframe. If a map is provided, attr values in keys will be inserted
        into the dataframe with the corresponding value as the column name.
    """

    # adapted from https://biopython.org/wiki/Phylo_cookbook
    def all_parents(tree):
        parents = defaultdict(lambda: None)
        for clade in tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents
    parents = all_parents(tree)

    # set up clade origin times if any edge lengths set
    if any(clade.branch_length is not None for clade in tree.find_clades()):
        if not hasattr(tree.root, 'origin_time'):
            tree.root.origin_time = opyt.or_value(
                tree.root.branch_length,
                0
            )
        else:
            assert tree.root.origin_time is not None
        for clade in tree.find_clades(order='level'):
            parent = parents[clade]
            if parent is not None and not hasattr(clade, 'origin_time'):
                if None not in (
                    getattr(parent, 'origin_time', None),
                    clade.branch_length,
                ) and not isanan(parent.origin_time):
                    clade.origin_time \
                        = parent.origin_time + clade.branch_length

    # attach ids to clades, if needed
    for fallback_id, clade in enumerate(tree.find_clades(order='preorder')):
        if not hasattr(clade, 'id'):
            clade.id = fallback_id
        else:
            assert isinstance(clade.id, int)

    assert not any((
        attr in opyt.or_value(exportattrs, [])
        for attr in ('origin_time', 'id', 'branch_length', 'name')
    ))
    return pd.DataFrame.from_records([
        {
            **{
                'id': clade.id,
                'ancestor_list': str(
                    [opyt.apply_if(parents[clade], lambda x: x.id)]
                ),
                'origin_time': getattr(clade, 'origin_time', None),
                'branch_length': clade.branch_length,
                'name': clade.name,
            },
            **{
                exportattrs[attr_name]
                if isinstance(exportattrs, Mapping)
                else attr_name: _rgetattr(clade, attr_name)
                for attr_name in opyt.or_value(exportattrs, [])
            },
        }
        for clade in tree.find_clades(order='preorder')
    ])
