from Bio import Phylo
from collections import defaultdict
from nanto import isanan
import opytional as opyt
import pandas as pd


def biopython_tree_to_alife_dataframe(tree: Phylo.BaseTree) -> pd.DataFrame:

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

    return pd.DataFrame.from_records([
        {
            'id': clade.id,
            'ancestor_list': str(
                [opyt.apply_if(parents[clade], lambda x: x.id)]
            ),
            'origin_time': getattr(clade, 'origin_time', None),
            'branch_length': clade.branch_length,
            'name': clade.name,
        }
        for clade in tree.find_clades(order='preorder')
    ])
