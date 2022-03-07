from Bio import Phylo
from iterpop import iterpop as ip
import pandas as pd
import typing

from .alife_dataframe_to_biopython_trees \
    import alife_dataframe_to_biopython_trees


def alife_dataframe_to_biopython_tree(
    df: pd.DataFrame,
    setup_branch_lengths: bool = False,
) -> typing.Optional[Phylo.BaseTree.Tree]:
    return ip.poursingleton(
        alife_dataframe_to_biopython_trees(
            df,
            setup_branch_lengths=setup_branch_lengths,
        ),
    )
