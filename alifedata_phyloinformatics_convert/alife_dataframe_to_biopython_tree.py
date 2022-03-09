from Bio import Phylo
from iterpop import iterpop as ip
import pandas as pd
import typing

from .alife_dataframe_to_biopython_trees \
    import alife_dataframe_to_biopython_trees


def alife_dataframe_to_biopython_tree(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_branch_lengths: bool = False,
) -> typing.Optional[Phylo.BaseTree.Tree]:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as a biopython tree.

    Returns None if df is empty.

    If two or more clades exist that do not share a common ancestor, an
    exception will be raised.

    See Also
    ----------
    alife_dataframe_to_biopython_trees
    """
    return ip.poursingleton(
        alife_dataframe_to_biopython_trees(
            df,
            setattrs=setattrs,
            setup_branch_lengths=setup_branch_lengths,
        ),
    )
