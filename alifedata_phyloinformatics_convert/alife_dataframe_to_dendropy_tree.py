import dendropy
from iterpop import iterpop as ip
import pandas as pd
import typing

from .alife_dataframe_to_dendropy_trees \
    import alife_dataframe_to_dendropy_trees


def alife_dataframe_to_dendropy_tree(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_edge_lengths: bool = False,
) -> typing.Optional[dendropy.Tree]:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as a dendropy tree.

    Returns None if df is empty.

    If two or more clades exist that do not share a common ancestor, an
    exception will be raised.

    See Also
    ----------
    alife_dataframe_to_dendropy_trees
    """

    return ip.poursingleton(
        alife_dataframe_to_dendropy_trees(
            df,
            setattrs=setattrs,
            setup_edge_lengths=setup_edge_lengths,
        ),
    )
