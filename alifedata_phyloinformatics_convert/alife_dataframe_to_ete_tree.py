from iterpop import iterpop as ip
import pandas as pd
import typing

from .alife_dataframe_to_ete_trees import alife_dataframe_to_ete_trees
from ._impl import ete3


def alife_dataframe_to_ete_tree(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_dists: bool = False,
) -> typing.Optional[ete3.Tree]:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as a ete tree.

    Returns None if df is empty.

    If two or more clades exist that do not share a common ancestor, an
    exception will be raised.

    See Also
    ----------
    alife_dataframe_to_ete_trees
    """

    return ip.poursingleton(
        alife_dataframe_to_ete_trees(
            df,
            setattrs=setattrs,
            setup_dists=setup_dists,
        ),
    )
