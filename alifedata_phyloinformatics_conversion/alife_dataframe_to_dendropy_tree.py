import dendropy
from iterpop import iterpop as ip
import pandas as pd
import typing

from .alife_dataframe_to_dendropy_trees \
    import alife_dataframe_to_dendropy_trees


def alife_dataframe_to_dendropy_tree(
    df: pd.DataFrame,
    setup_edge_lengths: bool = False,
) -> typing.Optional[dendropy.Tree]:
    return ip.poursingleton(
        alife_dataframe_to_dendropy_trees(
            df,
            setup_edge_lengths=setup_edge_lengths,
        ),
    )
