import opytional as opyt
import pandas as pd
import typing

from ._impl import parse_ancestor_list as _parse_ancestor_list


def alife_dataframe_to_dict_of_lists(
    df: pd.DataFrame,
) -> typing.Dict[int, typing.List[int]]:
    """Extract an adjacency representation for an alife standard dataframe.

    Returns
    -------
    dict of int : list of int
        A dictionary where the keys are organism IDs and the values are the
        lists of organisms' ancestors.
    """

    df = df.copy()

    df['parsed_ancestor_list'] \
        = df['ancestor_list'].apply(_parse_ancestor_list)

    return df.set_index('id', drop=False)['parsed_ancestor_list'].to_dict()
