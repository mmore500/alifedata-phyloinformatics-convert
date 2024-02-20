# This file has been copied from hstrat
# <https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_alifestd_make_ancestor_list_col.py>.
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.

import pandas as pd


def alifestd_make_ancestor_list_col(
    ids: pd.Series, ancestor_ids: pd.Series
) -> pd.Series:
    """Translate a column of integer ancestor id values into alife standard
    `ancestor_list` representation."""

    res = ancestor_ids.map("[{}]".format)
    res[ids == ancestor_ids] = "[none]"

    return res
