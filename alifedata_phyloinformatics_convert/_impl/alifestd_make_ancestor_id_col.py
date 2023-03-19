# This file has been copied from hstrat
# <https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_alifestd_make_ancestor_id_col.py>.
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.

import pandas as pd


def alifestd_make_ancestor_id_col(
    ids: pd.Series, ancestor_lists: pd.Series
) -> pd.Series:
    """Translate ancestor ids from a column of singleton `ancestor_list`s into a
    pure-integer series representation.

    Each organism must have one or zero ancestors (i.e., asexualasexual data).
    In the returned series, ancestor id will be assigned to own id for no-
    ancestor organisms.
    """
    ancestor_ids = (
        ancestor_lists.str.lower()
        .replace("[none]", "[-1]")
        .replace("[]", "[-1]")
        .str.strip("[]")
        .astype(int)
    )

    root_filter = ancestor_ids == -1
    ancestor_ids[root_filter] = ids[root_filter]

    return ancestor_ids
