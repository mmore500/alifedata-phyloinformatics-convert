# This file has been copied from hstrat (https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_alifestd_validate.py).
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.

import pandas as pd

from .alifestd_is_asexual import alifestd_is_asexual
from .alifestd_make_ancestor_id_col import alifestd_make_ancestor_id_col
from .alifestd_make_ancestor_list_col import alifestd_make_ancestor_list_col
from .parse_ancestor_list import parse_ancestor_list
from .all_unique import all_unique
from .is_subset import is_subset


def _validate_ancestors_asexual(
    phylogeny_df: pd.DataFrame, mutate: bool
) -> bool:
    if "ancestor_id" not in phylogeny_df:
        if not mutate:
            phylogeny_df = phylogeny_df.copy()
        phylogeny_df["ancestor_id"] = alifestd_make_ancestor_id_col(
            phylogeny_df["id"], phylogeny_df["ancestor_list"]
        )
    elif not (
        phylogeny_df["ancestor_list"].str.lower().replace("[]", "[none]")
        == alifestd_make_ancestor_list_col(
            phylogeny_df["id"], phylogeny_df["ancestor_id"]
        )
    ).all():
        return False

    return is_subset(
        phylogeny_df["ancestor_id"].to_numpy(),
        phylogeny_df["id"].to_numpy(),
    )


def _validate_ancestors_sexual(phylogeny_df: pd.DataFrame) -> bool:
    ids = set(phylogeny_df["id"])
    return all(
        ancestor_id in ids
        for ancestor_list_str in phylogeny_df["ancestor_list"]
        for ancestor_id in parse_ancestor_list(ancestor_list_str)
        if ancestor_id is not None
    )


def __try_alifestd_validate(
    phylogeny_df: pd.DataFrame,
    mutate: bool,
) -> bool:
    has_mandatory_columns = (
        "id" in phylogeny_df and "ancestor_list" in phylogeny_df
    )
    if not has_mandatory_columns:
        return False

    ids_valid = (
        pd.api.types.is_integer_dtype(phylogeny_df["id"])
        and (phylogeny_df["id"] >= 0).all()
        and all_unique(phylogeny_df["id"].to_numpy())
    )
    if not ids_valid:
        return False

    ancestor_lists_syntax_valid = (
        pd.api.types.is_string_dtype(phylogeny_df["ancestor_list"])
        and phylogeny_df["ancestor_list"].str.startswith("[").all()
        and phylogeny_df["ancestor_list"].str.endswith("]").all()
    )
    if not ancestor_lists_syntax_valid:
        return False

    return (
        _validate_ancestors_asexual(phylogeny_df, mutate)
        if alifestd_is_asexual(phylogeny_df)
        else _validate_ancestors_sexual(phylogeny_df)
    )


def _try_alifestd_validate(
    phylogeny_df: pd.DataFrame,
    mutate: bool = False,
) -> bool:
    """Is the phylogeny compliant to alife data standards?
    Input dataframe is not mutated by this operation unless `mutate` set True.
    """

    try:
        return __try_alifestd_validate(phylogeny_df, mutate)
    except:
        return False
