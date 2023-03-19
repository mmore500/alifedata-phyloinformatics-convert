# This file has been copied from hstrat
# <https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_alifestd_is_asexual.py>.
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.
import pandas as pd

from .alifestd_is_sexual import alifestd_is_sexual


def alifestd_is_asexual(phylogeny_df: pd.DataFrame) -> bool:
    """Do all organisms in the phylogeny have one or no immediate ancestor?
    Input dataframe is not mutated by this operation.
    """
    return "ancestor_id" in phylogeny_df or not alifestd_is_sexual(
        phylogeny_df
    )
