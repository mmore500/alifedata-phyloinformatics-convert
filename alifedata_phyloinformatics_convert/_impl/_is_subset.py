# This file has been copied from hstrat
# <https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_is_subset.py>.
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.
import numpy as np


def is_subset(subset: np.array, superset: np.array) -> bool:
    """Are all values in `subset` contained in `superset`?"""
    superset_lookup = set(superset)
    for val in subset:
        if val not in superset_lookup:
            return False

    return True
