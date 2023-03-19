# This file has been copied from hstrat
# <https://github.com/mmore500/hstrat/blob/47fe9048d6b327ca40f77424e3cf9392f7980689/hstrat/_auxiliary_lib/_all_unique.py>.
# It is up to date as of commit
# [1ce7a0e](https://github.com/mmore500/hstrat/commit/1ce7a0e).
# Eventually, it should be migrated to a standalone library.
import numpy as np


def all_unique(array: np.array) -> bool:
    """Are all values in `array` unique?"""
    seen_unique = set()
    for idx, val in enumerate(array):
        seen_unique.add(val)
        if idx + 1 != len(seen_unique):
            return False

    return True
