import numpy as np

def is_subset(subset: np.array, superset: np.array) -> bool:
    """Are all values in `subset` contained in `superset`?"""
    superset_lookup = set(superset)
    for val in subset:
        if val not in superset_lookup:
            return False

    return True
