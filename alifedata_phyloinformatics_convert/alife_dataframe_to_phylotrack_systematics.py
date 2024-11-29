import tempfile

import pandas as pd

from ._impl import phytrack_Systematics


def alife_dataframe_to_phylotrack_systematics(
    df: pd.DataFrame,
) -> phytrack_Systematics:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as a phylotrackpy Systematics object.

    Notes
    -----
    Edge length support is not yet implemented.
    """
    with tempfile.NamedTemporaryFile() as tmp:
        df.to_csv(tmp.name, index=False)
        res = phytrack_Systematics(lambda x: x)
        res.load_from_file(tmp.name, "id", True)
        return res
