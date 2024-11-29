import tempfile

import pandas as pd

from ._impl import phytrack_Systematics


def phylotrack_systematics_to_alife_dataframe(
    systematics: phytrack_Systematics,
) -> pd.DataFrame:
    """Convert a phylotrackpy Systematics object to a dataframe formatted to
    the artificial life community data format standards.

    Parameters
    ----------
    systematics:
        The phylotrackpy Systematics object to convert.
    """
    with tempfile.NamedTemporaryFile() as tmp:
        systematics.snapshot(tmp.name)
        res = pd.read_csv(tmp.name)
        res["ancestor_list"].replace('["NONE"]', "[None]", inplace=True)
        return res
