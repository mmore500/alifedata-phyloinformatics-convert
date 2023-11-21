import tempfile

import pandas as pd
from phylotrackpy.systematics import Systematics as phytrack_Systematics


def phylotrack_systematics_to_alife_dataframe(
    tree: phytrack_Systematics,
) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile() as tmp:
        tree.snapshot(tmp.name)
        res = pd.read_csv(tmp.name)
        res["ancestor_list"].replace('["NONE"]', "[None]", inplace=True)
        return res
