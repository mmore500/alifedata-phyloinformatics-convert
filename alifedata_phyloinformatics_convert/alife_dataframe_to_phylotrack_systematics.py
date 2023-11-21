import tempfile

import pandas as pd
from phylotrackpy.systematics import Systematics as phytrack_Systematics


def alife_dataframe_to_phylotrack_systematics(
    df: pd.DataFrame,
) -> phytrack_Systematics:
    with tempfile.NamedTemporaryFile() as tmp:
        df.to_csv(tmp.name, index=False)
        res = phytrack_Systematics(lambda x: x)
        res.load_from_file(tmp.name, "id", True)
        return res
