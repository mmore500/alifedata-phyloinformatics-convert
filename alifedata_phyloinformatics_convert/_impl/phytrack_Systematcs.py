from unittest.mock import MagicMock
import warnings

try:
    from phylotrackpy.systematics import Systematics as phytrack_Systematics
except (ImportError, ValueError):  # pragma: no cover
    warnings.warn(
        "ImportWarning: phylotrackpy.Systematics import failed; "
        "inserting a no-op mock for Systematics. "
        "This is likely because phylotrackpy is not installed.",
    )
    phytrack_Systematics = MagicMock()
