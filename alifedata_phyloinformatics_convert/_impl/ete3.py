from unittest.mock import MagicMock
import warnings

try:
    import ete3
except (ImportError, ValueError, ModuleNotFoundError):  # pragma: no cover
    warnings.warn(
        "ImportWarning: ete3 import failed; "
        "inserting a no-op mock for ete3. "
        "This is likely because ete3 is not installed or incompatible with "
        "Python >= 3.13.",
    )
    ete3 = MagicMock()
