
import numpy as np
from importlib import resources

try:
    from threeML.io.logging import setup_logger
except ImportError as e:
    raise ImportError("ThreeML not installed") from e

from threeML_extras.EFELike import EFELike
from ..addM87SED import get_SED_data

h = 6.62e-27  # cgs
log = setup_logger(__name__)


class M87_JWST_Plugin(EFELike):
    """
    Plugin for JWST averaged SED data.

    Parameters
    ----------
    name : str
    N_av : int
        Averaging window of the JWST SED file (default 300).
    """
    def __init__(self, name, N_av=300):

        available = get_available_N_av()
        if N_av not in available:
            raise ValueError(
                f"N_av={N_av} not available, choose from: {available}"
            )
        path = resources.files("M87data.data").joinpath(f"SED/JWST_average_{N_av}.csv")
        E_eV, EFE, unc_EFE = np.genfromtxt(path, delimiter=",", skip_header=1).T

        log.info(f"M87_JWST_Plugin: loaded JWST_average_{N_av}.csv")
        super().__init__(name, x_keV=E_eV * 1e-3, y_efe=EFE, y_unc=unc_EFE)


def get_available_N_av() -> list[int]:
    """Return available averaging windows based on existing JWST SED files."""
    data_dir = resources.files("M87data.data").joinpath("SED")
    return sorted([
        int(p.name.split("_")[-1].split(".")[0])
        for p in data_dir.iterdir()
        if p.name.startswith("JWST_average_") and p.name.endswith(".csv")
    ])