import numpy as np
from importlib import resources

try:
    from threeML.io.logging import setup_logger
except ImportError as e:
    raise ImportError("ThreeML not installed") from e

from threeML_extras.EFELike import EFELike, keV2erg
from ..addM87SED import get_SED_data

h = 6.62e-27  # cgs
log = setup_logger(__name__)


class M87_SED_Plugin(EFELike):
    """
    Plugin for waveband-sliced SED datasets loaded via get_SED_data.

    Parameters
    ----------
    name : str
    waveband : str
        One of 'radio', 'optical', 'xray', 'gev', 'tev'.
    dataset : str
        Dataset identifier passed to get_SED_data.
    D_Mpc, MBH_MSUN, theta_view_deg : float
    """

    WAVEBANDS = {}

    def __init__(self, name, waveband, dataset,
                 D_Mpc=16.8, MBH_MSUN=6.5e9, theta_view_deg=17):
        waveband = waveband.strip()
        if waveband not in self.WAVEBANDS:
            raise ValueError(f"waveband '{waveband}' not known, choose from: {', '.join(self.WAVEBANDS)}")

        log.info(
            f"{self.__class__.__name__}: dataset='{dataset}', waveband='{waveband}', "
            f"D_Mpc={D_Mpc}, MBH_MSUN={MBH_MSUN:.2e}, theta_view_deg={theta_view_deg}"
        )

        _, df_data, _, _ = get_SED_data(dataset=dataset, D_Mpc=D_Mpc,
                                         MBH_MSUN=MBH_MSUN, theta_view_deg=theta_view_deg)
        wb = self.WAVEBANDS[waveband]
        self.df = df_data.loc[wb["min"]:wb["max"]].copy()

        x_keV = self.df["Frequency_Hz"].values * h / keV2erg
        y_efe = self.df["nuFnu_1e-12ergscm2"].values * 1e-12
        y_unc = self.df["sigma_nuFnu_1e-12ergscm2"].values * 1e-12

        super().__init__(name, x_keV=x_keV, y_efe=y_efe, y_unc=y_unc)

    def add2ax_SED_eV_ergscm2(self, ax, color="k", 
                              plot_effective_area_correction=False, **kwargs):
        """Plot using the original df columns directly. 
            plot_effective_area_correction is only rescaling data and uncertainty,
            but not propagating the uncertainties from the eff. area corr.
        """
        effarcorr = 1
        if plot_effective_area_correction:
            effarcorr = self._nuisance_parameter.value
        ax.errorbar(
            self.df["Frequency_Hz"] * h / 1.60218e-12,  # erg → eV
            effarcorr * self.df["nuFnu_1e-12ergscm2"] * 1e-12,
            yerr=effarcorr * self.df["sigma_nuFnu_1e-12ergscm2"] * 1e-12,
            ls="", marker=".", c=color, **kwargs,
        )
        return ax

