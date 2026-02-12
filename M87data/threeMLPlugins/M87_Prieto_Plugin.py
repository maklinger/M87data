import numpy as np
try:
    from threeML import PluginPrototype
    from threeML.io.logging import setup_logger
    from threeML.plugins.XYLike import _chi2_like
except ImportError as e:
    raise ImportError(
        "ThreeML not installed"
    ) from e
try:
    from astromodels.core.parameter import Parameter
    from astromodels.functions.priors import Uniform_prior
except ImportError as e:
    raise ImportError(
        "astromodels not installed"
    ) from e
import collections
from typing import Any, Dict, List, Optional, Tuple, Union
import astropy.units as u

from ..addM87SED import get_SED_data

eV2erg = 1.60218e-12
erg2eV = 1/eV2erg
keV2erg = 1e3 * eV2erg
erg2keV = 1/keV2erg
h = 6.62e-27 # cgs

wavebands = {
    "radio": {"min": 0, "max": 11},
    "ir": {"min": 12, "max": 13},
    "optical": {"min": 14, "max": 26},
    "xray": {"min": 27, "max": 28}
}

log = setup_logger(__name__)


class M87_Prieto_Plugin(PluginPrototype):
    def __init__(self, name, waveband,
        D_Mpc=16.8, MBH_MSUN=6.5e9, theta_view_deg=17):

        waveband = waveband.strip()
        if waveband not in wavebands:
            valid = ", ".join(wavebands)
            raise ValueError(
                f"waveband '{waveband}' not known, choose from: {valid}"
            )

        self.waveband = waveband

        df, df_data, df_VM, df_UL = get_SED_data(dataset="M87SED_Prieto_quiet",
            D_Mpc=D_Mpc, MBH_MSUN=MBH_MSUN, theta_view_deg=theta_view_deg)
        self.df = df_data.loc[wavebands[waveband]["min"]:wavebands[waveband]["max"]].copy()

        self.x = self.df["Frequency_Hz"].values * h * u.erg
        self.x_keV = self.x.to(u.keV).value
        self.x_erg = self.x.to(u.erg).value

        self.y = self.df["nuFnu_1e-12ergscm2"].values * 1e-12 * u.erg/u.cm**2/u.s
        self.yerr = self.df["sigma_nuFnu_1e-12ergscm2"].values * 1e-12 * u.erg/u.cm**2/u.s
        # self.lgy = np.log(self.y.to(u.erg/u.cm**2/u.s).value)
        # self.lgyerr_up = np.log(self.yu / self.y)
        # self.lgyerr_low = - np.log(self.yl / self.y)
        # self.lgyerr = np.maximum(self.lgyerr_low, self.lgyerr_up)

        # create the hash for the nuisance parameters
        # nuisance_parameters = collections.OrderedDict()
        self._nuisance_parameter: Parameter = Parameter(
            "cons_%s" % name,
            1.0,
            min_value=0.8,
            max_value=1.2,
            delta=0.05,
            free=False,
            desc="Effective area correction for %s" % name,
        )

        nuisance_parameters: Dict[str, Parameter] = collections.OrderedDict()
        nuisance_parameters[
            self._nuisance_parameter.name] = self._nuisance_parameter

        # call the prototype constructor
        super().__init__(name, nuisance_parameters)

    def get_number_of_data_points(self):
        """
        This returns the number of data points that are used to evaluate the likelihood.
        For binned measurements, this is the number of active bins used in the fit. For
        unbinned measurements, this would be the number of photons/particles that are
        evaluated on the likelihood
        """
        return len(self.x)

    def set_model(self, model):

        # attach the model to the object
        self._model = model
    
    def get_log_like(self):
        n_point_sources = self._model.get_number_of_point_sources()

        assert (
            n_point_sources > 0
        ), "You need to have at least one point source defined"
        assert (
            self._model.get_number_of_extended_sources() == 0
        ), "MAGICLike does not support extended sources"

        # astropy implementation is slower
        # expectation = np.sum([
        #     (self.x**2 * source(self.x, tag=self._tag)).to(
        #         u.erg / u.cm ** 2 / u.s).value
        #     for source in list(self._model.point_sources.values())
        # ], axis=0)
        
        # Make a function which will stack all point sources
        expectation = np.sum([
            self.x_erg * self.x_keV * source(self.x_keV, tag=self._tag)
            for source in list(self._model.point_sources.values())
        ], axis=0)
        if self._nuisance_parameter.value == 1:
            return _chi2_like(self.y, self.yerr, expectation)
        else:
            return _chi2_like(
                self._nuisance_parameter.value * self.y, 
                self._nuisance_parameter.value * self.yerr, 
                expectation
            )

    def inner_fit(self):

        return self.get_log_like()
    
    def use_effective_area_correction(self,
                                      min_value: Union[int, float] = 0.8,
                                      max_value: Union[int, float] = 1.2) -> None:
        """
        Activate the use of the effective area correction, which is a multiplicative factor in front of the model which
        might be used to mitigate the effect of intercalibration mismatch between different instruments.

        NOTE: do not use this is you are using only one detector, as the multiplicative constant will be completely
        degenerate with the normalization of the model.

        NOTE2: always keep at least one multiplicative constant fixed to one (its default value), when using this
        with other OGIPLike-type detectors

        :param min_value: minimum allowed value (default: 0.8, corresponding to a 20% - effect)
        :param max_value: maximum allowed value (default: 1.2, corresponding to a 20% + effect
        :return:
        """
        log.info(
            f"{self._name} is using effective area correction (between {min_value} and {max_value})")
        self._nuisance_parameter.free = True
        self._nuisance_parameter.bounds = (min_value, max_value)

        # Use a uniform prior by default

        self._nuisance_parameter.set_uninformative_prior(Uniform_prior)


    def fix_effective_area_correction(self,
                                      value: Union[int, float] = 1) -> None:
        """
        Fix the multiplicative factor (see use_effective_area_correction) to the provided value (default: 1)

        :param value: new value (default: 1, i.e., no correction)
        :return:
        """
        log.info(
            f"{self._name} is using a fixed effective area correction with {value}")
        self._nuisance_parameter.value = value
        self._nuisance_parameter.fix = True

    def add2ax_SED_eV_ergscm2(self, ax, color="k"):
        ax.errorbar(self.x_keV*1e3, self.y, yerr=self.yerr, ls="", marker=".", c=color)
        return ax
        

