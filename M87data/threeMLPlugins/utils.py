import numpy as np
import os
import astropy.units as u
try:
    from threeML import Log_uniform_prior, Gaussian, Uniform_prior
except ImportError as e:
    raise ImportError(
        "ThreeML not installed"
    ) from e

def get_Value_Bounds_Prior(conf):
    """ calculate value, bounds and prior of parameter from config dict and
        make sure units are in units of keV

    Args:
        conf (dict): dict with value, unit, min, max and prior

    Returns:
        [float, [float, float], Prior]: value, bounds, prior
    """
    # check dimension
    unit = u.Unit(conf["unit"])
    target_unit = 1
    if u.get_physical_type(unit) == "energy":
        # energy
        target_unit = (1*unit).to(u.keV).value
    elif u.get_physical_type(unit*u.cm**2*u.s) == "energy":
        # energy flux
        target_unit = (1*unit).to(u.keV/u.cm**2/u.s).value
    value = conf["value"] * target_unit
    lower = None
    if "min" in conf:
        lower = conf["min"] * target_unit
    upper = None
    if "max" in conf:
        upper = conf["max"] * target_unit
    prior = None
    if "prior" in conf:
        if conf["prior"] == "Uniform":
            prior = Uniform_prior(
                lower_bound=lower, upper_bound=upper)
        elif conf["prior"] == "LogUniform":
            prior = Log_uniform_prior(
                lower_bound=lower, upper_bound=upper)
        elif conf["prior"] == "Gaussian":
            if not "mu" in conf:
                print("no mu in conf for parameter prior!")
            prior = Gaussian(
                mu=conf["mu"], sigma=conf["sigma"])

    return value, [lower, upper], prior


def parse_effAreaCorrection(det, configs, time):
    if "effAreaCor" in configs:
        configs = configs["effAreaCor"]
        if configs["do"]:
            det.use_effective_area_correction(*configs["bounds"])
            if configs["fix"]:
                if isinstance(configs["value"], list):
                    val = configs["value"][time]
                else:
                    val = configs["value"]
                det.fix_effective_area_correction(val)
    return det