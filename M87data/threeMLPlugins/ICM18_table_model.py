import collections 
import numpy as np
import astropy.units as u
from astromodels.functions.function import (
    Function1D,
    FunctionMeta
)

from importlib import resources


class ICM18TableModel(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        A simple wrapper for the X-ray ICM model from the 2018 EHT MWL paper
    latex : $ICM$
    parameters :
        norm :
            desc : norm factor
            initial value : 1
            min : 0
            max : 1e5

    """


    def _setup(self):

        path = resources.files("M87data.data").joinpath(f"ICMTable/ICM_table.csv")
        self.Elo, self.Eup, self.flux = np.genfromtxt(
            path, skip_header=1, delimiter=",").T
        self.Emid = (self.Elo + self.Eup) / 2

    # -------------------------------------------------
    # Units
    # -------------------------------------------------
    def _set_units(self, x_unit, y_unit):
        self.norm.unit = u.dimensionless_unscaled

    # noinspection PyPep8Naming
    def evaluate(self, 
        x, # energy in keV
        norm
    ):

        # x_ in keV
        if isinstance(x, u.Quantity):
            # keep the values
            x_ = x.to(u.keV).value
            # keep the unit
            unit_ = self.y_unit
        else:

            # we do not need to do anything here
            x_ = x

            # this will basically be ignored
            unit_ = 1.0


        return unit_ * norm * np.exp(
            np.interp(
                np.log(x_),
                np.log(self.Emid),
                np.log(self.flux +1e-100),
                left=-1e100, right=-1e100))

