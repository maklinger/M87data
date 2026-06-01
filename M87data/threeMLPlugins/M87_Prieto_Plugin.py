from .M87_SED_Plugin import M87_SED_Plugin


class M87_Prieto_Plugin(M87_SED_Plugin):
    """Prieto 2016 dataset. See M87_SED_Plugin for parameters."""

    WAVEBANDS = {
        "radio":   {"min":  0, "max": 11},
        "ir":      {"min": 12, "max": 13},
        "optical": {"min": 14, "max": 26},
        "xray":    {"min": 27, "max": 28},
    }

    def __init__(self, name, waveband, D_Mpc=16.8, MBH_MSUN=6.5e9, 
                 theta_view_deg=17, systematic_fraction=0.0):
        super().__init__(name, waveband, "M87SED_Prieto_quiet", D_Mpc, MBH_MSUN, 
                         theta_view_deg, systematic_fraction)