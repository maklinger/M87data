from .M87_SED_Plugin import M87_SED_Plugin


class M87_2017_Plugin(M87_SED_Plugin):
    """EHT MWL 2017 dataset."""

    WAVEBANDS = {
        "radio":   {"min": 0,  "max": 16},
        "optical": {"min": 17, "max": 25},
        "xray":    {"min": 26, "max": 45},
        "gev":     {"min": 46, "max": 49},
        "tev":     {"min": 50, "max": 66},
    }

    def __init__(self, name, waveband, D_Mpc=16.8, MBH_MSUN=6.5e9, 
                 theta_view_deg=17, systematic_fraction=0.0):
        super().__init__(name, waveband, "M87SED_EHTMWL2017", D_Mpc, MBH_MSUN, 
                         theta_view_deg, systematic_fraction)


