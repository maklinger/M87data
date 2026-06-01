from .M87_SED_Plugin import M87_SED_Plugin


class M87_2018_Plugin(M87_SED_Plugin):
    """EHT MWL 2018 dataset."""

    WAVEBANDS = {
        "radio":   {"min": 0,  "max": 13},
        "optical": {"min": 13, "max": 26},
        "xray":    {"min": 26, "max": 47},
        "gev":     {"min": 47, "max": 51},
        "tev":     {"min": 51, "max": 61},
    }

    def __init__(self, name, waveband, D_Mpc=16.8, MBH_MSUN=6.5e9, theta_view_deg=17):
        super().__init__(name, waveband, "M87SED_EHTMWL2018", D_Mpc, MBH_MSUN, theta_view_deg)


