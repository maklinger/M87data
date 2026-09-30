from .M87_SED_Plugin import M87_SED_Plugin


class M87_2018_Plugin(M87_SED_Plugin):
    """EHT MWL 2018 dataset."""

    WAVEBANDS = {
        "radio":   {"min": 0,  "max": 13},
        "radio_VERA_EAVN_KAVA":   {"min": 0,  "max": 3},
        "radio_VLBA2.4":   {"min": 3,  "max": 3},
        "radio_VLBA4.3":   {"min": 4,  "max": 4},
        "radio_VLBI":   {"min": 5,  "max": 5},
        "radio_GMVA":   {"min": 6,  "max": 6},
        "radio_KVN":   {"min": 7,  "max": 8},
        "radio_ALMA93":   {"min": 9,  "max": 9},
        "radio_ALMA221":   {"min": 10,  "max": 10},
        "radio_SMA":   {"min": 11,  "max": 11},
        "radio_EHT":   {"min": 12,  "max": 12},
        "optical": {"min": 13, "max": 26},
        "xray":    {"min": 26, "max": 47},
        "gev":     {"min": 47, "max": 51},
        "tev":     {"min": 51, "max": 61},
    }

    def __init__(self, name, waveband, D_Mpc=16.8, MBH_MSUN=6.5e9, 
                 theta_view_deg=17, systematic_fraction=0.0):
        super().__init__(name, waveband, "M87SED_EHTMWL2018", D_Mpc, MBH_MSUN, 
                         theta_view_deg, systematic_fraction)


