import pandas as pd
from importlib import resources


def get_lightcurves_VHE_flare_2018():

    path = resources.files("M87data.data").joinpath(f"lightcurves/Algaba24_tablec14.dat")
    df = pd.read_csv(path, delimiter=",", header=None)
    df.columns = ["instrument","MJD_start", "MJD_end", "Obs time", "zenith angle", "significance", "flux E>350GeV", "flux error upper", "flux error lower"]
    df["MJD_mid"] = (df["MJD_start"] + df['MJD_end'])/2
    dfh = df.loc[df["instrument"]=="HESS"].copy() 
    dfm = df.loc[df["instrument"]=="MAGIC"].copy() 
    dfv = df.loc[df["instrument"]=="VERITAS"].copy() 
    return df, dfh, dfm, dfv

def add2ax_lightcurves_VHE_flare_2018(
        ax, mjd_ref = 58229, 
        args_hess={"ls": "", "marker": "."}, 
        args_magic={"ls": "", "marker": "."},
        args_veritas={"ls": "", "marker": "."}):

    df, dfh, dfm, dfv = get_lightcurves_VHE_flare_2018()
    conv = 1e-12
    ax.errorbar(dfh["MJD_mid"]-mjd_ref, conv*dfh['flux E>350GeV'], 
                yerr=[conv*dfh["flux error lower"], conv*dfh["flux error upper"]], 
                **args_hess)

    ax.errorbar(dfm["MJD_mid"]-mjd_ref, conv*dfm['flux E>350GeV'], 
                yerr=[conv*dfm["flux error lower"], conv*dfm["flux error upper"]], 
                **args_magic)

    ax.errorbar(dfv["MJD_mid"]-mjd_ref, conv*dfv['flux E>350GeV'], 
                yerr=[conv*dfv["flux error lower"], conv*dfv["flux error upper"]], 
                **args_veritas)
    return ax

