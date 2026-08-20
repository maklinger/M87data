import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pandas as pd
from importlib import resources

deg2rad = np.pi/180
arcmin2rad = deg2rad/60
arcsec2rad = arcmin2rad/60
marcsec2rad = arcsec2rad/1000
muarcsec2rad = marcsec2rad/1000

pc = 3.086e18 # cm
msun = 1.989e33  # g
G = 6.67259e-8 # cgs
c = 2.99792458e10 # cm/s
h = 6.62e-27 # ergs

eV2erg = 1.60218e-12
erg2eV = 1/eV2erg


def get_SED_data(dataset="M87SED_EHTMWL2018",
        D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17):
        
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    path = resources.files("M87data.data").joinpath(f"SED/{dataset}.csv")
    df = pd.read_csv(path, delimiter=",")
    df["radius [rg]"] = df['Angular_Scale_arcsec'] * distance/rg * arcsec2rad / np.sin(inclination) / 2
    df["energy [eV]"] = df['Frequency_Hz'] * h * erg2eV
    df["UL"] = df["sigma_nuFnu_1e-12ergscm2"]==0
    # df["sigma_nuFnu_1e-12ergscm2"][df["UL"]==True] = df["nuFnu_1e-12ergscm2"][df["UL"]==True] *0.2
    df.loc[df["UL"], "sigma_nuFnu_1e-12ergscm2"] = (
        df.loc[df["UL"], "nuFnu_1e-12ergscm2"] * 0.3
    )
    df_data = df.loc[~df["UL"]].copy() 
    df_UL = df.loc[df["UL"]].copy()
    df_VM = df_data[(df_data["Instrument"]=="VERITAS") | (df_data["Instrument"]=="MAGIC")].copy()
    if dataset=="M87SED_EHTMWL2018":
        df_data = df_data.drop(list(np.arange(61, 65, dtype=int)))
        df_data = df_data.drop(list(np.arange(66, 69, dtype=int)))

    return df, df_data, df_VM, df_UL

def add_SED(ax, dataset="M87SED_EHTMWL2018", angular_scale_color=True, 
            cmap = plt.get_cmap("viridis"), D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17,
            convert_to_rg=True, colorbar=True, fixed_color="k", marker=".", ms=6, alphaUL=1,
            frequencymarkers=True, absorbed=True, ymin=1e-20, ymax=1e-10, label=None,
            alpha=1, cb_shrink=0.75, cb_pad=0.01):
    
    df, df_data, df_VM, df_UL = get_SED_data(dataset, D_Mpc, MBH_MSUN, theta_view_deg)

    if angular_scale_color:
        if convert_to_rg:
            distance = D_Mpc* 1e6 * pc
            MBH = MBH_MSUN * msun
            rg = MBH*G/c**2
            inclination = np.deg2rad(theta_view_deg) 
            fac_ang_dist_z = distance/rg * arcsec2rad / np.sin(inclination)
            col_label= r'aperture $\log_{10}(z/r_g)$'
            norm = Normalize(vmin=1, vmax=9)
        else:
            fac_ang_dist_z = 1
            col_label = r'aperture $\log_{10}(z/1\,arcsec)$'
            norm = Normalize(vmin=np.log10(df_data['Angular_Scale_arcsec'].min()), 
                            vmax=np.log10(df_data['Angular_Scale_arcsec'].max()))
        if dataset=="M87SED_EHTMWL2018":
            ax.errorbar(df_VM["Frequency_Hz"]*h*erg2eV, df_VM["nuFnu_1e-12ergscm2"]*1e-12, 
                yerr=df_VM["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=".", 
                c='grey', ms=ms, alpha=0.3)
        for i in range(len(df_data)):
            ax.errorbar([df_data["Frequency_Hz"].iat[i]*h*erg2eV], [df_data["nuFnu_1e-12ergscm2"].iat[i]*1e-12], 
                        yerr=[df_data["sigma_nuFnu_1e-12ergscm2"].iat[i]*1e-12], ls="", marker=marker, ms=ms,
                        color=cmap(norm(np.log10(df_data['Angular_Scale_arcsec'].iat[i]*fac_ang_dist_z))),
                        label=label, alpha=alpha,
                        zorder=10-np.log10(df_data['Angular_Scale_arcsec'].iat[i]*fac_ang_dist_z))
            label=None
        for i in range(len(df_UL)):
            ax.errorbar([df_UL["Frequency_Hz"].iat[i]*h*erg2eV], [df_UL["nuFnu_1e-12ergscm2"].iat[i]*1e-12], 
                        yerr=[df_UL["sigma_nuFnu_1e-12ergscm2"].iat[i]*1e-12], uplims=True, ls="", marker="_", 
                        ms=ms,
                        color=cmap(norm(np.log10(df_UL['Angular_Scale_arcsec'].iat[i]*fac_ang_dist_z))), alpha=alphaUL)
        
        if colorbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, label=col_label, shrink=cb_shrink, pad=cb_pad)
    else:
        ax.errorbar(df_data["Frequency_Hz"]*h*erg2eV, df_data["nuFnu_1e-12ergscm2"]*1e-12, 
                    yerr=df_data["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=marker, 
                    c=fixed_color, ms=ms, label=label, alpha=alpha)
        if dataset=="M87SED_EHTMWL2018":
            ax.errorbar(df_VM["Frequency_Hz"]*h*erg2eV, df_VM["nuFnu_1e-12ergscm2"]*1e-12, ms=ms,
                yerr=df_VM["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=".", c='grey', alpha=0.3, label=label)
        ax.errorbar(df_UL["Frequency_Hz"]*h*erg2eV, df_UL["nuFnu_1e-12ergscm2"]*1e-12, ms=ms,
                    yerr=df_UL["sigma_nuFnu_1e-12ergscm2"]*1e-12,uplims=True, ls="", marker="_", c="grey", alpha=alphaUL)

    if frequencymarkers:
        for i in [7,8, 10, 11, 13, 14, 16, 17]:
            ax.plot(h * erg2eV *10**i, 1.1*ymax, marker="|", c="grey")
        ax.plot(h * erg2eV *1e9, 0.9 * ymax, marker="|", c="grey")
        ax.plot(h * erg2eV *1e12, 0.9 * ymax, marker="|", c="grey")
        ax.plot(h * erg2eV *1e15, 0.9 * ymax, marker="|", c="grey")
        ax.plot(h * erg2eV *1e18, 0.9 * ymax, marker="|", c="grey")
        ax.text(h * erg2eV *1e9, ymax, "1GHz", va="bottom", ha="center", c="grey")
        ax.text(h * erg2eV *1e12, ymax, "1THz", va="bottom", ha="center", c="grey")

        ax.plot(c * h * erg2eV, 0.9 * ymax, marker="|", c="lightgrey")
        ax.text(c * h * erg2eV, 0.6 * ymax, "cm", va="top", ha="center", c="lightgrey")
        ax.text(c /0.1 * h * erg2eV, 0.6 * ymax, "mm", va="top", ha="center", c="lightgrey")
        ax.plot(c /0.1 * h * erg2eV, 0.9 * ymax, marker="|", c="lightgrey")
        ax.plot(c /0.01 * h * erg2eV, 0.9 * ymax, marker="|", c="lightgrey")
        ax.plot(c /0.001 * h * erg2eV, 0.9 * ymax, marker="|", c="lightgrey")
        ax.text(c /1e-4 * h * erg2eV, 0.6 * ymax, r"$\mu$m", va="top", ha="center", c="lightgrey")
        ax.plot(c /1e-4 * h * erg2eV, 0.9 * ymax, marker="|", c="lightgrey")

    if absorbed:
        ax.axvspan(1e1, 1e3, color="k", alpha=0.3, ls=":")
        ax.text(1e2, 1e1*ymin, "dust&\nphotoel.\nabs.", ha="center", va="bottom", color="grey")
        ax.axvspan(1e13, 1e14, color="k", alpha=0.3, ls=":")
        ax.text(3e13, 1e1*ymin, "EBL.\nabs.", ha="center", va="bottom", color="grey")

    # ax.set_title("M87 VHE flare (2018)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(10**np.arange(-6., 14, 3), minor=0)
    ax.set_xticks(10**np.arange(-6., 14.5, 1), minor=1)
    ax.set_ylabel(r"$\nu F_\nu$ [erg/cm²s]")
    ax.set_xlabel("Energy [eV]", fontsize=14)

    ax.set_xlim(1e-6, 1e14)
    ax.set_ylim(ymin, ymax)
    # ax.legend(framealpha=1, ncol=1, loc="upper right")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    if angular_scale_color:
        return ax, norm
    else:
        return ax
    

def add_LHAASO(ax, color="k", marker="."):

    path = resources.files("M87data.data").joinpath(f"SED/LHAASO22.csv")
    E, EFE = np.genfromtxt(path, skip_header=1, delimiter=",").T
    path = resources.files("M87data.data").joinpath(f"SED/LHAASO22_min.csv")
    E_min, EFE_min = np.genfromtxt(path, skip_header=1, delimiter=",").T
    path = resources.files("M87data.data").joinpath(f"SED/LHAASO22_max.csv")
    E_max, EFE_max = np.genfromtxt(path, skip_header=1, delimiter=",").T
    
    ax.errorbar(E, EFE, yerr=[EFE-EFE_min, EFE_max-EFE], color=color, marker=marker, ls="")
    return ax