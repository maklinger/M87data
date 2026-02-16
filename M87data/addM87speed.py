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


def get_jet_speed_data(D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17):
    
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    vel_files = [
        "EVN-1.6GHz-A14-app-z.d",
            "VLBI-1.6GHz-R89-app-z.d",
            "VLBA-1.7GHz-C07-app-z.d",
            "VLBA-1.7GHz+EVN-5GHz-G12-app-z.d",
            "VLBA-15GHz-K04-app-z.d",
            "VLBA-15GHz-K07-app-z.d",
        "VLBA-15GHz-M16-app-z.d",
        "VLA-15GHz-B95-app-z.d",
            "VLBA-43GHz-A09-app-z.d",
            "VLBA-43GHz-L07-app-z.d",
        "VLBA-43GHz-M16-app-z.d",
        "HST-B99-app-z.d",
        "HST-M13-app-z.d",
    ]
    vel_files_4 = [
        "KaVA-22GHz-H17-app-z.d",
        "HSA-86GHz-H16-app-z.d"
    ]

    fac_ang_dist_z = distance/rg * arcsec2rad / np.sin(inclination)

    df = pd.DataFrame(columns=["z [rg]", "betaGamma", "sigma_betaGamma"])

    for i, filename in enumerate(vel_files):
        path = resources.files("M87data.data").joinpath(f"speed/{filename}")
        data = np.genfromtxt(path).T
        ang_dist = data[0]
        z = ang_dist*fac_ang_dist_z
        beta_app = data[1]
        sigma_beta_app = data[2]
        beta = beta_app/(beta_app*np.cos(inclination) + np.sin(inclination))
        bG = beta/np.sqrt(1-beta**2)
        fac_betaapp_beta = beta/beta_app * np.abs(1-beta*np.cos(inclination))
        fac_bG_beta = (1-beta**2)**(-3/2)
        sigma_bG = fac_betaapp_beta * fac_bG_beta * sigma_beta_app

        if isinstance(bG, np.ndarray):
            for j in range(len(bG)):
                df.loc[len(df)] = [z[j], bG[j], sigma_bG[j]]
        else:
            df.loc[len(df)] = [z, bG, sigma_bG]
        

    for i, filename in enumerate(vel_files_4):
        path = resources.files("M87data.data").joinpath(f"speed/{filename}")
        data = np.genfromtxt(path).T
        ang_dist = data[0]
        z = ang_dist*fac_ang_dist_z
        beta_app = data[2]
        sigma_beta_app = data[3]
        beta = beta_app/(beta_app*np.cos(inclination) + np.sin(inclination))
        bG = beta/np.sqrt(1-beta**2)
        fac_betaapp_beta = beta/beta_app * np.abs(1-beta*np.cos(inclination))
        fac_bG_beta = (1-beta**2)**(-3/2)
        sigma_bG = fac_betaapp_beta * fac_bG_beta * sigma_beta_app

        if isinstance(bG, np.ndarray):
            for j in range(len(bG)):
                df.loc[len(df)] = [z[j], bG[j], sigma_bG[j]]
        else:
            df.loc[len(df)] = [z, bG, sigma_bG]
    
    return df


    

DEFAULT_JET_SPEED_COLS = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255',
       '#44AA99', '#999933', '#6699CC', '#004488', '#EECC66', '#994455',
       '#997700', '#EE99AA', '#000000', '#555555']

def add_jet_speed(ax, D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17, 
                    cols=DEFAULT_JET_SPEED_COLS, legalpha=0, number_leg_1=7):
    
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    vel_files = [
        "EVN-1.6GHz-A14-app-z.d",
        "VLBI-1.6GHz-R89-app-z.d",
        "VLBA-1.7GHz-C07-app-z.d",
        "VLBA-1.7GHz+EVN-5GHz-G12-app-z.d",
        "VLBA-15GHz-K04-app-z.d",
        "VLBA-15GHz-K07-app-z.d",
        "VLBA-15GHz-M16-app-z.d",
        "VLA-15GHz-B95-app-z.d",
        "VLBA-43GHz-A09-app-z.d",
        "VLBA-43GHz-L07-app-z.d",
        "VLBA-43GHz-M16-app-z.d",
        "HST-B99-app-z.d",
        "HST-M13-app-z.d",
        "KaVA-22+43GHz-P19-app-z.d"
    ]
    vel_files_4 = [
        "KaVA-22GHz-H17-app-z.d",
        "HSA-86GHz-H16-app-z.d"
    ]

    fac_ang_dist_z = distance/rg * arcsec2rad / np.sin(inclination)
    ax.set_xlabel(xlabel = r"jet axial distance $z$ [$r_g$]")
    ax.set_ylabel(r"jet speed $\beta \Gamma$")
    ax.set_ylim(1e-3, 1e2)
    ax.set_xlim(0.3, 1e8)

    handles_group1 = []
    handles_group2 = []

    for i, filename in enumerate(vel_files):
        path = resources.files("M87data.data").joinpath(f"speed/{filename}")
        data = np.genfromtxt(path).T
        ang_dist = data[0]
        if filename=="KaVA-22+43GHz-P19-app-z.d":
            ang_dist*=1e-3
        beta_app = data[1]
        sigma_beta_app = data[2]
        beta = beta_app/(beta_app*np.cos(inclination) + np.sin(inclination))
        bG = beta/np.sqrt(1-beta**2)
        fac_betaapp_beta = beta/beta_app * np.abs(1-beta*np.cos(inclination))
        fac_bG_beta = (1-beta**2)**(-3/2)
        sigma_bG = fac_betaapp_beta * fac_bG_beta * sigma_beta_app

        label = filename.split("-app")[0].replace("-", " ")
        parts = label.split()
        # wrap last part in parentheses
        parts[-1] = f"({parts[-1]})"
        label = " ".join(parts)
        if np.sum(sigma_beta_app) >0:
            h = ax.errorbar(ang_dist*fac_ang_dist_z, bG, yerr=sigma_bG, 
                        ls="", marker=".", alpha=0.5, color=cols[i],
                        label=label.replace("+", " &\n"))
        else:
            h, = ax.loglog([ang_dist*fac_ang_dist_z], [bG], 
                        ls="", marker="s", alpha=0.5, color=cols[i],
                        label=label)
        if len(handles_group1) < number_leg_1:
            handles_group1.append(h)
        else:
            handles_group2.append(h)

    for i, filename in enumerate(vel_files_4):
        path = resources.files("M87data.data").joinpath(f"speed/{filename}")
        data = np.genfromtxt(path).T
        ang_dist = data[0]
        sigma_ang_dist = data[1]
        beta_app = data[2]
        sigma_beta_app = data[3]
        beta = beta_app/(beta_app*np.cos(inclination) + np.sin(inclination))
        bG = beta/np.sqrt(1-beta**2)
        fac_betaapp_beta = beta/beta_app * np.abs(1-beta*np.cos(inclination))
        fac_bG_beta = (1-beta**2)**(-3/2)
        sigma_bG = fac_betaapp_beta * fac_bG_beta * sigma_beta_app

        label = filename.split("-app")[0].replace("-", " ")
        parts = label.split()
        # wrap last part in parentheses
        parts[-1] = f"({parts[-1]})"
        label = " ".join(parts)

        h = ax.errorbar(ang_dist*fac_ang_dist_z, bG, yerr=sigma_bG, xerr=sigma_ang_dist*fac_ang_dist_z,
                    ls="", marker=".", alpha=0.5, color=cols[i + len(vel_files)],
                    label=label)
        
        if len(handles_group1) < number_leg_1:
            handles_group1.append(h)
        else:
            handles_group2.append(h)

    leg1 = ax.legend(handles=handles_group1, loc="upper left", 
                     framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=handles_group2, loc="lower left", ncol=1, 
                     framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
        
    # ax.legend(loc="center left", ncol=1, framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    return ax


