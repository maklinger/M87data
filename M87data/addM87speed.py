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


def load_core_data():
    datasets = []
    current_label = None
    block_data = []

    path = (
        resources.files("M87data.data")
        .joinpath("collimation/M87coresizes.txt")
    )

    with path.open("r") as f:
        for line in f:
            line = line.strip()

            # Skip empty lines between blocks
            if line == "":
                # If we finished a block, store it
                if current_label is not None and block_data:
                    arr = np.array(block_data, dtype=float)
                    datasets.append({
                        "label": current_label,
                        "data": {
                            "col1": arr[:,0],
                            "col2": arr[:,1],
                            "col3": arr[:,2],
                            "col4": arr[:,3],
                        }
                    })
                    block_data = []
                continue

            # Comment line → start new block label
            if line.startswith("#"):
                current_label = line[1:].strip()   # remove "#"
                continue

            # Numeric data line → part of a block
            parts = line.split()
            if len(parts) == 4:
                block_data.append([float(x) for x in parts])

        # Append last block if file doesn't end with blank line
        if current_label is not None and block_data:
            arr = np.array(block_data, dtype=float)
            datasets.append({
                "label": current_label,
                "data": {
                    "col1": arr[:,0],
                    "col2": arr[:,1],
                    "col3": arr[:,2],
                    "col4": arr[:,3],
                }
            })

    return datasets

def get_collimation_data(D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17):
    
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    fac_ang_dist_z = distance/rg * marcsec2rad / np.sin(inclination)
    fac_ang_dist_r = distance/rg * marcsec2rad 

    df = pd.DataFrame(columns=["z [rg]", "r [rg]", "sigma_r [rg]"])

    core_datasets = load_core_data()
    for ds in core_datasets[::-1]:
        for i in range(len(ds["data"]["col1"])):
            df.loc[len(df)] = [ds["data"]["col1"][i]*fac_ang_dist_z, 
                               ds["data"]["col2"][i]*fac_ang_dist_r/2, 
                               ds["data"]["col4"][i]*fac_ang_dist_r/2]

    # these files contain the width = 2 * radius
    w_files = [
        "EVN-1.6GHz-AN12-w-z.d",
        "VLBA-2.3GHz-H13-w-z.d",
        "VLBA-5.0GHz-H13-w-z.d",
        "VLBA-8.4GHz-H13-w-z.d",
        "VLBA-15GHz-AN12-w-z.d",
        "VLBA-15GHz-H13-w-z.d",
        "VLBA-22GHz-H13-w-z.d",
        "VLBA-43GHz-AN12-w-z.d",
        "VLBA-43GHz-H13-w-z.d",
        "HSA-86GHz-H16-w-z.d"
    ]


    for i, filename in enumerate(w_files):
        path = resources.files("M87data.data").joinpath(f"collimation/{filename}")
        a, r, sr = np.genfromtxt(path).T
        for j in range(len(a)):
            df.loc[len(df)] = [a[j]*fac_ang_dist_z, r[j]*fac_ang_dist_r/2, sr[j]*fac_ang_dist_r/2]

    # these files contain the radius = width/2
    r_files = [
        "MERLIN-1.8GHz-AN12-r-z.d"
    ]

    for i, filename in enumerate(r_files):
        path = resources.files("M87data.data").joinpath(f"collimation/{filename}")
        a, r, sr = np.genfromtxt(path).T
        for j in range(len(a)):
            df.loc[len(df)] = [a[j]*fac_ang_dist_z, r[j]*fac_ang_dist_r, sr[j]*fac_ang_dist_r]

    return df






def add_collimation(ax, D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17, 
                    cmap=plt.get_cmap("viridis"), cmap_core=plt.get_cmap("magma"),
                    plot_BH=True, legalpha=0):
    
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    w_files = [
        "EVN-1.6GHz-AN12-w-z.d",
        "VLBA-2.3GHz-H13-w-z.d",
        "VLBA-5.0GHz-H13-w-z.d",
        "VLBA-8.4GHz-H13-w-z.d",
        "VLBA-15GHz-AN12-w-z.d",
        "VLBA-15GHz-H13-w-z.d",
        "VLBA-22GHz-H13-w-z.d",
        "VLBA-43GHz-AN12-w-z.d",
        "VLBA-43GHz-H13-w-z.d",
        "HSA-86GHz-H16-w-z.d"
    ]
    r_files = [
        "MERLIN-1.8GHz-AN12-r-z.d"
    ]

    cols = cmap(np.linspace(0,1, len(w_files) + len(r_files)))
    cols2 = cmap_core(np.linspace(0,1, 9))

    fac_ang_dist_z = distance/rg * marcsec2rad / np.sin(inclination)
    fac_ang_dist_r = distance/rg * marcsec2rad 
    ax.set_xlabel(xlabel = r"jet axial distance $z$ [$r_g$]")
    ax.set_ylabel(ylabel = r"jet radius $r$ [$r_g$]")
    ax.set_ylim(0.3, 1e6)
    ax.set_xlim(0.3, 1e8)
    
    if plot_BH:
        blackhole = np.logspace(-1, 0, 300)
        ax.fill_between(2 * blackhole, 0* blackhole, 2 * np.sqrt(1 - blackhole**2), color="k")

    handles_group1 = []
    handles_group2 = []

    core_datasets = load_core_data()
    icol = 0
    for ds in core_datasets[::-1]:
        h = ax.errorbar(ds["data"]["col1"]*fac_ang_dist_z, ds["data"]["col2"]*fac_ang_dist_r/2, 
                    xerr=ds["data"]["col3"]*fac_ang_dist_z, yerr=ds["data"]["col4"]*fac_ang_dist_r/2, 
                    ls="", marker="s", alpha=0.5, color=cols2[icol],
                    label="".join(ds["label"].split("core ")))
        icol += 1
        handles_group1.append(h)

    for i, filename in enumerate(r_files):
        path = resources.files("M87data.data").joinpath(f"collimation/{filename}")
        ang_dist, ang_radius, ang_rad_sigma = np.genfromtxt(path).T
        label = filename.split("-r")[0].replace("-", " ")
        parts = label.split()
        # wrap last part in parentheses
        parts[-1] = f"({parts[-1]})"
        label = " ".join(parts)
        h = ax.errorbar(ang_dist*fac_ang_dist_z, ang_radius*fac_ang_dist_r, yerr=ang_rad_sigma*fac_ang_dist_r, 
                    ls="", marker=".", alpha=0.5, color=cols[i],
                    label=label)
        handles_group2.append(h)

    for i, filename in enumerate(w_files):
        path = resources.files("M87data.data").joinpath(f"collimation/{filename}")
        ang_dist, ang_radius, ang_rad_sigma = np.genfromtxt(path).T
        label = filename.split("-w")[0].replace("-", " ")
        parts = label.split()
        # wrap last part in parentheses
        parts[-1] = f"({parts[-1]})"
        label = " ".join(parts)
        h = ax.errorbar(ang_dist*fac_ang_dist_z, ang_radius*fac_ang_dist_r/2, yerr=ang_rad_sigma*fac_ang_dist_r/2, 
                    ls="", marker=".", alpha=0.5, color=cols[i+len(r_files)],
                    label=label)
        handles_group2.append(h)
        
    leg1 = ax.legend(handles=handles_group1, title="Core size", loc="upper left", 
                     framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=handles_group2, loc="lower right", ncol=1, 
                     framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
        
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    return ax


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
                    cols=DEFAULT_JET_SPEED_COLS, legalpha=0):
    
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
        if np.sum(sigma_beta_app) >0:
            ax.errorbar(ang_dist*fac_ang_dist_z, bG, yerr=sigma_bG, 
                        ls="", marker=".", alpha=0.5, color=cols[i],
                        label=filename.split("-app")[0].replace("-", " ").replace("+", " &\n"))
        else:
            ax.loglog(ang_dist*fac_ang_dist_z, bG, 
                        ls="", marker="s", alpha=0.5, color=cols[i],
                        label=filename.split("-app")[0].replace("-", " "))

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
        ax.errorbar(ang_dist*fac_ang_dist_z, bG, yerr=sigma_bG, xerr=sigma_ang_dist*fac_ang_dist_z,
                    ls="", marker=".", alpha=0.5, color=cols[i + len(vel_files)],
                    label=filename.split("-app")[0].replace("-", " "))

    ax.legend(loc="center left", ncol=1, framealpha=legalpha, labelspacing=0.1, handletextpad=-0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    return ax


