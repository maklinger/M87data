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

    filename = resources.files("M87data.data").joinpath("collimation/M87coresizes.txt")
    with open(filename, "r") as f:
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


def get_jet_speed_data(D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17, 
                    files_path="/Users/marc/Work/M87/M87_data/velocity/"):
    
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 

    vel_files = [
        # "EVN-1.6GHz-A14-app-z.d",
            # "VLBI-1.6GHz-R89-app-z.d",
            # "VLBA-1.7GHz-C07-app-z.d",
            # "VLBA-1.7GHz+EVN-5GHz-G12-app-z.d",
            # "VLBA-15GHz-K04-app-z.d",
            # "VLBA-15GHz-K07-app-z.d",
        # "VLBA-15GHz-M16-app-z.d",
        # "VLA-15GHz-B95-app-z.d",
            # "VLBA-43GHz-A09-app-z.d",
            # "VLBA-43GHz-L07-app-z.d",
        "VLBA-43GHz-M16-app-z.d",
        # "HST-B99-app-z.d",
        # "HST-M13-app-z.d",
    ]
    vel_files_4 = [
        # "KaVA-22GHz-H17-app-z.d",
        # "HSA-86GHz-H16-app-z.d"
    ]

    fac_ang_dist_z = distance/rg * arcsec2rad / np.sin(inclination)

    df = pd.DataFrame(columns=["z [rg]", "betaGamma", "sigma_betaGamma"])

    for i, filename in enumerate(vel_files):
        data = np.genfromtxt(f"{files_path}{filename}").T
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
        data = np.genfromtxt(f"{files_path}{filename}").T
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
                    files_path="/Users/marc/Work/M87/M87_data/velocity/",
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
        data = np.genfromtxt(f"{files_path}{filename}").T
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
        data = np.genfromtxt(f"{files_path}{filename}").T
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



def add_mag_field(ax, cols=DEFAULT_JET_SPEED_COLS, legalpha=0):
    ax.errorbar([7.5], [np.sqrt(50*124)], xerr=[[7.5 - 4.07], [13.6-7.5]], yerr=[[124-np.sqrt(50*124)], [np.sqrt(50*124) - 50]], 
                marker="", capsize=3, ls="", color=cols[0],
                label="230 GHz core (K15)")
    ax.errorbar([18.9], [8.3], xerr=[[18.9 - 11.2], [31.5-18.9]], yerr=[[8.3-2.38], [22.45-8.3]], 
                marker="", capsize=3, ls="", color=cols[1],
                label="86 GHz core (H16)")

    ax.errorbar([36.28], [np.sqrt(1*15)], xerr=[[36.28 - 22.9], [56.7-36.28]], yerr=[[np.sqrt(1*15)-1], [15 - np.sqrt(1*15)]], 
                marker="", capsize=3, ls="", color=cols[2],
                label="43 GHz core (K14)")
    ax.errorbar([36.28], [2], xerr=[[36.28 - 22.9], [56.7-36.28]], yerr=[[2-0.574], [5.41 - 2]], 
                marker="", capsize=3, ls="", color=cols[3],
                label="43 GHz core (H12)")

    ax.errorbar([274.2], [0.2], xerr=[[274.2 - 210.1], [353.2 -274.2]], yerr=[[0.2-0.0574], [0.541 - 0.2]], 
                marker="", capsize=3, ls="", color=cols[4],
                label="5 GHz core (R96)")


    ax.errorbar([5.4], [5], yerr=[10], lolims=True,
                marker="", capsize=3, ls="", color=cols[5],
                label="EHT MWL21 (model 1)")

    ax.loglog([2], [1e3], marker="x", ls="", color=cols[6], label=r"$P_{BZ}$ constraint (B19)")


    ax.errorbar([10], [np.sqrt(30)], yerr=[[30**0.5 - 1], [30 - 30**0.5]], 
                marker="", capsize=3, ls="", color=cols[7],
                label="EHT results (EHTC 19, 21)")

    ax.errorbar([18.9], [np.sqrt(61*210)], yerr=[[np.sqrt(61*210) - 61], [210 - np.sqrt(61*210)]], 
                marker="", capsize=3, ls="", color=cols[8],
                label=r"$T_B$ of 86 GHz core (K18)")

    ax.loglog([36.28], [1.04], marker="d",ls="", color=cols[9], label="Core-shift (Z14)")

    ax.errorbar([18.5], [4.8], yerr=[[4.8 - 2.4], [7.2 - 4.8]], 
                marker="", capsize=3, ls="", color=cols[10],
                label=r"Core-shift (W21)")


    ax.loglog([764779.88, 764779.88], [1e-3, 6e-4], marker="s",ls="", color=cols[11], label="HST-1 X-ray \nvariability (H03, H09)")

    ax.set_xlabel(xlabel = r"jet axial distance $z$ [$r_g$]")
    ax.set_ylabel(ylabel = r"magnetic field $B$ [$G$]")
    ax.set_ylim(1e-4, 1e4)
    ax.set_xlim(0.3, 1e8)
    ax.legend(ncol=1, loc="upper right", framealpha=legalpha, labelspacing=0.2, handletextpad=-0.3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    return ax

def add_SED(ax, dataset="M87SED_EHTMWL2018",
            files_path="/Users/marc/Work/M87/M87_data/SED/", angular_scale_color=True, 
            cmap = plt.get_cmap("viridis"), D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17,
            convert_to_rg=True, colorbar=True, fixed_color="k", marker=".", alphaUL=1,
            frequencymarkers=True, absorbed=True, ymin=1e-20, ymax=1e-10):
    
    df, df_data, df_VM, df_UL = get_SED_data(dataset, files_path, D_Mpc, MBH_MSUN, theta_view_deg)

    if angular_scale_color:
        if convert_to_rg:
            distance = D_Mpc* 1e6 * pc
            MBH = MBH_MSUN * msun
            rg = MBH*G/c**2
            inclination = np.deg2rad(theta_view_deg) 
            fac_ang_dist_z = distance/rg * arcsec2rad / np.sin(inclination)
            col_label= r'up to $\log_{10}(z/r_g)$'
            norm = Normalize(vmin=1, vmax=9)
        else:
            fac_ang_dist_z = 1
            col_label = r'up to $\log_{10}(z/1\,arcsec)$'
            norm = Normalize(vmin=np.log10(df_data['Angular_Scale_arcsec'].min()), 
                            vmax=np.log10(df_data['Angular_Scale_arcsec'].max()))
        if dataset=="M87SED_EHTMWL2018":
            ax.errorbar(df_VM["Frequency_Hz"]*h*erg2eV, df_VM["nuFnu_1e-12ergscm2"]*1e-12, 
                yerr=df_VM["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=".", c='grey', alpha=0.3)
        for i in range(len(df_data)):
            ax.errorbar([df_data["Frequency_Hz"].iat[i]*h*erg2eV], [df_data["nuFnu_1e-12ergscm2"].iat[i]*1e-12], 
                        yerr=[df_data["sigma_nuFnu_1e-12ergscm2"].iat[i]*1e-12], ls="", marker=marker, ms=8,
                        color=cmap(norm(np.log10(df_data['Angular_Scale_arcsec'].iat[i]*fac_ang_dist_z))))
        for i in range(len(df_UL)):
            ax.errorbar([df_UL["Frequency_Hz"].iat[i]*h*erg2eV], [df_UL["nuFnu_1e-12ergscm2"].iat[i]*1e-12], 
                        yerr=[df_UL["sigma_nuFnu_1e-12ergscm2"].iat[i]*1e-12], uplims=True, ls="", marker="_", 
                        color=cmap(norm(np.log10(df_UL['Angular_Scale_arcsec'].iat[i]*fac_ang_dist_z))), alpha=alphaUL)
        
        if colorbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, label=col_label, shrink=0.75, pad=0.01)
    else:
        ax.errorbar(df_data["Frequency_Hz"]*h*erg2eV, df_data["nuFnu_1e-12ergscm2"]*1e-12, 
                    yerr=df_data["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=marker, c=fixed_color)
        if dataset=="M87SED_EHTMWL2018":
            ax.errorbar(df_VM["Frequency_Hz"]*h*erg2eV, df_VM["nuFnu_1e-12ergscm2"]*1e-12, 
                yerr=df_VM["sigma_nuFnu_1e-12ergscm2"]*1e-12, ls="", marker=".", c='grey', alpha=0.3)
        ax.errorbar(df_UL["Frequency_Hz"]*h*erg2eV, df_UL["nuFnu_1e-12ergscm2"]*1e-12, 
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
    ax.set_ylabel("$\\nu F_\\nu$ (erg/cm2/s)")
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
    


def get_SED_data(dataset="M87SED_EHTMWL2018",
        files_path="/Users/marc/Work/M87/M87_data/SED/",
        D_Mpc=16.7, MBH_MSUN=6.6e9, theta_view_deg=17):
    distance = D_Mpc* 1e6 * pc
    MBH = MBH_MSUN * msun
    rg = MBH*G/c**2
    inclination = np.deg2rad(theta_view_deg) 
    df = pd.read_csv(f"{files_path}{dataset}.csv", delimiter=",")
    df["radius [rg]"] = df['Angular_Scale_arcsec'] * distance/rg * arcsec2rad / np.sin(inclination)
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