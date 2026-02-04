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

DEFAULT_MAG_COLS = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255',
       '#44AA99', '#999933', '#6699CC', '#004488', '#EECC66', '#994455',
       '#997700', '#EE99AA', '#000000', '#555555']


def add_mag_field(ax, cols=DEFAULT_MAG_COLS, legalpha=0):
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
