import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import sys
sys.path.append("../")

from structure_factor import get_h, get_Sq_t_rho, get_q1
from add_line import add_line

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

def plot_Rt_varied_L(ax, which="CCB+G"):
    if which == "CCB+G":
        phi_A, phi_B = 5.5, 3.85
    elif which == "CCB+LA":
        phi_A, phi_B = 22, 6
    elif which == "CCB+LA+G":
        phi_A, phi_B = 7.5, 3
    elif which == "LB+G":
        phi_A, phi_B = 3, 6
    folder  = "/home/yduan/code/Quorum_sensing/scripts/data"
    fins = [
        # f"{folder}/Sqt_rho/L320_320_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
        f"{folder}/Sqt_rho/L640_640_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
        f"{folder}/Sqt_rho/L1280_1280_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
        ]
    
    # if which == "CCB+G":
    #     Lx = [640, 1280]
    # else:
    #     Lx = [320, 640, 1280]
    Lx = [640, 1280]
    for j, fin in enumerate(fins):
        h = get_h(fin)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fin, dt=h)
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
        # ax1.plot(t, 2 * np.pi /q1_A, "o", fillstyle="none", label=r"$L=%g$" % Lx[j])
        if which == "CCB+G":
            mask = t < 2.5e6
            t = t[mask]
            q1_B = q1_B[mask]
        ax.plot(t, 2 * np.pi /q1_B, "o", fillstyle="none", label=r"$L=%g$" % Lx[j])
        ax.legend(loc="lower right", fontsize="large")

if __name__ == "__main__":
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.5), constrained_layout=True, sharey="row")


    snaps = ["fig/snaps/CCB_G.png", "fig/snaps/CCB_LA_G.png", "fig/snaps/CCB_LA.png", "fig/snaps/LB_G.png"]

    for i, ax in enumerate(axes[1]):
        im = mpimg.imread(snaps[i])
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
    
    modes = ["CCB+G", "CCB+LA+G", "CCB+LA", "LB+G"]
    titles = ["(a) CCB+G", "(b) CCB+LA+G", "(c) CCB+LA", "(d) LB+G"]
    for i, ax in enumerate(axes[0]):
        plot_Rt_varied_L(ax, which=modes[i])
        ax.set_xscale("log")
        ax.set_yscale("log")
        axes[0, i].set_title(titles[i], fontsize="xx-large")
        ax.set_xlabel(r"$t$", fontsize="xx-large")
    

    axes[0, 0].set_ylabel(r"Typical domain size", fontsize="xx-large")
    axes[1, 0].set_ylabel(r"Snapshots at $t=10^6$", fontsize="xx-large")
    

    add_line(axes[0, 0], 0.3, 0.5, 0.95, 0.13, label="0.13", xl=0.50, yl=0.4, fontsize="large")
    add_line(axes[0, 1], 0.3, 0.5, 0.95, 0.13, label="0.13", xl=0.50, yl=0.4, fontsize="large")
    add_line(axes[0, 2], 0.3, 0.45, 0.95, 0.13, label="0.13", xl=0.50, yl=0.4, fontsize="large")
    add_line(axes[0, 3], 0.3, 0.42, 0.95, 0.333, label="0.33", xl=0.50, yl=0.7, fontsize="large")

    # plt.show()
    plt.savefig("fig/Fig1.pdf", dpi=150)
    plt.close()