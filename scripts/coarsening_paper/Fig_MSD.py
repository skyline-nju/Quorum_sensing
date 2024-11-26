import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
sys.path.append("../MSD/")


from add_line import add_line
from scripts.coarsening_paper.plot_reduced_bubble import show_R

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def untangle(x, Lx):
    x_new = np.copy(x)
    phase_shift = 0 
    for i in range(1, x.size):
        dx = x[i] - x[i-1]
        if dx < -Lx/2:
            phase_shift += Lx
        elif dx > Lx/2:
            phase_shift -= Lx
        x_new[i] += phase_shift
    return x_new


def read_traj(f_traj, Lx, Ly, n_cut=100, which="B"):
    with open(f_traj, "r") as f:
        lines = f.readlines()[n_cut:]
        t, x, y = np.zeros((3, len(lines)))
        for i, line in enumerate(lines):
            s = line.rstrip("\n").split("\t")
            t[i] = float(s[0]) / 3
            if which == "A":
                x[i] = float(s[1])
                y[i] = float(s[2])
            else:
                x[i] = float(s[3])
                y[i] = float(s[4])
        x += Lx/2
        y += Ly/2
        x[x>=Lx/2] -= Lx
        y[y>=Ly/2] -= Ly
        x = untangle(x, Lx)
        y = untangle(y, Ly)
    return t, x, y


def show_MSD_CCB_gas(ax, which="CCB+G"):
    if which == "CCB+G":
        L_arr = np.array([240, 320, 400, 480, 640])
        phi_A = 4.5
        phi_B = 1.53
    elif which == "LB+G":
        L_arr = np.array([160, 240, 320, 400, 480])
        phi_A = 0
        phi_B = 6
    for i, L in enumerate(L_arr):
        if which == "CCB+G":
            if L != 320:
                seed = 2000
            else:
                seed = 2001
        else:
            seed = 2000
        folder = f"{root_tahmineh}/QS_MPI_noself/MC_L{L:d}_r10_10_10_h0.333"
        fin = f"{folder}/mc_L{L:d}_{L:d}_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.dat"
        n = 1000
        ncut = 500
        t, x, y = read_traj(fin, L_arr[i], L_arr[i], ncut)
        t_new = np.arange(1, n+1) * 100
        rr_new = np.zeros(n)
        count = 0
        for j in range(t.size - n):
            x0, y0 = x[j], y[j]
            rr_new += (x[j+1:j+n+1] - x0)**2 + (y[j+1: j+n+1] - y0)**2
            count += 1
        rr_new /= count
        ax.plot(t_new, rr_new, lw=2, label=r"$L=%d$" % L)
    ax.legend(fontsize="large")
    ax.set_xscale("log")
    ax.set_yscale("log")
    add_line(ax, 0.3, 0.1, 0.95, 1, label="slope=1", xl=0.60, yl=0.2, fontsize="large")


def D_vs_L_varied_composition(ax):
    L_arr = [240, 320, 400, 480, 640]
    MSD_arr = np.array([230.9, 122.7, 73.85, 52.19, 29.24]) / 4e4
    ax.plot(L_arr, MSD_arr, "-s", label="CCB+G")

    # L_arr = [160, 240, 320, 480]
    # MSD_arr = np.array([21.47, 11.05, 6.8, 3.08]) / 4e4
    # ax.plot(L_arr, MSD_arr, "->", label=r"$\bar{\rho}_A=3, \bar{\rho}_B=6$")

    L_arr = [160, 240, 320, 400, 480]
    MSD_arr = np.array([21, 10.67, 6.15, 3.98, 2.86]) / 4e4
    ax.plot(L_arr, MSD_arr, "->", label=r"LB+G")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    add_line(ax, 0.05, 1, 0.99, -2, label="slope=-2", xl=0.60, yl=0.4, fontsize="large")
    add_line(ax, 0.05, 0.6, 0.9, -2, label="slope=-2", xl=0.60, yl=0.2, fontsize="large")


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)

    ax1, ax2, ax3 = axes

    show_MSD_CCB_gas(ax1)
    # show_MSD_CCB_gas(ax1, which="LB+G")

    D_vs_L_varied_composition(ax2)

    show_R(ax3)

    fs = "xx-large"
    ax1.set_xlabel(r"$t$", fontsize=fs)
    ax1.set_ylabel(r"MSD", fontsize=fs)
    ax2.set_xlabel(r"$L$", fontsize=fs)
    ax2.set_ylabel(r"$D$", fontsize=fs)
    ax3.set_xlabel(r"$t$", fontsize=fs)
    ax3.set_ylabel(r"Typical domain size", fontsize=fs)
    ax1.set_title(r"(a) CCB+G", fontsize=fs)
    ax2.set_title(r"(b) Domain diffusivity", fontsize=fs)
    ax3.set_title(r"(c) Reduced bubble model", fontsize=fs)

    plt.show()
    # plt.savefig("fig/Fig5.pdf")
    plt.close()
