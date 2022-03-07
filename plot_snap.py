import os
from read_gsd import read_one_frame
import matplotlib.pyplot as plt
import numpy as np


def plot_one_panel(fname, i_frame=-1, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    
    snap = read_one_frame(fname, i_frame)
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    x = snap.particles.position[:, 0] + Lx / 2
    y = snap.particles.position[:, 1] + Ly / 2
    type_id = snap.particles.typeid

    mask0 = type_id == 0
    mask1 = type_id == 1
    ax.plot(x[mask0], y[mask0], ".", ms=0.05, c="b", alpha=1)
    ax.plot(x[mask1], y[mask1], ".", ms=0.05, c="r", alpha=0.7)

    ax.axis("scaled")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)

    if flag_show:
        plt.show()
        plt.close()


def PD_pA_pB(L, J, rho0, Dr=0.1, eta=-2):

    if L == 40:
        if rho0 == 40:
            pA_arr = np.arange(1, 8) * 10
            pB_arr = np.arange(1, 8) * 10
        elif rho0 == 20:
            pA_arr = np.arange(1, 8) * 5
            pB_arr = np.arange(1, 8) * 5
    elif L == 80:
        if rho0 == 40:
            pA_arr = np.array([20, 25, 30])
            pB_arr = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60])
        elif rho0 == 20:
            pA_arr = np.array([20, 25, 30]) / 2
            pB_arr = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60]) / 2
    
    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/L{L}_k0.7_r{rho0:g}_e{eta:g}_J{J:g}"
    if not os.path.exists(folder):
        folder = f"{prefix}/L{L}_k0.7_r{rho0:g}_e{eta:g}_J{J:.1f}"
    file_pat = f"{folder}/L{L}_{L}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_" \
        f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_1002.gsd"
    fout_name = f"fig/PD/L{L}_Dr{Dr:g}_r{rho0:g}_e{eta:g}_J{J:.2f}.png"
    nrows = pB_arr.size
    ncols = pA_arr.size
    figsize = (ncols * 2, nrows * 2)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, pB in enumerate(pB_arr[::-1]):
        for i, pA in enumerate(pA_arr):
            fin = file_pat % (pA, pB)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i])
            axes[j, i].axis("off")

    x = (np.arange(ncols) + 0.5) / ncols
    for i, pA in enumerate(pA_arr):
        fig.text(x[i], 0.002, r"$\phi_A=%g$" % pA, fontsize="x-large")
    y = (np.arange(nrows) + 0.5) / nrows
    for j, pB in enumerate(pB_arr):
        fig.text(0.001, y[j], r"$\phi_B=%g$" % pB, rotation="vertical",
            fontsize="x-large")

    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1)
    if fout_name is None:
        plt.show()
    else:
        plt.savefig(fout_name)
    plt.close()


if __name__ == "__main__":
    # fin = r"D:\tmp\L40_k0.7_r40_e-2_J0.5\L40_40_Dr0.100_k0.70_p30_40_r40_e-2.000_J0.500_-0.500_1002.gsd"
    # plot_one_panel(fin)
    
    PD_pA_pB(L=40, J=0.5, rho0=20)