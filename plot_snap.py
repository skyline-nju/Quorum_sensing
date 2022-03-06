import enum
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


def PD_pA_pB(pA_arr, pB_arr, file_pat):
    nrows = pB_arr.size
    ncols = pA_arr.size
    figsize = (15, 15)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, pB in enumerate(pB_arr[::-1]):
        for i, pA in enumerate(pA_arr):
            fin = file_pat % (pA, pB)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i])
            axes[j, i].axis("off")

    
    x = (np.arange(7) + 0.5) / 7
    for i, pA in enumerate(pA_arr):
        fig.text(x[i], 0.002, r"$\phi_A=%g$" % pA, fontsize="x-large")
    
    y = (np.arange(nrows) + 0.5) / nrows
    for j, pB in enumerate(pB_arr):
        fig.text(0.001, y[j], r"$\phi_B=%g$" % pB, rotation="vertical", fontsize="x-large")



    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1)
    # plt.show()
    plt.savefig("PD/a.png")
    plt.close()


if __name__ == "__main__":
    # fin = r"D:\tmp\L40_k0.7_r40_e-2_J0.5\L40_40_Dr0.100_k0.70_p30_40_r40_e-2.000_J0.500_-0.500_1002.gsd"
    # plot_one_panel(fin)

    pA_arr = np.arange(1, 8) * 10
    pB_arr = np.arange(1, 8) * 10
    file_pat = r"D:\tmp\L40_k0.7_r40_e-2_J0.9\L40_40_Dr0.100_k0.70_p%g_%g_r40_e-2.000_J0.900_-0.900_1002.gsd"
    PD_pA_pB(pA_arr, pB_arr, file_pat)