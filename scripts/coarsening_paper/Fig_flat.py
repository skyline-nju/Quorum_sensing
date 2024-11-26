import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import sys
sys.path.append("../")

from sk_1D import plot_R_t_1D

root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


def plot_spacing(ax1, ax2):
    Lx = 1280
    seed = 1000
    folder = "/home/yduan/code/Quorum_sensing/scripts/data"
    f_interface = f"{folder}/interface_locs/L{Lx}_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed}.npz"
    with np.load(f_interface, "r") as data:
        t, xs = data["t"], data["xs"]
        d1, d2 = np.zeros((2, t.size))
        w1, w2 = np.zeros((2, t.size))

        for i in range(t.size):
            d1[i] = xs[i][0] + Lx - xs[i][3]
            d2[i] = xs[i][2] - xs[i][1]

            w1[i] = xs[i][1] - xs[i][0]
            w2[i] = xs[i][3] - xs[i][2]
        
    
        ncut = 125
        ax1.plot(t[:ncut], d1[:ncut], "s", fillstyle="none", ms=3, c="tab:green")
        print(t.size)
        ax1.plot(t[:ncut], d2[:ncut], "o", fillstyle="none", ms=3, c="tab:pink")

        ax2.plot(t[:ncut], w1[:ncut], "s", fillstyle="none", ms=3, c="tab:green")
        ax2.plot(t[:ncut], w2[:ncut], "o", fillstyle="none", ms=3, c="tab:pink")


def Fig4():
    fig =  plt.figure(layout="constrained", figsize=(10, 8))

    subfigs = fig.subfigures(2, 1, height_ratios=[1, 8.], hspace=0.001)
    subfigsnest = subfigs[1].subfigures(1, 2, width_ratios=[1, 0.6])

    ax_snap = subfigs[0].subplots(1, 1)

    im = mpimg.imread("fig/snaps/flat_i170.png")
    ax_snap.imshow(im)
    ax_snap.set_xticks([])
    ax_snap.set_yticks([])
    ax_snap.set_title(r"$L_x=1280, L_y=80, t=6.5\times 10^6$", fontsize="x-large")


    ax_st = subfigsnest[0].subplots(2, 1, sharex=True)

    fname = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]
    
    rhoB = fields[:, 1] / 10

    rect = [0, 1280, 1, 6]
    im = ax_st[1].imshow(rhoB[:61], origin="lower", aspect="auto", extent=rect, vmax=3)

    rect = [0, 1280, 1e6, t[270]/3]
    im = ax_st[0].imshow(rhoB[61:271], origin="lower", aspect="auto", extent=rect, vmax=3)

    subfigsnest[0].colorbar(im, shrink=0.6, ax=ax_st, location="right", extend="max")
    subfigsnest[0].text(0.895, 0.84, r"$\frac{\langle \rho_B\rangle_{y}}{\rho_0}$", fontsize=24)

    ax_st[0].set_ylabel(r"$t$", fontsize="xx-large")
    ax_st[1].set_ylabel(r"$\log_{10} t$", fontsize="xx-large")
    ax_st[1].set_xlabel(r"$x$", fontsize="xx-large")

    ax_st[0].axhline(6.5e6, linestyle="dashed", color="tab:red", lw=2)


    ax_right = subfigsnest[1].subplots(3, 1, height_ratios=[1, 1, 1.8])

    plot_R_t_1D(ax_right[2])
    plot_spacing(ax_right[1], ax_right[0])
    ax_right[2].set_ylabel(r"mean domain size", fontsize="xx-large")
    ax_right[1].set_ylabel(r"domain spacing", fontsize="xx-large")
    ax_right[0].set_ylabel(r"domain size", fontsize="xx-large")
    ax_right[2].set_xlabel(r"$t$", fontsize="xx-large")


    ax_right[1].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)
    ax_right[2].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)
    ax_right[0].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)

    # ax_right[1].set_xlim(6.5e6)
    # ax_right[2].set_xlim(6.5e6)
    plt.show()
    # plt.savefig("fig/Fig4.pdf")
    plt.close()


if __name__ == "__main__":
    fig =  plt.figure(layout="constrained", figsize=(10, 6.5))

    subfigs = fig.subfigures(1, 2, width_ratios=[1, 0.5])
    subfigsnest = subfigs[0].subfigures(2, 1, height_ratios=[1, 6.5])

    ax_snap = subfigsnest[0].subplots(1, 1)

    im = mpimg.imread("fig/snaps/flat_i170.png")
    ax_snap.imshow(im)
    ax_snap.set_xticks([])
    ax_snap.set_yticks([])
    ax_snap.set_title(r"(a) $L_x=1280, L_y=80, t=6.5\times 10^6$", fontsize="xx-large")

    ax_st = subfigsnest[1].subplots(2, 1, sharex=True)

    fname = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]
    
    rhoB = fields[:, 1] / 10

    rect = [0, 1280, 1, 6]
    im = ax_st[1].imshow(rhoB[:61], origin="lower", aspect="auto", extent=rect, vmax=3)

    rect = [0, 1280, 1e6, t[270]/3]
    im = ax_st[0].imshow(rhoB[61:271], origin="lower", aspect="auto", extent=rect, vmax=3)

    subfigsnest[1].colorbar(im, shrink=0.6, ax=ax_st, location="right", extend="max")
    subfigsnest[1].text(0.895, 0.84, r"$\frac{\langle \rho_B\rangle_{y}}{\rho_0}$", fontsize=24)

    ax_st[0].set_ylabel(r"$t$", fontsize="xx-large")
    ax_st[1].set_ylabel(r"$\log_{10} t$", fontsize="xx-large")
    ax_st[1].set_xlabel(r"$x$", fontsize="xx-large")

    ax_st[0].axhline(6.5e6, linestyle="dashed", color="tab:red", lw=2)
    ax_st[0].set_title(r"(b)", fontsize="xx-large")


    ax_right = subfigs[1].subplots(3, 1, height_ratios=[1, 1, 1.8])

    ax_right[0].set_title(r"(c) domain size", fontsize="xx-large")
    ax_right[1].set_title(r"(d) domain spacing", fontsize="xx-large")
    ax_right[2].set_title(r"(e) mean domain size", fontsize="xx-large")

    plot_R_t_1D(ax_right[2])
    plot_spacing(ax_right[1], ax_right[0])
    # ax_right[2].set_ylabel(r"mean domain size", fontsize="xx-large")
    # ax_right[1].set_ylabel(r"domain spacing", fontsize="xx-large")
    # ax_right[0].set_ylabel(r"domain size", fontsize="xx-large")
    ax_right[2].set_xlabel(r"$t$", fontsize="xx-large")


    ax_right[1].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)
    ax_right[2].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)
    ax_right[0].axvline(6.5e6, linestyle="dashed", color="tab:red", lw=2)

    # ax_right[1].set_xlim(6.5e6)
    # ax_right[2].set_xlim(6.5e6)
    plt.show()
    # plt.savefig("fig/Fig4.pdf")
    plt.close()