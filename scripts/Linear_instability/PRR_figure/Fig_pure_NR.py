import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import image as mpimg
import os
import sys
sys.path.append("..")

from cal_omega_barV import get_bar_v_omega
from general_PD import find_long_instabi, find_short_instabi, find_contours


def get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, bar_vA=1, bar_vB=1, kappa=0.7,
                            extent=[0, 4, 0, 4], qmin=1e-5, qmax=1, Nq=200, resolution=1000, overwrite=False, ll=3./20, simple_gamma=False):
    if simple_gamma:
        folder = "../data/PD_pA_pB_SG/"
    elif ll == 0.:
        folder = "../data/PD_pA_pB_woS/"
    else:
        folder = "../data/PD_pA_pB/"
    fout = "%se%g_%g_%g_%g_D%g_%g_r%g_%g_v%g_%g_qm%g_r%g.npz" % (
        folder, etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, bar_vA, bar_vB, qmax, resolution
    )
    if not overwrite and os.path.exists(fout):
        with np.load(fout, "rb") as data:
            return data["extent"], data["state"], data["q_range"]
    else:
        nrows, ncols = resolution, resolution
        pA = np.linspace(extent[0], extent[1], ncols)
        pB = np.linspace(extent[2], extent[3], nrows)
        pA_2D, pB_2D = np.meshgrid(pA, pB)
        vA_0, vB_0, wAA, wAB, wBA, wBB = get_bar_v_omega(
            etaAA, etaAB, etaBA, etaBB, pA_2D, pB_2D, bar_rho_A, bar_rho_B, bar_vA, bar_vB, kappa)
        w1 = 1 + wAA
        w2 = 1 + wBB
        wc = wAB * wBA
        sigma_D = Dr_A / Dr_B
        sigma_v = vA_0 / vB_0
        Pe = vB_0 / Dr_B

        mask_LS, mask_LOI, mask_LSI = find_long_instabi(sigma_D, sigma_v, Pe, w1, w2, wc, q0=qmin, ll=ll, simple_gamma=simple_gamma)
        q0_a4, q0_Delta3 = np.zeros((2, nrows, ncols))
        for row in range(nrows):
            for col in range(ncols):
                if mask_LS[row, col]:
                    q0_a4[row, col], q0_Delta3[row, col] = find_short_instabi(
                        sigma_D, sigma_v[row, col], Pe[row, col], w1[row, col], w2[row, col], wc[row, col], qmin=qmin, qmax=qmax, Nq=Nq, ll=ll, simple_gamma=simple_gamma)
        state = np.zeros((nrows, ncols), dtype=np.byte)
        state[mask_LSI] = 1
        state[mask_LOI] = 2
        state[q0_a4 > 0] = 3
        state[q0_Delta3 > 0] = 4
        q_range = np.array([qmin, qmax, Nq])
        np.savez_compressed(fout, extent=extent, state=state, q_range=q_range)
        return extent, state, q_range


def plot_state(state, extent, xlim=None, ylim=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    nrows, ncols = state.shape

    colors = np.zeros((nrows, ncols, 3), int)
    colors[state==0] = np.array([255,255,255])
    colors[state==1] = np.array([142,186,217])
    colors[state==2] = np.array([241,187,224])
    colors[state==3] = np.array([255,190,134])
    colors[state==4] = np.array([149,207,149])


    im = ax.imshow(colors, origin="lower", extent=extent)

    patches = [mpatches.Patch(color='tab:blue', label='LSI', alpha=0.5),
               mpatches.Patch(color='tab:pink', label='LOI',alpha=0.5),
               mpatches.Patch(color='tab:orange', label='SSI', alpha=0.5),
               mpatches.Patch(color='tab:green', label='SOI', alpha=0.5),
            ]
    ax.legend(handles=patches, loc="upper right", fontsize="large", frameon=False, labelspacing=0.25)

    # ax.plot(3/40, 56/40, "o")
    # ax.plot([5.92/40, 1.95/40], [27.13/40, 66.6/40], ":s")

    
    if xlim is None:
        xlim = [extent[0], extent[1]]
    if ylim is None:
        ylim = [extent[2], extent[3]]
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    if flag_show:
        plt.show()
        plt.close()


def plot_PD_composition(state, extent, xlim=None, ylim=None, ax=None, fill=True, legend_fs="x-large"):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    nrows, ncols = state.shape
    # ax.imshow(state, origin="lower", extent=extent)
    contours = find_contours(state)
    for contour in contours["LWS"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if fill:
            ax.fill(x, y, c="tab:blue", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:blue", lw=1)
    for contour in contours["LWO"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if fill:
            ax.fill(x, y, c="tab:pink", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:pink", lw=1)
    
    for contour in contours["SWS"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if fill:
            ax.fill(x, y, c="tab:orange", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:orange", lw=1)


    for contour in contours["SWO"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if fill:
            ax.fill(x, y, c="tab:green", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:green", lw=2.5, linestyle="dashed", label="SOI")
    ax.legend(fontsize=legend_fs)
    # # ax.plot(1, 1, "o", ms=2, c="k")


def panel_4_2_2():
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    fig = plt.figure(figsize=(8, 4))
    subfigs = fig.subfigures(1, 2, wspace=0.001, hspace=0.001, width_ratios=[1.05, 1])

    ax_left = subfigs[0].subplots(1, 1, gridspec_kw=dict(hspace=0, wspace=0, left=0.04, right=0.995, bottom=0.04, top=0.995))
    im = mpimg.imread(f"fig/L20_20_Dr0.1_r80_e0_J1.00_-1.00.jpeg")
    extent = [0.5-0.125, 2+0.125, 0.25-0.125, 1.75+0.125]
    ax_left.imshow(im, extent=extent)
    ax_left.set_xticks([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2])
    ax_left.set_yticks([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75])
    for tick in ax_left.get_yticklabels():
        tick.set_rotation(90)

    Dr_A = Dr_B = 0.1
    etaAA = 0
    etaBB = 0

    etaAB = 1
    etaBA = -1
    bar_rho_A = 1
    bar_rho_B = 1
    bar_vA = 1
    bar_vB = 1
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, bar_vA, bar_vB, qmax=2.5, Nq=400, resolution=2000)
    plot_PD_composition(state, extent, ax=ax_left, fill=False)

    ax_left.text(0.82, 0.02, r"$\bar{\rho}_A/\rho_0$", fontsize="xx-large", transform=ax_left.transAxes)
    ax_left.text(0.02, 0.77, r"$\bar{\rho}_B/\rho_0$", fontsize="xx-large", rotation=90, transform=ax_left.transAxes)
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    ax_left.text(0.008, 0.945, "(a)", fontsize="xx-large", transform=ax_left.transAxes, bbox=bbox)

    ax_right = subfigs[1].subplots(2, 2, sharex=True, sharey=True, gridspec_kw=dict(hspace=0, wspace=0, left=0.013, right=1, bottom=0.038, top=0.995))
    # fins = ["fig/L80_p60_60.jpg", "fig/L80_p80_80.jpg", "fig/L80_p100_100.jpg", "fig/L80_p80_80_ori.jpg"]
    fins = ["fig/p60.jpg", "fig/p80.jpg", "fig/p100.jpg", "fig/p80_ori.jpg"]

    labels = ["(b)", "(c)", "(d)", "(e)"]
    for i, ax in enumerate(ax_right.flat):
        im = mpimg.imread(fins[i])
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.text(0.02, 0.885, labels[i], fontsize="xx-large", transform=ax.transAxes, bbox=bbox)

    ax_cb = ax_right[1, 1].inset_axes([0.75, 0., 0.25, 0.25])
    im = mpimg.imread("fig/circle2.png")
    ax_cb.set_title(r"$\theta_i$", fontsize=20)
    ax_cb.imshow(im)
    ax_cb.set_xticks([])
    ax_cb.set_yticks([])

    ax_left.text(0.75, 0.75, "(b)", fontsize="xx-large", c="k", ha="center", va="center")
    ax_left.text(1, 1, "(c)", fontsize="xx-large", c="k", ha="center", va="center")
    ax_left.text(1.25, 1.25, "(d)", fontsize="xx-large", c="k", ha="center", va="center")

    plt.show()
    # plt.savefig("fig/snaps_pure_NR.pdf", dpi=150)
    plt.close()



if __name__ == "__main__":
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    fig = plt.figure(figsize=(8, 10))
    subfigs = fig.subfigures(2, 1, wspace=0.001, hspace=0.001, height_ratios=[8, 2])

    ax_top = subfigs[0].subplots(1, 1, gridspec_kw=dict(hspace=0, wspace=0, left=0.03, right=0.995, bottom=0.04, top=0.998))
    im = mpimg.imread(f"fig/L20_20_Dr0.1_r80_e0_J1.00_-1.00.jpeg")
    extent = [0.5-0.125, 2+0.125, 0.25-0.125, 1.75+0.125]
    ax_top.imshow(im, extent=extent)
    ax_top.set_xticks([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2])
    ax_top.set_yticks([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75])
    for tick in ax_top.get_yticklabels():
        tick.set_rotation(90)

    Dr_A = Dr_B = 0.1
    etaAA = 0
    etaBB = 0

    etaAB = 1
    etaBA = -1
    bar_rho_A = 1
    bar_rho_B = 1
    bar_vA = 1
    bar_vB = 1
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, bar_vA, bar_vB, qmax=2.5, Nq=400, resolution=2000)
    plot_PD_composition(state, extent, ax=ax_top, fill=False, legend_fs=22)

    ax_top.text(0.8, 0.02, r"$\bar{\rho}_A/\rho_0$", fontsize=26, transform=ax_top.transAxes)
    ax_top.text(0.02, 0.8, r"$\bar{\rho}_B/\rho_0$", fontsize=26, rotation=90, transform=ax_top.transAxes)
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    ax_top.text(0.006, 0.968, "(a)", fontsize=20, transform=ax_top.transAxes, bbox=bbox)

    ax_top.tick_params(axis="x", labelsize="large")
    ax_top.tick_params(axis="y", labelsize="large")

    ax_bottom = subfigs[1].subplots(1, 4, sharex=True, sharey=True, gridspec_kw=dict(hspace=0, wspace=0, left=0., right=1, bottom=0.005, top=1))
    # fins = ["fig/L80_p60_60.jpg", "fig/L80_p80_80.jpg", "fig/L80_p100_100.jpg", "fig/L80_p80_80_ori.jpg"]
    fins = ["fig/p60.jpg", "fig/p80.jpg", "fig/p100.jpg", "fig/p80_ori.jpg"]

    labels = ["(b)", "(c)", "(d)", "(e)"]
    for i, ax in enumerate(ax_bottom):
        im = mpimg.imread(fins[i])
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.text(0.025, 0.875, labels[i], fontsize=20, transform=ax.transAxes, bbox=bbox)

    ax_cb = ax_bottom[-1].inset_axes([0.75, 0., 0.25, 0.25])
    im = mpimg.imread("fig/circle2.png")
    ax_cb.set_title(r"$\theta_i$", fontsize=26)
    ax_cb.imshow(im)
    ax_cb.set_xticks([])
    ax_cb.set_yticks([])

    ax_top.text(0.75, 0.75, "(b)", fontsize=20, c="k", ha="center", va="center")
    ax_top.text(1, 1, "(c)", fontsize=20, c="k", ha="center", va="center")
    ax_top.text(1.25, 1.25, "(d)", fontsize=20, c="k", ha="center", va="center")

    # plt.show()
    plt.savefig("fig/snaps_pure_NR_v2.pdf", dpi=100)
    plt.close()