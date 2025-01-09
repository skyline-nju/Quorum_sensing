import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import matplotlib.patches as mpatches
import sys
from binodal_simulation import plot_PD

sys.path.append("../")
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

from cal_omega_barV import get_bar_v_omega
from general_PD import find_long_instabi, find_short_instabi, find_contours


def get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, bar_vA=1, bar_vB=1, kappa=0.7,
                            extent=[0, 4, 0, 4], qmin=1e-5, qmax=1, Nq=200, resolution=1000, overwrite=False, ll=3./20, simple_gamma=False):
    if simple_gamma:
        folder = "data/PD_pA_pB_SG/"
    elif ll == 0.:
        folder = "data/PD_pA_pB_woS/"
    else:
        folder = "data/PD_pA_pB/"
    if not os.path.exists(folder):
        os.mkdir(folder)
    fout = "%se%g_%g_%g_%g_D%g_%g_r%g_%g_qm%g_r%g.npz" % (
        folder, etaAA, etaAB, etaBA, etaBB, Dr_A, Dr_B, bar_rho_A, bar_rho_B, qmax, resolution
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


def plot_linear_stability_diagram(state, extent, xlim=None, ylim=None, ax=None, mode="fill"):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    nrows, ncols = state.shape
    contours = find_contours(state)
    for contour in contours["LWS"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if mode == "fill":
            ax.fill(x, y, c="tab:blue", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:blue", lw=1, linestyle="dashed")
            pass

    for contour in contours["LWO"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if mode == "fill":
            ax.fill(x, y, c="tab:pink", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:pink", lw=1, linestyle="dashed")
            pass
    
    # for contour in contours["SWS"]:
    #     x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
    #     y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
    #     ax.fill(x, y, c="tab:orange", alpha=0.5)


    for contour in contours["SWO"]:
        x = (contour[:, 1] / ncols) * (extent[1] - extent[0]) + extent[0]
        y = (contour[:, 0] / nrows) * (extent[3] - extent[2]) + extent[2]
        if mode == "fill":
            ax.fill(x, y, c="tab:green", alpha=0.5)
        else:
            ax.plot(x, y, c="tab:green", lw=1, linestyle="dashed")

    if xlim is None:
        xlim = [extent[0], extent[1]]
    if ylim is None:
        ylim = [extent[2], extent[3]]
    if mode == "fill":
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
    if flag_show:
        plt.show()
        plt.close()


def load_composition_plane(ax=None, label_font_size="xx-large"):
    from scipy.interpolate import CubicSpline
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    image = mpimg.imread("fig/L40_40_Dr0.1_r20_e-2_J0.50_reduced.jpeg")
    extent = [-0.125, 3.625, -0.125, 3.625]
    ax.imshow(image, extent=extent)
    ax.text(0.85, 0.02, r"$\bar{\rho}_A/\rho_0$", fontsize=label_font_size, transform=ax.transAxes)
    ax.text(0.02, 0.82, r"$\bar{\rho}_B/\rho_0$", fontsize=label_font_size, rotation=90, transform=ax.transAxes)

    x_u =np.array([0.3, 0.315, 0.34, 0.38, 0.43, 0.5, 0.65, 0.8, 0.9, 1, 1.1])
    y_u =np.array([3.625, 3.375, 3.125, 2.875, 2.625, 2.5, 2.625, 2.875, 3.125, 3.375, 3.625])
    # ax.plot(x_u, y_u, "o")
    
    cs = CubicSpline(x_u, y_u)
    xs = np.linspace(x_u[0], x_u[-1], 100)
    ys = cs(xs)

    ax.plot(xs, ys, "--", lw=4, c="k", alpha=0.5)
    
    x_l = np.array([-0.125, 0.125, 0.375, 0.375, 0.375, 0.375, 0.375, 0.375, 0.125, -0.125])
    y_l = np.array([2.225, 2.2, 2.15, 1.875, 1.655, 1.0, 0.75, 0.45, 0.4, 0.375])

    ax.plot(x_l, y_l, "--", lw=4, c="k", alpha=0.5)
    # cs = CubicSpline(y_l[::-1], x_l[::-1])
    # ys = np.linspace(y_l[-1], y_l[0], 100)
    # xs = cs(ys)

    x_d = np.array([0.3, 0.5, 0.7, 2.1, 2.2, 2.23])
    y_d = np.array([-0.125, 0.125, 0.375, 0.375, 0.125, -0.125])
    ax.plot(x_d, y_d, "--", lw=4, c="k", alpha=0.5)

    # ax.plot([0, 3], [0, 3], c="w", linestyle="dotted", lw=4)

    x_r = np.array([3.625, 3.375, 3.125, 2.875, 2.655, 2.55, 2.655, 2.875, 3.125, 3.375, 3.625])
    y_r = np.array([0.6, 0.75, 0.95, 1.2, 1.5, 1.675, 1.8, 1.875, 1.925, 1.95, 1.96])
    # ax.plot(x_r, y_r, "o", c="tab:red")
    ax.plot(x_r, y_r, "--", lw=4, c="k", alpha=0.5)
    # ax.set_yticklabels(ax.get_yticks(), rotation=90)
    for tick in ax.get_yticklabels():
        tick.set_rotation(90)
    # ax.plot(0.3, 0.6, "o", c="tab:blue", ms=8)
    # ax.plot(0.55, 0.385, "o", c="tab:pink", ms=8)
    # ax.plot(2, 2, "o", c="tab:green", ms=8)

    ax.set_ylim(ymax=3.375)
    # pA1 = [7.6525/10, 8.625/10, 9.597/10]
    # pB1 = [10.5675/10, 13.36/10, 16.152/10]
    # c = ["tab:cyan", "tab:red", "tab:orange"]
    # for j in range(3):
    #     ax.plot(pA1[j], pB1[j], "s", ms=8, c=c[j])

    # ax.plot(1, 8.7 * 2/10, "o", c="tab:brown", ms=8)

    # gas_binodals = [0.3886, 0.111]
    # liquid_binodals = [1.1689, 2.353]
    # rhoA_b = [gas_binodals[0], liquid_binodals[0]]
    # rhoB_b = [gas_binodals[1], liquid_binodals[1]]
    # ax.plot(gas_binodals[0], gas_binodals[1], "v", c="k")
    # ax.plot(liquid_binodals[0], liquid_binodals[1], "^", c="k")
    # ax.plot(rhoA_b, rhoB_b, "k:")

    # k = (rhoB_b[1] - rhoB_b[0]) / (rhoA_b[1] - rhoA_b[0])
    # x = 0.47
    # y = (x-0.55) * k + 0.385
    # ax.plot(x, y, "o", fillstyle="none")
    # print(x, y)

    # x = 0.5
    # y = (x-0.55) * k + 0.385
    # ax.plot(x, y, "o", fillstyle="none")
    # print(x, y)

    # x = 0.48
    # y = (x-0.55) * k + 0.385
    # ax.plot(x, y, "o", fillstyle="none")
    # print(x, y)

    # x = 0.475
    # y = (x-0.55) * k + 0.385
    # ax.plot(x, y, "o", fillstyle="none")
    # print(x, y)

    if flag_show:
        plt.show()
        plt.close()


def plot_PD_J05():
    from binodal_simulation import plot_PD

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7.5), width_ratios=[10, 7], constrained_layout=True)
    load_composition_plane(ax1)
    ax1_in = ax1.inset_axes([0.69, 0.69, 0.31, 0.31])
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    plot_linear_stability_diagram(state, extent, xlim=[0, 3.5], ylim=[0, 3.25], ax=ax1_in)

    patches = [mpatches.Patch(color='tab:blue', label='LSI', alpha=0.5),
            mpatches.Patch(color='tab:pink', label='LOI',alpha=0.5),
            mpatches.Patch(color='tab:green', label='SOI', alpha=0.5),
            ]
    ax1_in.legend(handles=patches, loc="upper right", fontsize="x-large", borderpad=0.2, labelspacing=0.2, handlelength=1.5)
    ax1_in.text(0.75, 0.06, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax1_in.transAxes)
    ax1_in.text(0.05, 0.74, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax1_in.transAxes)
    # ax1_in.plot([0, 2], [0, 1.92], c="tab:grey", linestyle="dotted", lw=2)
    # ax1_in.set_yticks([0, 1, 2, 3])
    plot_PD(ax2)
    # ax2.plot(3/4, 5/4, "o")

    Lx_in = 0.4 * 1.1
    Ly_in = 0.3 * 1.1
    ax2_in = ax2.inset_axes([1-Lx_in, 1-Ly_in, Lx_in, Ly_in])
    plot_linear_stability_diagram(state, extent, ax=ax2_in, mode="line")
    plot_PD(ax2_in, show_tie_line=False)
    ax2_in.set_xlim(0, 3.5)
    ax2_in.set_ylim(0, 3.4)
    ax2_in.text(0.75, 0.03, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax2_in.transAxes)
    ax2_in.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax2_in.transAxes)
    fig.text(0.02, 0.965, "(a)",fontsize="xx-large")
    fig.text(0.6, 0.965, "(b)",fontsize="xx-large")
    # ax2.plot(1, 2.8, "o")
    # ax2.plot(0.45, 1.47, "o")
    # ax2.plot(0.4, 1.0, "o")
    # ax2.plot(0.35, 0.8, "o")
    # ax2.plot(0.43, 0.79, "o")


    # binodals_Jp05(ax2)
    plt.show()
    # plt.savefig("fig/PD_J05.pdf")
    plt.close()


def PD_small_large():
    from binodal_simulation import plot_PD

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7.5), width_ratios=[10, 7], constrained_layout=True)
    load_composition_plane(ax1, label_font_size=20)
    ax1_in = ax1.inset_axes([0.69, 0.69, 0.31, 0.31])
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    plot_linear_stability_diagram(state, extent, xlim=[0, 3.5], ylim=[0, 3.25], ax=ax1_in)

    patches = [mpatches.Patch(color='tab:blue', label='LSI', alpha=0.5),
            mpatches.Patch(color='tab:pink', label='LOI',alpha=0.5),
            mpatches.Patch(color='tab:green', label='SOI', alpha=0.5),
            ]
    ax1_in.legend(handles=patches, loc="upper right", fontsize="x-large", borderpad=0.2, labelspacing=0.2, handlelength=1.5)
    ax1_in.text(0.75, 0.06, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax1_in.transAxes)
    ax1_in.text(0.05, 0.74, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax1_in.transAxes)
    # ax1_in.plot([0, 2], [0, 1.92], c="tab:grey", linestyle="dotted", lw=2)
    # ax1_in.set_yticks([0, 1, 2, 3])
    plot_PD(ax2, label_font_size=20)
    # ax2.plot(3/4, 5/4, "o")

    Lx_in = 0.4 * 1.1
    Ly_in = 0.3 * 1.1
    ax2_in = ax2.inset_axes([1-Lx_in, 1-Ly_in, Lx_in, Ly_in])
    plot_linear_stability_diagram(state, extent, ax=ax2_in, mode="line")
    plot_PD(ax2_in, show_tie_line=False)
    ax2_in.set_xlim(0, 3.5)
    ax2_in.set_ylim(0, 3.4)
    ax2_in.text(0.75, 0.03, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax2_in.transAxes)
    ax2_in.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax2_in.transAxes)
    fig.text(0.02, 0.965, "(a)",fontsize=20)
    fig.text(0.6, 0.965, "(b)",fontsize=20)
    # ax2.plot(1, 2.8, "o")
    # ax2.plot(0.45, 1.47, "o")
    # ax2.plot(0.4, 1.0, "o")
    # ax2.plot(0.35, 0.8, "o")
    # ax2.plot(0.43, 0.79, "o")


    # binodals_Jp05(ax2)
    plt.show()
    # plt.savefig("fig/PD_J05_2.pdf")
    plt.close()


def PD_small():
    fig, ax = plt.subplots(1, 1, figsize=(8, 7.5), constrained_layout=True)
    load_composition_plane(ax, label_font_size=20)
    ax_in = ax.inset_axes([0.69, 0.69, 0.31, 0.31])
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    plot_linear_stability_diagram(state, extent, xlim=[0, 3.5], ylim=[0, 3.25], ax=ax_in)

    patches = [mpatches.Patch(color='tab:blue', label='LSI', alpha=0.5),
            mpatches.Patch(color='tab:pink', label='LOI',alpha=0.5),
            mpatches.Patch(color='tab:green', label='SOI', alpha=0.5),
            ]
    ax_in.legend(handles=patches, loc="upper right", fontsize="x-large", borderpad=0.2, labelspacing=0.2, handlelength=1.5)
    ax_in.text(0.75, 0.06, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax_in.transAxes)
    ax_in.text(0.05, 0.74, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax_in.transAxes)
    plt.show()
    # plt.savefig("fig/PD_L40_J05.pdf")
    plt.close()


def PD_large():
    fig = plt.figure(figsize=(6, 7.5), layout="constrained")
    subfig = fig.subfigures(1, 1, wspace=0.0001, hspace=0.0001)

    ax1 = subfig.subplots(1, 1)
    plot_PD(ax1, label_font_size=20)
    Lx_in = 0.4 * 1.1
    Ly_in = 0.3 * 1.1
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    ax1_in = ax1.inset_axes([1-Lx_in, 1-Ly_in, Lx_in, Ly_in])
    plot_linear_stability_diagram(state, extent, ax=ax1_in, mode="line")
    plot_PD(ax1_in, show_tie_line=False)
    ax1_in.set_xlim(0, 3.5)
    ax1_in.set_ylim(0, 3.4)
    ax1_in.text(0.75, 0.03, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax1_in.transAxes)
    ax1_in.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax1_in.transAxes)


    plt.show()
    # plt.savefig("fig/PD_large_J05.pdf")
    plt.close()


def PD_and_snaps():
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7.5), width_ratios=[7, 10], constrained_layout=True)
    fig = plt.figure(figsize=(9, 7.5), layout="constrained")
    subfigs = fig.subfigures(1, 2, wspace=0.0001, hspace=0.0001, width_ratios=[7, 4.6])

    ax1 = subfigs[0].subplots(1, 1)
    plot_PD(ax1, label_font_size=20)
    Lx_in = 0.4 * 1.1
    Ly_in = 0.3 * 1.1
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    ax1_in = ax1.inset_axes([1-Lx_in, 1-Ly_in, Lx_in, Ly_in])
    plot_linear_stability_diagram(state, extent, ax=ax1_in, mode="line")
    plot_PD(ax1_in, show_tie_line=False)
    ax1_in.set_xlim(0, 3.5)
    ax1_in.set_ylim(0, 3.4)
    ax1_in.text(0.75, 0.03, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax1_in.transAxes)
    ax1_in.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax1_in.transAxes)

    ax_snaps = subfigs[1].subplots(6, 1, sharex=True, sharey=True, gridspec_kw=dict(hspace=0, wspace=0, left=0, right=1, bottom=0., top=1))

    fnames = ["G_CCB.png", "G_LA_CCB.png", "LA_CCB.png",
              "G_LB_CCB.png", "LB_CCB_0.png", "LB_CCB_1.png",
              "G_LB_CCB_L640.png", "LB_CCB_L640.png"]
    titles = ["(a) G+CCB", "(b) G+LA+CCB", "(c) LA+CCB", "(d) G+LB+CCB",
              "(e) LB+CCB: initial configuration", "(f) LB+CCB: stready state",
              "(g) G+LB+CCB", "(h) LB+CCB"]
    for i, ax in enumerate(ax_snaps):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        # ax.set_title(titles[i], fontsize="x-large")

    # plt.show()
    plt.savefig("fig/PD_snaps_J05.pdf")
    plt.close()


def PD_L40_instability():
    fig = plt.figure(figsize=(16, 7))

    subfigs = fig.subfigures(1, 2, wspace=0.01, width_ratios=[1, 1.25])

    gridspec_kw=dict(hspace=0, wspace=0, left=0.03, right=0.98, bottom=0.01, top=0.98)
    ax = subfigs[0].subplots(1, 1, gridspec_kw=gridspec_kw)

    label_fs = "x-large"
    ax.set_title(r"(a) $L_x=L_y=40$", fontsize=label_fs)

    subfigsnest = subfigs[1].subfigures(2, 1, height_ratios=[1, 1], hspace=0.05)

    gridspec_kw=dict(hspace=0, wspace=0, left=0.01, right=0.999, bottom=0.0, top=0.95)
    ax_SB = subfigsnest[0].subplots(4, 3, width_ratios=[1, 2, 4], gridspec_kw=gridspec_kw)

    gridspec_kw=dict(hspace=0.05, wspace=0, left=0.01, right=0.999, bottom=-0.025, top=1)
    ax_SP = subfigsnest[1].subplots(2, 6, gridspec_kw=gridspec_kw)
    # fig, ax = plt.subplots(1, 1, figsize=(8, 7.5), constrained_layout=True)
    load_composition_plane(ax, label_font_size=20)
    ax_in = ax.inset_axes([0.69, 0.69, 0.31, 0.31])
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    plot_linear_stability_diagram(state, extent, xlim=[0, 3.5], ylim=[0, 3.25], ax=ax_in)

    mk1, = ax.plot(0.75, 1.25, "^", c="tab:cyan", ms=8)
    mk2, = ax.plot(0.575, 0.275, "s", c="tab:olive", ms=7)
    mk3, = ax.plot(0.75, 0.5, "o", c="tab:orange", ms=7)


    patches = [mpatches.Patch(color='tab:blue', label='LSI', alpha=0.5),
            mpatches.Patch(color='tab:pink', label='LOI',alpha=0.5),
            mpatches.Patch(color='tab:green', label='SOI', alpha=0.5),
            ]
    ax_in.legend(handles=patches, loc="upper right", fontsize="x-large", borderpad=0.2, labelspacing=0.2, handlelength=1.5)
    ax_in.text(0.75, 0.06, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax_in.transAxes)
    ax_in.text(0.05, 0.74, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax_in.transAxes)
    
    
    Lx_arr = [80, 160, 320]
    t_labels = [r"$t=0$", r"$t=4e3$", r"$t=5e3$", r"$t=1e4$"]
    labels = ["(b)", "(c)", "(d)"]
    for row in range(4):
        for col in range(3):
            fname = "fig/SB_instab/%d_%d.png" % (Lx_arr[col], row)
            im = mpimg.imread(fname)
            ax = ax_SB[row, col]
            ax.imshow(im)
            ax.set_xticks([])
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(t_labels[row], fontsize="x-large")
            if row == 0:
                ax.set_title(r"%s $L_x=%d$" % (labels[col], Lx_arr[col]), fontsize="x-large")
    dx = 0.05
    ax_in = ax_SB[0, -1].inset_axes([1-dx, 1-dx*8, dx, dx*8])
    ax_in.set_xticklabels([])
    ax_in.set_yticklabels([])
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.plot(0, 0, "^", c=mk1.get_c(), ms=8)

    t_arr = [160, 200, 300, 400, 1000, 20000]
    titles = [r"(e) $t=160$", r"$t=200$", r"$t=300$", r"$t=400$", r"$t=10^3$", r"$t=2\times 10^4$"]
    # titles = [r"(e) $t=160$", r"(f) $t=200$", r"(g) $t=300$", r"(h) $t=400$", r"(i) $t=10^3$", r"(j) $t=2\times 10^4$"]

    for i, ax in enumerate(ax_SP[0]):
        fname = "fig/spiral_instab/%d.png" % t_arr[i]
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[i], fontsize="x-large")
    ax_SP[0,0].set_ylabel(r"$L_x=L_y=160$", fontsize="x-large")

    dx = 0.17
    ax_in = ax_SP[0, -1].inset_axes([1-dx, 1-dx, dx, dx])
    ax_in.set_xticklabels([])
    ax_in.set_yticklabels([])
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.plot(0, 0, "s", c=mk2.get_c(), ms=8)

    t_arr = [0, 500, 1000, 2020, 3000, 5000]
    titles = [r"(f) $t=0$", r"$t=500$", r"$t=1000$", r"$t=2020$", r"$t=3000$", r"$t=5000$"]
    # titles = [r"(k) $t=0$", r"(l) $t=500$", r"(m) $t=1000$", r"(n) $t=2020$", r"(o) $t=3000$", r"(p) $t=5000$"]
    for i, ax in enumerate(ax_SP[1]):
        fname = "fig/spiral_instab2/%d.png" % t_arr[i]
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[i], fontsize="x-large")
    ax_SP[1,0].set_ylabel(r"$L_x=L_y=80$", fontsize="x-large")

    dx = 0.17
    ax_in = ax_SP[1, -1].inset_axes([1-dx, 1-dx, dx, dx])
    ax_in.set_xticklabels([])
    ax_in.set_yticklabels([])
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.plot(0, 0, "o", c=mk3.get_c(), ms=8)

    plt.show()
    # plt.savefig("fig/PD_L40_J05_instability.pdf", dpi=150)
    plt.close()


def PD_large_snap_profile():
    from matplotlib.patches import ConnectionPatch
    fs ="x-large"
    fig = plt.figure(figsize=(13, 8), layout="constrained")
    subfigs = fig.subfigures(1, 2, wspace=0.0001, hspace=0.0001, width_ratios=[1, 1.3])
    ax1 = subfigs[0].subplots(1, 1)
    plot_PD(ax1, label_font_size=20, ms=5)
    Lx_in = 0.4 * 1.1
    Ly_in = 0.3 * 1.1
    Dr = 0.1
    etaAA = etaBB = -2
    etaAB = 0.5
    etaBA = -etaAB
    extent, state, q_range = get_PD_composition_data(etaAA, etaAB, etaBA, etaBB, Dr, Dr, 1, 1, qmax=2.5, Nq=400, resolution=2000, extent=[0, 4.5, 0, 6.5])
    ax1_in = ax1.inset_axes([1-Lx_in, 1-Ly_in, Lx_in, Ly_in])
    plot_linear_stability_diagram(state, extent, ax=ax1_in, mode="line")
    plot_PD(ax1_in, show_tie_line=False)
    ax1_in.set_xlim(0, 3.5)
    ax1_in.set_ylim(0, 3.4)
    ax1_in.text(0.75, 0.03, r"$\bar{\rho}_A/\rho_0$", fontsize="x-large", transform=ax1_in.transAxes)
    ax1_in.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize="x-large", rotation=90, transform=ax1_in.transAxes)
    ax1.set_title(r"(a) $\eta^0_{AB}=-\eta^0_{BA}=0.5$, $\eta^0_{AA}=\eta^0_{BB}=-2$, $D_r=0.1,\rho_0=10$", fontsize=fs)
    ax_right = subfigs[1].subplots(5, 2)
    ms = 8
    clist = []
    mks = ["o", "P", "p", "s", "*"]
    line, = ax1.plot(0.75, 1, "o", ms=ms)
    clist.append(line.get_c())
    line, = ax1.plot(1.5, 1, "P", ms=ms)
    clist.append(line.get_c())
    line, = ax1.plot(1.8, 1.45, "p", ms=ms)
    clist.append(line.get_c())
    line, = ax1.plot(0.45, 1.47, "s", ms=ms)
    clist.append(line.get_c())
    line, = ax1.plot(1, 2.8, "*", ms=ms)
    clist.append(line.get_c())

    ax1.plot(1, 2.5, "v", ms=ms, fillstyle="none", c="tab:purple")
    ax1.plot(1, 3.4, ">", ms=ms, fillstyle="none", c="tab:purple")
    ax1.plot(1, 3.8, "^", ms=ms, fillstyle="none", c="b")
    ax1.plot(0.825, 3.4, "<", ms=ms, fillstyle="none", c="tab:purple")


    fnames = ["G_CCB.png", "G_LA_CCB.png",
              "LA_CCB.png", "none",
              "G_LB_CCB.png", "none", 
              "LB_CCB_L120_40.png", "none",
              "LB_CCB_0.png", "LB_CCB_1.png"]
    titles = ["(b) G+CCB", "(c) G+LA+CCB",
              "(d) LA+CCB", "none",
              "(f) G+LB+CCB", "none",
              r"(h) LB+CCB: $L_x=3L_y=120$", "none",
              "(j) LB+CCB: initial configuration", "(k) LB+CCB: stready state"]
    count = 0
    for i, ax in enumerate(ax_right.flat):
        if fnames[i] != "none":
            fname = f"fig/snap/{fnames[i]}"
            im = mpimg.imread(fname)
            ax.imshow(im)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            ax.set_title(titles[i], fontsize=fs)

            dx = 0.05
            ax_in = ax.inset_axes([1-dx, 1-dx*3, dx, dx*3])
            ax_in.set_xticklabels([])
            ax_in.set_yticklabels([])
            ax_in.set_xticks([])
            ax_in.set_yticks([])
            ax_in.plot(0, 0, mks[count], c=clist[count], ms=8)
            if count < 4:
                count +=1
    
    ax1, ax2, ax3 = ax_right[1:4, 1]
    with np.load("data/time_ave_profile/LA_CCB_profiles_p18_14.5_L480_160.npz", "rb") as data:
        x, fx = data["x"], data["fx"]
        # x /= 480
        ax1.plot(x, fx[0] / 10, label=r"$S=A$", c="tab:blue")
        ax1.plot(x, fx[1] / 10, label=r"$S=B$", c="tab:red")
        ax1.set_xlim(0, 480)
        # ax1.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax1.set_title(r"(e) $\langle\rho_S\rangle_{y,t}/\rho_0 $ for LA+CCB", fontsize=fs)
        ax1.text(0.95, 0.08, r"$x$", fontsize=fs, transform=ax1.transAxes)
    
    with np.load("data/instant_profile/L480_160_Dr0.100_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_1001.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax2.plot(x, fields[0, 0]/5, c="tab:blue")
        ax2.plot(x, fields[0, 1]/5, c="tab:red")
        # ax2.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax2.set_title(r"(g) $\langle\rho_S\rangle_{y}/\rho_0 $ for G+LB+CCB", fontsize=fs)
        ax2.set_xlim(0, 480)
        ax2.text(0.95, 0.01, r"$x$", fontsize=fs, transform=ax2.transAxes)


    with np.load("data/instant_profile/L120_40_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.100_1000.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax3.plot(x, fields[-1, 0]/10, c="tab:blue")
        ax3.plot(x, fields[-1, 1]/10, c="tab:red")
        # ax3.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax3.set_title(r"(i) $\langle\rho_S\rangle_{y}/\rho_0 $ for LB+CCB", fontsize=fs)
        ax3.set_xlim(0, 120)
        ax3.set_xticks([0, 50, 100])
        ax3.text(0.95, 0.09, r"$x$", fontsize=fs, transform=ax3.transAxes)



    ylim = ax_right[4, 0].get_ylim()
    xlim = ax_right[4, 0].get_xlim()
    y = (ylim[1] - ylim[0]) * 0.75 + ylim[0]
    x1 = (xlim[1]-xlim[0]) / 4 * 2 + xlim[0]
    x2 = (xlim[1]-xlim[0]) / 4 * 3 + xlim[0]
    ax_right[4, 0].axhline(y, 0.5, 0.75, linestyle="dashed", c="tab:grey")
    ax_right[4, 0].axvline(x1, 0.75, 1, linestyle="dashed", c="tab:grey")
    ax_right[4, 0].axvline(x2, 0.75, 1, linestyle="dashed", c="tab:grey")

    con = ConnectionPatch(xyA=(0, 0), coordsA=ax_right[3, 0].transAxes, xyB=(0.5, 1), coordsB=ax_right[4, 0].transAxes,
                          linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)
    subfigs[1].add_artist(con)

    con = ConnectionPatch(xyA=(1, 0), coordsA=ax_right[3, 0].transAxes, xyB=(0.75, 1), coordsB=ax_right[4, 0].transAxes,
                          linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)
    subfigs[1].add_artist(con)

    ax_right[2, 0].arrow(0.75, 0.5, 0.12, 0, transform=ax_right[2,0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_right[2, 1].arrow(0.75, 0.5, 0.12, 0, transform=ax_right[2,1].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_right[3, 0].arrow(0.3, 0.5, 0.06, 0, transform=ax_right[3,0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_right[3, 0].arrow(0.7, 0.5, 0.06, 0, transform=ax_right[3,0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_right[3, 1].arrow(0.25, 0.75, 0.06, 0, transform=ax_right[3,1].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_right[3, 1].arrow(0.65, 0.75, 0.06, 0, transform=ax_right[3,1].transAxes, width=0.04, color="k", ec="k", head_length=0.04)

    ax_right[4, 0].arrow(0.075, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.2, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.325, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.45, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.575, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.7, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.825, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_right[4, 0].arrow(0.935, 0.5, 0.04, 0, transform=ax_right[4, 0].transAxes, width=0.02, color="k", ec="k", head_length=0.02)

    ax_right[1,1].legend(fontsize="large", loc=(0.3, 0.01), frameon=False)

    x = np.array([0.9, 0.64, 0.4, 0.8])
    y = np.array([0.5, 0.48, 0.5, 0.7])
    theta = np.array([-30, -40, 20, 55], float) * np.pi / 180
    for i in range(x.size):
        ax_right[4, 1].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]),
                             transform=ax_right[4, 1].transAxes, width=0.015, color="k", ec="k", head_length=0.03)

    plt.show()
    # plt.savefig("fig/PD_large_J05_snap_profile.pdf", dpi=150)
    plt.close()

if __name__ == "__main__":
    # PD_and_snaps()
    # PD_large_snap_profile()

    PD_L40_instability()