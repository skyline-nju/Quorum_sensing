import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib import image as mpimg

import os
import sys
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

def varied_Ly(axes=None):

    fig, axes = plt.subplots(2, 3, figsize=(4, 6), sharex=True, sharey=True, constrained_layout=True)

    Ly_arr = [40, 80, 160]
    Lx = 240

    for col, Ly in enumerate(Ly_arr):
        fname = "data/space_time/L240_%d_Dr0.100_k0.70_p8_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_2000.npz" % Ly
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = 1e5
                rhoA = fields0[80:130, 0] / 10
                rhoB = fields0[80:130, 1] /10
            else:
                T = (t0[-1] - t0[0]) * 0.1
                rhoA = fields0[1:, 0] / 10
                rhoB = fields0[1:, 1] / 10
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            axes[0, col].imshow(rhoA, origin="lower", extent=extent, aspect="auto", vmin=3.5/10, vmax=22/10)
            axes[1, col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=24/10, vmax=60/10)
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0, 0].ticklabel_format(axis="y", scilimits=[-5, 4])
    plt.show()
    plt.close()


def varied_rB():
    fig, axes = plt.subplots(2, 3, figsize=(4, 6), sharex=True, sharey="col", constrained_layout=True)

    pB_arr = [76, 68, 50]
    Lx = 120

    for col, pB in enumerate(pB_arr):
        fname = "data/space_time/L120_40_Dr0.100_k0.70_p20_%d_r20_20_20_e-2.000_J0.500_-0.500_h0.100_2000.npz" % pB
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = (t0[-1] - t0[150]) * 0.1
                rhoA = fields0[150:, 0] / 20
                rhoB = fields0[150:, 1] / 20
            elif col == 1:
                T = (t0[-1] - t0[9]) * 0.1
                rhoA = fields0[9:, 0] / 20
                rhoB = fields0[9:, 1] / 20
            else:
                # T = (t0[99] - t0[49]) * 0.1
                # rhoA = fields0[49:100, 0] / 20
                # rhoB = fields0[49:100, 1] / 20
                T = (t0[-1] - t0[99]) * 0.1
                rhoA = fields0[99:, 0] / 20
                rhoB = fields0[99:, 1] / 20
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            axes[0, col].imshow(rhoA, origin="lower", extent=extent, aspect="auto", vmin=3.5/20, vmax=22/20)
            axes[1, col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=1.5, vmax=4.5)
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0, 0].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[0, 1].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[0, 2].ticklabel_format(axis="y", scilimits=[-5, 3])

    plt.show()
    plt.close()


def plot_osi_old():
    fig = plt.figure(figsize=(12, 8), layout='constrained')

    subfigs = fig.subfigures(1, 2, wspace=0.01, hspace=0.001, width_ratios=[7.5, 2.5])

    axes_left = subfigs[0].subplots(2, 3)

    pB_arr = [76, 68, 50]
    Lx = 120
    axes = axes_left[0]
    titles = [r"(a) $\bar{\rho}_B/\rho_0=3.8$", r"(b) $\bar{\rho}_B/\rho_0=3.4$", r"(c) $\bar{\rho}_B/\rho_0=2.5$"]
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    for col, pB in enumerate(pB_arr):
        fname = "data/space_time/L120_40_Dr0.100_k0.70_p20_%d_r20_20_20_e-2.000_J0.500_-0.500_h0.100_2000.npz" % pB
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = (t0[-1] - t0[0]) * 0.1
                rhoA = fields0[0:, 0] / 20
                rhoB = fields0[0:, 1] / 20
            elif col == 1:
                T = (t0[-1] - t0[9]) * 0.1
                rhoA = fields0[9:, 0] / 20
                rhoB = fields0[9:, 1] / 20
            else:
                # T = (t0[99] - t0[49]) * 0.1
                # rhoA = fields0[49:100, 0] / 20
                # rhoB = fields0[49:100, 1] / 20
                T = (t0[-1] - t0[199]) * 0.1
                rhoA = fields0[199:, 0] / 20
                rhoB = fields0[199:, 1] / 20
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            im = axes[col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=2, vmax=5)
            axes[col].set_title(titles[col], fontsize="x-large")
    # axes[0].text(0.02, 0.935, "(a)", fontsize="large", transform=axes[0].transAxes, bbox=bbox)
    # axes[1].text(0.02, 0.935, "(b)", fontsize="large", transform=axes[1].transAxes, bbox=bbox)
    # axes[2].text(0.02, 0.935, "(c)", fontsize="large", transform=axes[2].transAxes, bbox=bbox)

    axes[1].axhline(54300-9000, color="w", linestyle="dashed")
    axes[1].axhline(56400-9000, color="w", linestyle="dashed")


    
    cb = subfigs[0].colorbar(im, shrink=0.5, ax=axes_left, location="right",  extend="both")
    axes[0].set_ylabel(r"$t$", fontsize="x-large")
    # cb.set_label(r"$\langle \rho_B(\mathbf{r},t)\rangle_y /\rho_0$", fontsize="x-large")
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[1].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[2].ticklabel_format(axis="y", scilimits=[-5, 3])

    axes = axes_left[1]
    Ly_arr = [40, 80, 160]
    Lx = 240
    titles = [r"(d) $L_y=40$", r"(e) $L_y=80$", r"(f) $L_y=160$"]

    for col, Ly in enumerate(Ly_arr):
        fname = "data/space_time/L240_%d_Dr0.100_k0.70_p8_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_2000.npz" % Ly
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = 1e5
                rhoA = fields0[80:130, 0] / 10
                rhoB = fields0[80:130, 1] /10
            else:
                T = (t0[-1] - t0[0]) * 0.1
                rhoA = fields0[1:, 0] / 10
                rhoB = fields0[1:, 1] / 10
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            axes[col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=2, vmax=5)
            axes[col].set_title(titles[col], fontsize="x-large")
            axes[col].set_xlabel(r"$x$", fontsize="x-large")
    axes[0].set_ylabel(r"$t$", fontsize="x-large")
    # axes[0].text(0.02, 0.935, "(d)", fontsize="large", transform=axes[0].transAxes, bbox=bbox)
    # axes[1].text(0.02, 0.935, "(e)", fontsize="large", transform=axes[1].transAxes, bbox=bbox)
    # axes[2].text(0.02, 0.935, "(f)", fontsize="large", transform=axes[2].transAxes, bbox=bbox)
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[1].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[2].ticklabel_format(axis="y", scilimits=[-5, 4])
    subfigs[0].text(0.913, 0.765, r"$\frac{\langle \rho_B(\mathbf{r},t)\rangle_y}{\rho_0}$", fontsize="xx-large")


    ax = subfigs[1].subplots(1, 1)
    im = mpimg.imread("fig/osi_snaps.jpg")
    ax.imshow(im)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    # ax.axis("off")
    ax.set_title("(g)", fontsize="x-large")
    # ax.text(0.013, 0.973, "(g)", fontsize="large", transform=ax.transAxes, bbox=bbox)

    # con = ConnectionPatch(xyA=(T, 54300-9000), coordsA=axes_left[0, 1].transData, xyB=, coordsB=ax.transAxes)
    con = ConnectionPatch(xyA=(0.025, 1), coordsA=ax.transAxes, xyB=(120, 56400-9000), coordsB=axes_left[0, 1].transData, linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)
    fig.add_artist(con)


    con = ConnectionPatch(xyA=(0.025, 0), coordsA=ax.transAxes, xyB=(120, 54300-9000), coordsB=axes_left[0, 1].transData, linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)

    fig.add_artist(con)

    plt.show()
    # plt.savefig("fig/osi_to_travelling.pdf")
    plt.close()


def plot_space_time_snaps():
    fig = plt.figure(figsize=(8, 6.5), layout='constrained')
    subfigs = fig.subfigures(1, 2, wspace=0.01, hspace=0.001, width_ratios=[7.5, 2.5])

    axes_left = subfigs[0].subplots(2, 3)

    pB_arr = [76, 68, 50]
    Lx = 120
    axes = axes_left[0]
    titles = [r"(a) $\bar{\rho}_B/\rho_0=3.8$", r"(b) $\bar{\rho}_B/\rho_0=3.4$", r"(c) $\bar{\rho}_B/\rho_0=2.5$"]
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    for col, pB in enumerate(pB_arr):
        fname = "data/space_time/L120_40_Dr0.100_k0.70_p20_%d_r20_20_20_e-2.000_J0.500_-0.500_h0.100_2000.npz" % pB
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = (t0[-1] - t0[0]) * 0.1
                rhoA = fields0[0:, 0] / 20
                rhoB = fields0[0:, 1] / 20
            elif col == 1:
                T = (t0[-1] - t0[9]) * 0.1
                rhoA = fields0[9:, 0] / 20
                rhoB = fields0[9:, 1] / 20
            else:
                # T = (t0[99] - t0[49]) * 0.1
                # rhoA = fields0[49:100, 0] / 20
                # rhoB = fields0[49:100, 1] / 20
                T = (t0[-1] - t0[199]) * 0.1
                rhoA = fields0[199:, 0] / 20
                rhoB = fields0[199:, 1] / 20
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            im = axes[col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=2, vmax=5)
            axes[col].set_title(titles[col], fontsize="xx-large")
    # axes[0].text(0.02, 0.935, "(a)", fontsize="large", transform=axes[0].transAxes, bbox=bbox)
    # axes[1].text(0.02, 0.935, "(b)", fontsize="large", transform=axes[1].transAxes, bbox=bbox)
    # axes[2].text(0.02, 0.935, "(c)", fontsize="large", transform=axes[2].transAxes, bbox=bbox)

    axes[1].axhline(54300-9000, color="w", linestyle="dashed")
    axes[1].axhline(56400-9000, color="w", linestyle="dashed")
    
    axes[0].set_ylabel(r"$t$", fontsize="xx-large")
    axes[0].set_yticks([0, 0.5e5, 1e5, 1.5e5, 2e5])
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[1].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[2].ticklabel_format(axis="y", scilimits=[-5, 3])

    axes = axes_left[1]
    Ly_arr = [40, 80, 160]
    Lx = 240
    titles = [r"(e) $L_y=40$", r"(f) $L_y=80$", r"(g) $L_y=160$"]

    for col, Ly in enumerate(Ly_arr):
        fname = "data/space_time/L240_%d_Dr0.100_k0.70_p8_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_2000.npz" % Ly
        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]
            print(t0.size)

            if col == 0:
                T = 1e5
                rhoA = fields0[80:130, 0] / 10
                rhoB = fields0[80:130, 1] /10
            else:
                T = (t0[-1] - t0[0]) * 0.1
                rhoA = fields0[1:, 0] / 10
                rhoB = fields0[1:, 1] / 10
            print(rhoA.min(), rhoA.max(), rhoB.min(), rhoB.max())
            extent = [0, Lx, 0, T]
            axes[col].imshow(rhoB, origin="lower", extent=extent, aspect="auto", vmin=2, vmax=5)
            axes[col].set_title(titles[col], fontsize="xx-large")
            axes[col].set_xlabel(r"$x$", fontsize="xx-large")
    axes[0].set_ylabel(r"$t$", fontsize="xx-large")
    # axes[0].text(0.02, 0.935, "(d)", fontsize="large", transform=axes[0].transAxes, bbox=bbox)
    # axes[1].text(0.02, 0.935, "(e)", fontsize="large", transform=axes[1].transAxes, bbox=bbox)
    # axes[2].text(0.02, 0.935, "(f)", fontsize="large", transform=axes[2].transAxes, bbox=bbox)
    
    # axes[0, 0].set_ylim(200000, 300000)
    axes[0].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[1].ticklabel_format(axis="y", scilimits=[-5, 4])
    axes[2].ticklabel_format(axis="y", scilimits=[-5, 4])
    # subfigs[0].text(0.913, 0.765, r"$\frac{\langle \rho_B(\mathbf{r},t)\rangle_y}{\rho_0}$", fontsize="xx-large")


    (ax, ax_cb) = subfigs[1].subplots(2, 1, height_ratios=[9.5, 0.5])
    cb = subfigs[0].colorbar(im, shrink=0.5, ax=axes_left, cax=ax_cb, location="bottom",  extend="both")
    cb.set_label(r"$\langle \rho_B(\mathbf{r},t)\rangle_y /\rho_0$", fontsize="xx-large")

    im = mpimg.imread("fig/osi_snaps.jpg")
    ax.imshow(im)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    # ax.axis("off")
    ax.set_title("(d)", fontsize="xx-large")
    # ax.text(0.013, 0.973, "(g)", fontsize="large", transform=ax.transAxes, bbox=bbox)

    # con = ConnectionPatch(xyA=(T, 54300-9000), coordsA=axes_left[0, 1].transData, xyB=, coordsB=ax.transAxes)
    con = ConnectionPatch(xyA=(0.025, 1), coordsA=ax.transAxes, xyB=(120, 56400-9000), coordsB=axes_left[0, 1].transData, linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)
    fig.add_artist(con)


    con = ConnectionPatch(xyA=(0.025, 0), coordsA=ax.transAxes, xyB=(120, 54300-9000), coordsB=axes_left[0, 1].transData, linestyle=":", color="tab:grey", lw=1.5)
    con.set_annotation_clip(False)

    fig.add_artist(con)

    mk = ["^", ">", "v", "<", "<", "<"]
    color = ["b", "tab:purple", "tab:purple", "tab:purple", "tab:purple", "tab:purple"]
    for i, ax in enumerate(axes_left.flat):
        dx = 0.15
        aspect = 0.65
        ax_in = ax.inset_axes([1-dx, 0, dx, dx*aspect])
        ax_in.set_xticklabels([])
        ax_in.set_yticklabels([])
        ax_in.set_xticks([])
        ax_in.set_yticks([])
        ax_in.plot(0, 0, mk[i], c=color[i], ms=9, fillstyle="none")
    plt.show()
    # plt.savefig("fig/osi_to_travelling.pdf")
    plt.close()


if __name__ == "__main__":
    
    plot_space_time_snaps()