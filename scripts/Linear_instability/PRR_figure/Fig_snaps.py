import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import image as mpimg
import os
import sys


def CCB_phase_separation():
    fig = plt.figure(figsize=(12, 8.7), layout="constrained")
    subfigs = fig.subfigures(3, 1, wspace=0.001, hspace=0.001, height_ratios=[160, 160, 640])

    axes0 = subfigs[0].subplots(1, 3, sharex=True, sharey=True)
    axes1 = subfigs[1].subplots(1, 3, sharex=True, sharey=True)
    axes2 = subfigs[2].subplots(1, 2, sharex=True, sharey=True)

    axes = [i for i in axes0]
    for i in axes1:
        axes.append(i)
    for i in axes2:
        axes.append(i)

    fnames = ["G_CCB.png", "G_LA_CCB.png", "LA_CCB.png",
              "G_LB_CCB.png", "LB_CCB_0.png", "LB_CCB_1.png",
              "G_LB_CCB_L640.png", "LB_CCB_L640.png"]
    titles = ["(a) G+CCB", "(b) G+LA+CCB", "(c) LA+CCB", "(d) G+LB+CCB",
              "(e) LB+CCB: initial configuration", "(f) LB+CCB: stready state",
              "(g) G+LB+CCB", "(h) LB+CCB"]

    for i, ax in enumerate(axes):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        # ax.set_title(titles[i], fontsize="x-large")

    ax_in = axes1[2].inset_axes([1-0.25, 1-0.25, 0.25, 0.25])
    fname = f"fig/snap/LB_CCB_L120_40.png"
    im = mpimg.imread(fname)
    ax_in.imshow(im)
    ax_in.set_xticklabels([])
    ax_in.set_yticklabels([])
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.set_aspect("equal")

    fs = "xx-large"
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    axes0[0].text(0.01, 0.83, "(a) G+CCB", fontsize=fs, transform=axes0[0].transAxes, backgroundcolor="w", bbox=bbox)
    axes0[1].text(0.01, 0.83, "(b) G+LA+CCB", fontsize=fs, transform=axes0[1].transAxes, backgroundcolor="w", bbox=bbox)
    axes0[2].text(0.01, 0.83, "(c) LA+CCB", fontsize=fs, transform=axes0[2].transAxes, backgroundcolor="w", bbox=bbox)
    axes1[0].text(0.01, 0.83, "(d) G+LB+CCB", fontsize=fs, transform=axes1[0].transAxes, backgroundcolor="w", bbox=bbox)
    axes1[1].text(0.01, 0.83, r"(e) LB+CCB: $t=0$", fontsize=fs, transform=axes1[1].transAxes, backgroundcolor="w", bbox=bbox)
    axes1[2].text(0.01, 0.8, r"(f) LB+CCB: $t=10^5$", fontsize=fs, transform=axes1[2].transAxes, backgroundcolor="w", bbox=bbox)
    axes2[0].text(0.01, 0.955, "(g) G+LB+CCB", fontsize=fs, transform=axes2[0].transAxes, backgroundcolor="w", bbox=bbox)
    axes2[1].text(0.01, 0.955, "(h) LB+CCB", fontsize=fs, transform=axes2[1].transAxes, backgroundcolor="w", bbox=bbox)

    axes1[0].arrow(0.75, 0.5, 0.15, 0, transform=axes1[0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    axes1[1].arrow(0.075, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.2, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.325, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.45, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.575, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.7, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.825, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.935, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_in.arrow(0.25, 0.5, 0.16, 0, transform=ax_in.transAxes, width=0.08, color="k", ec="k", head_length=0.08)
    ax_in.arrow(0.7, 0.5, 0.16, 0, transform=ax_in.transAxes, width=0.08, color="k", ec="k", head_length=0.08)

    x = np.array([0.9, 0.64, 0.4])
    y = np.array([0.5, 0.48, 0.5])
    theta = np.array([-30, -40, 20,], float) * np.pi / 180
    for i in range(x.size):
        axes1[2].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=axes1[2].transAxes, width=0.015, color="k", ec="k", head_length=0.03)

    x = np.array([0.68, 0.35, 0.75, 0.1, 0.87, 0.38, 0.65, 0.4])
    y = np.array([0.65, 0.45, 0.25, 0.65, 0.45, 0.85, 0.92, 0.1])
    theta = np.array([-145, -120, 160, -40, 140, -30, 160, 145], float) * np.pi / 180
    for i in range(x.size):
        axes2[1].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=axes2[1].transAxes, width=0.0075, color="k", ec="k", head_length=0.03)
    
    x = np.array([0.4, 0.4, 0.6, 0.8, 0.9, 0.75, 0.9])
    y = np.array([0.9, 0.65, 0.5, 0.45, 0.2, 0.9, 0.75])
    theta = np.array([0, 0, 35, 180, 0, 180, -90], float) * np.pi / 180
    for i in range(x.size):
        axes2[0].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=axes2[0].transAxes, width=0.0075, color="k", ec="k", head_length=0.03)
    
    plt.show()
    # plt.savefig("fig/snaps_PS_CCB.pdf")
    plt.close()


def CCB_phase_separation_profiles():
    fig = plt.figure(figsize=(12, 5.5), layout="constrained")
    subfigs = fig.subfigures(3, 1, wspace=0.001, hspace=0.001, height_ratios=[160, 160, 190])

    axes0 = subfigs[0].subplots(1, 3, sharex=True, sharey=True)
    axes1 = subfigs[1].subplots(1, 3, sharex=True, sharey=True)
    axes2 = subfigs[2].subplots(1, 3, sharex=True)

    axes = [i for i in axes0]
    for i in axes1:
        axes.append(i)
    # for i in axes2:
        # axes.append(i)

    fnames = ["G_CCB.png", "G_LA_CCB.png", "LA_CCB.png",
              "G_LB_CCB.png", "LB_CCB_0.png", "LB_CCB_1.png"]
    titles = ["(a) G+CCB", "(b) G+LA+CCB", "(c) LA+CCB", "(d) G+LB+CCB",
              "(e) LB+CCB: initial configuration", "(f) LB+CCB: stready state"]

    for i, ax in enumerate(axes):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        # ax.set_title(titles[i], fontsize="x-large")

    ax_in = axes1[2].inset_axes([1-0.25, 1-0.25, 0.25, 0.25])
    fname = f"fig/snap/LB_CCB_L120_40.png"
    im = mpimg.imread(fname)
    ax_in.imshow(im)
    ax_in.set_xticklabels([])
    ax_in.set_yticklabels([])
    ax_in.set_xticks([])
    ax_in.set_yticks([])
    ax_in.set_aspect("equal")

    fs = "xx-large"
    axes0[0].set_title("(a) G+CCB", fontsize=fs)
    axes0[1].set_title("(b) G+LA+CCB", fontsize=fs)
    axes0[2].set_title("(c) LA+CCB", fontsize=fs)
    axes1[0].set_title("(d) G+LB+CCB", fontsize=fs)
    axes1[1].set_title(r"(e) LB+CCB: $t=0$", fontsize=fs)
    axes1[2].set_title(r"(f) LB+CCB: $t=10^5$", fontsize=fs)
    # bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    # axes0[0].text(0.01, 0.83, "(a) G+CCB", fontsize=fs, transform=axes0[0].transAxes, backgroundcolor="w", bbox=bbox)
    # axes0[1].text(0.01, 0.83, "(b) G+LA+CCB", fontsize=fs, transform=axes0[1].transAxes, backgroundcolor="w", bbox=bbox)
    # axes0[2].text(0.01, 0.83, "(c) LA+CCB", fontsize=fs, transform=axes0[2].transAxes, backgroundcolor="w", bbox=bbox)
    # axes1[0].text(0.01, 0.83, "(d) G+LB+CCB", fontsize=fs, transform=axes1[0].transAxes, backgroundcolor="w", bbox=bbox)
    # axes1[1].text(0.01, 0.83, r"(e) LB+CCB: $t=0$", fontsize=fs, transform=axes1[1].transAxes, backgroundcolor="w", bbox=bbox)
    # axes1[2].text(0.01, 0.8, r"(f) LB+CCB: $t=10^5$", fontsize=fs, transform=axes1[2].transAxes, backgroundcolor="w", bbox=bbox)
    # axes2[0].text(0.01, 0.955, "(g) G+LB+CCB", fontsize=fs, transform=axes2[0].transAxes, backgroundcolor="w", bbox=bbox)
    # axes2[1].text(0.01, 0.955, "(h) LB+CCB", fontsize=fs, transform=axes2[1].transAxes, backgroundcolor="w", bbox=bbox)

    axes1[0].arrow(0.75, 0.5, 0.15, 0, transform=axes1[0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    axes1[1].arrow(0.075, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.2, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.325, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.45, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.575, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.7, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.825, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    axes1[1].arrow(0.935, 0.5, 0.04, 0, transform=axes1[1].transAxes, width=0.02, color="k", ec="k", head_length=0.02)
    ax_in.arrow(0.25, 0.5, 0.16, 0, transform=ax_in.transAxes, width=0.08, color="k", ec="k", head_length=0.08)
    ax_in.arrow(0.7, 0.5, 0.16, 0, transform=ax_in.transAxes, width=0.08, color="k", ec="k", head_length=0.08)

    x = np.array([0.9, 0.64, 0.4])
    y = np.array([0.5, 0.48, 0.5])
    theta = np.array([-30, -40, 20,], float) * np.pi / 180
    for i in range(x.size):
        axes1[2].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=axes1[2].transAxes, width=0.015, color="k", ec="k", head_length=0.03)

    ax1, ax2, ax3 = axes2
    with np.load("data/time_ave_profile/LA_CCB_profiles_p18_14.5_L480_160.npz", "rb") as data:
        x, fx = data["x"], data["fx"]
        x /= 480
        ax1.plot(x, fx[0] / 10, label=r"$S=A$", c="tab:blue")
        ax1.plot(x, fx[1] / 10, label=r"$S=B$", c="tab:red")
        ax1.set_xlim(0, 1)
        ax1.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax1.set_title(r"(g) $\langle\rho_S\rangle_{y,t}/\rho_0 $ for LA+CCB", fontsize=fs)
    
    with np.load("data/instant_profile/L480_160_Dr0.100_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_1001.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax2.plot(x/480, fields[0, 0]/5, c="tab:blue")
        ax2.plot(x/480, fields[0, 1]/5, c="tab:red")
        ax2.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax2.set_title(r"(h) $\langle\rho_S\rangle_{y}/\rho_0 $ for G+LB+CCB", fontsize=fs)


    with np.load("data/instant_profile/L120_40_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.100_1000.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax3.plot(x/120, fields[-1, 0]/10, c="tab:blue")
        ax3.plot(x/120, fields[-1, 1]/10, c="tab:red")
        ax3.set_xlabel(r"$x/L_x$", fontsize="xx-large")
        ax3.set_title(r"(i) $\langle\rho_S\rangle_{y}/\rho_0 $ for LB+CCB", fontsize=fs)

    ax2.arrow(0.78, 0.5, 0.08, 0, transform=ax2.transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax3.arrow(0.2, 0.8, 0.08, 0, transform=ax3.transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax3.arrow(0.65, 0.8, 0.08, 0, transform=ax3.transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax1.legend(loc="lower center", ncol=2, fontsize="large", borderpad=0.2, handlelength=1.5)

    # plt.show()
    plt.savefig("fig/snaps_PS_CCB_profiles.pdf")
    plt.close()


def Eq_phase_separation():
    fig, axes = plt.subplots(2, 3, figsize=(10, 2.8), sharex=True, sharey=True, constrained_layout=True)

    fnames = ["G_LA.png", "G_LAB.png", "G_LA_LAB.png", "G_LA_LB.png", "LA_LB.png", "LA_LB_LAB.png"]
    titles = ["(a) G+LA", "(b) G+LAB", "(c) G+LA+LAB", "(d) G+LA+LB",
              "(e) LA+LB", "(f) LA+LB+LAB"]

    for i, ax in enumerate(axes.flat):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(titles[i], fontsize="x-large")
    
    plt.show()
    # plt.savefig("fig/snaps_PS_Eq.pdf")
    plt.close()


def w_wo_rep():
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.35), sharex=True, sharey=True, constrained_layout=True)
    fnames = ["wo_rep.png", "w_rep.png"]
    titles = ["(a) Without repulsion", "(b) With repulsion"]
    
    for i, ax in enumerate(axes):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(titles[i], fontsize="xx-large")
    # plt.show()
    plt.savefig("fig/snaps_wo_w_rep.pdf")
    plt.close()


def plot_snaps_Dr1():
    fig, axes = plt.subplots(1, 5, sharex=True, sharey=True, constrained_layout=True, figsize=(8, 1.85))

    t_arr = np.array([0, 200, 300, 600, 4800])
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    for i, t in enumerate(t_arr):
        fname = f"fig/snap_Dr1/t={t:d}.png"
        ax = axes[i]
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        title = r"%s $t=%d$" % (labels[i], t)
        ax.set_title(title, fontsize="xx-large")
    plt.show()
    # plt.savefig("fig/SB_instability_Dr1.pdf")
    plt.close()


if __name__ == "__main__":
    # CCB_phase_separation_profiles()
    # Eq_phase_separation()
    # plot_snaps_Dr1()

    # w_wo_rep()
    plot_snaps_Dr1()


