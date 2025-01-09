import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

f2 = "../data/profile/L10_4_Dr0.100_r200_e0.000_0.000_J1.800000_-1.800000_h0.100_1002.npz"
f0 = "../data/profile/L40_16_Dr0.100_0.100_k0.70_p80_80_r80_80_80_e0.000_-2.000_J0.10000_-0.10000_h0.100_2001.npz"
f1 = "../data/profile/L80_32_Dr0.100_0.100_k0.70_p40_40_r40_40_40_e-2.000_-2.000_J0.10000_-0.10000_h0.100_2001.npz"


def plot_rho_v_p():
    fig = plt.figure(figsize=(6*1.5, 4*1.25))
    subfigs = fig.subfigures(2, 1, wspace=0.001, hspace=0.001, height_ratios=[1, 3])
    ax_snap = subfigs[0].subplots(1, 3)
    ax_profile = subfigs[1].subplots(3, 3, sharex="col", sharey="row")

    for i, ax in enumerate(ax_snap):
        im = mpimg.imread("../data/profile/s%d.png" % i)
        # dx = Lx[col]/rho_A.size * 0.5
        # dy = Ly[col]/frame[0].shape[0] * 0.5
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
    ax_snap[1].arrow(0.5, 0.5, -0.15, 0, transform=ax_snap[1].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_snap[2].arrow(0.7, 0.5, -0.1, 0, transform=ax_snap[2].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_snap[2].arrow(0.3, 0.5, 0.1, 0, transform=ax_snap[2].transAxes, width=0.04, color="k", ec="k", head_length=0.04)

    # subfigs[0].subplots_adjust(wspace=0, hspace=0, left=-0.2, bottom=0, right=1, top=1)

    fname = [f0, f1, f2]
    frames = [-1, -1, -1]
    shifts = [0, 0, 0]
    r0 = [80, 40, 200]
    Lx = np.array([40, 80, 10])
    Ly = Lx * 0.4
    for col in range(3):
        with np.load(fname[col], "r") as data:
            x = data["x"]
            frame = data["fields"][frames[col]]
            rho_A = np.mean(frame[0], axis=0)
            rho_B = np.mean(frame[1], axis=0)
            # p_A = np.mean(frame[2], axis=0) / rho_A
            # p_B = np.mean(frame[3], axis=0) / rho_B
            p_A = np.mean(frame[2], axis=0) / r0[col]
            p_B = np.mean(frame[3], axis=0) / r0[col]
            v_A = np.mean(frame[6], axis=0)
            v_B = np.mean(frame[7], axis=0)
            vrho_A = np.mean(frame[6] * frame[0], axis=0)
            vrho_B = np.mean(frame[7] * frame[1], axis=0)
        ax1, ax2, ax3 = ax_profile[0, col], ax_profile[1, col], ax_profile[2, col]
        dx = shifts[col]
        ax1.plot(x, np.roll(rho_A/r0[col], dx), c="tab:blue", label="S=A")
        ax1.plot(x, np.roll(rho_B/r0[col], dx), c="tab:red", label="S=B")
        # ax1_twin = ax1.twinx()
        # ax1_twin.plot(x, np.roll(p_A, dx), c="tab:blue", linestyle="dashed")
        # ax1_twin.plot(x, np.roll(p_B, dx), c="tab:red", linestyle="dashed")
        # ax1.set_ylim(0)
        ax2.plot(x, np.roll(v_A, dx), c="tab:blue")
        ax2.plot(x, np.roll(v_B, dx), c="tab:red")

        ax3.plot(x, np.roll(p_A, dx), c="tab:blue")
        ax3.plot(x, np.roll(p_B, dx), c="tab:red")

        ax1.set_xlim(0, Lx[col])
        # ax1.set_xlim(0, 1)

        ax3.set_xlabel(r"$x$", fontsize="x-large")

    ax_profile[0, 0].legend(fontsize="large")


    fs = "x-large"
    ax_profile[0, 0].set_ylabel(r"$\langle\rho_S\rangle_y/\rho_0$", fontsize=fs)
    ax_profile[1, 0].set_ylabel(r"$\langle v_S\rangle_y$", fontsize=fs)
    # ax_profile[2, 0].set_ylabel(r"$\langle p_{x,S}\rangle_y/\langle\rho_S\rangle_y$", fontsize=fs)
    ax_profile[2, 0].set_ylabel(r"$\langle p_{x,S}\rangle_y/\rho_0$", fontsize=fs)


    subfigs[0].subplots_adjust(wspace=0.08, hspace=0., left=0.065, bottom=0.05, right=0.99, top=1.05)
    subfigs[1].subplots_adjust(wspace=0.08, hspace=0, left=0.065, bottom=0.11, right=0.99, top=1.01)

    label_font_size ="x-large"
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    ax_snap[0].text(0.015, 0.82, "(a)", fontsize=label_font_size, transform=ax_snap[0].transAxes, backgroundcolor="w", bbox=bbox)
    ax_snap[1].text(0.015, 0.82, "(b)", fontsize=label_font_size, transform=ax_snap[1].transAxes, backgroundcolor="w", bbox=bbox)
    ax_snap[2].text(0.015, 0.82, "(c)", fontsize=label_font_size, transform=ax_snap[2].transAxes, backgroundcolor="w", bbox=bbox)

    ax_profile[0, 0].text(0.015, 0.82, "(d)", fontsize=label_font_size, transform=ax_profile[0, 0].transAxes)
    ax_profile[0, 1].text(0.015, 0.82, "(e)", fontsize=label_font_size, transform=ax_profile[0, 1].transAxes)
    ax_profile[0, 2].text(0.015, 0.82, "(f)", fontsize=label_font_size, transform=ax_profile[0, 2].transAxes)



    plt.show()
    # plt.savefig("fig/snap_profile_quasi_1D.pdf")
    plt.close()



def plot_rho_p():
    fig = plt.figure(figsize=(6*1.5, 3*1.25))
    subfigs = fig.subfigures(2, 1, wspace=0.001, hspace=0.001, height_ratios=[1, 2])
    ax_snap = subfigs[0].subplots(1, 3)
    ax_profile = subfigs[1].subplots(2, 3, sharex="col", sharey="row")

    for i, ax in enumerate(ax_snap):
        im = mpimg.imread("../data/profile/s%d.png" % i)
        # dx = Lx[col]/rho_A.size * 0.5
        # dy = Ly[col]/frame[0].shape[0] * 0.5
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
    ax_snap[1].arrow(0.5, 0.5, -0.15, 0, transform=ax_snap[1].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_snap[2].arrow(0.7, 0.5, -0.1, 0, transform=ax_snap[2].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax_snap[2].arrow(0.3, 0.5, 0.1, 0, transform=ax_snap[2].transAxes, width=0.04, color="k", ec="k", head_length=0.04)

    # subfigs[0].subplots_adjust(wspace=0, hspace=0, left=-0.2, bottom=0, right=1, top=1)

    fname = [f0, f1, f2]
    frames = [-1, -1, -1]
    shifts = [0, 0, 0]
    r0 = [80, 40, 200]
    Lx = np.array([40, 80, 10])
    Ly = Lx * 0.4
    for col in range(3):
        with np.load(fname[col], "r") as data:
            x = data["x"]
            frame = data["fields"][frames[col]]
            rho_A = np.mean(frame[0], axis=0)
            rho_B = np.mean(frame[1], axis=0)
            # p_A = np.mean(frame[2], axis=0) / rho_A
            # p_B = np.mean(frame[3], axis=0) / rho_B
            p_A = np.mean(frame[2], axis=0) / r0[col]
            p_B = np.mean(frame[3], axis=0) / r0[col]
            v_A = np.mean(frame[6], axis=0)
            v_B = np.mean(frame[7], axis=0)
            vrho_A = np.mean(frame[6] * frame[0], axis=0)
            vrho_B = np.mean(frame[7] * frame[1], axis=0)
        ax1, ax2 = ax_profile[0, col], ax_profile[1, col]
        dx = shifts[col]
        ax1.plot(x, np.roll(rho_A/r0[col], dx), c="tab:blue", label="S=A")
        ax1.plot(x, np.roll(rho_B/r0[col], dx), c="tab:red", label="S=B")
        # ax1_twin = ax1.twinx()
        # ax1_twin.plot(x, np.roll(p_A, dx), c="tab:blue", linestyle="dashed")
        # ax1_twin.plot(x, np.roll(p_B, dx), c="tab:red", linestyle="dashed")
        # ax1.set_ylim(0)

        ax2.plot(x, np.roll(p_A, dx), c="tab:blue")
        ax2.plot(x, np.roll(p_B, dx), c="tab:red")

        ax1.set_xlim(0, Lx[col])
        # ax1.set_xlim(0, 1)

        ax2.set_xlabel(r"$x$", fontsize="x-large")

    ax_profile[0, 0].legend(fontsize="large")


    fs = "x-large"
    ax_profile[0, 0].set_ylabel(r"$\langle\rho_S\rangle_y/\rho_0$", fontsize=fs)
    # ax_profile[1, 0].set_ylabel(r"$\langle v_S\rangle_y$", fontsize=fs)
    # ax_profile[2, 0].set_ylabel(r"$\langle p_{x,S}\rangle_y/\langle\rho_S\rangle_y$", fontsize=fs)
    ax_profile[1, 0].set_ylabel(r"$\langle p_{x,S}\rangle_y/\rho_0$", fontsize=fs)


    subfigs[0].subplots_adjust(wspace=0.08, hspace=0., left=0.065, bottom=0.05, right=0.99, top=1.05)
    subfigs[1].subplots_adjust(wspace=0.08, hspace=0, left=0.065, bottom=0.16, right=0.99, top=1.01)

    label_font_size ="x-large"
    bbox=dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.08")
    ax_snap[0].text(0.015, 0.82, "(a)", fontsize=label_font_size, transform=ax_snap[0].transAxes, backgroundcolor="w", bbox=bbox)
    ax_snap[1].text(0.015, 0.82, "(b)", fontsize=label_font_size, transform=ax_snap[1].transAxes, backgroundcolor="w", bbox=bbox)
    ax_snap[2].text(0.015, 0.82, "(c)", fontsize=label_font_size, transform=ax_snap[2].transAxes, backgroundcolor="w", bbox=bbox)

    ax_profile[0, 0].text(0.015, 0.82, "(d)", fontsize=label_font_size, transform=ax_profile[0, 0].transAxes)
    ax_profile[0, 1].text(0.015, 0.82, "(e)", fontsize=label_font_size, transform=ax_profile[0, 1].transAxes)
    ax_profile[0, 2].text(0.015, 0.82, "(f)", fontsize=label_font_size, transform=ax_profile[0, 2].transAxes)
    ax_profile[1, 0].set_yticks([-4, 0, 4])



    plt.show()
    # plt.savefig("fig/snap_profile_quasi_1D.pdf")
    plt.close()


def plot_rho_p_space_time():
    fig = plt.figure(figsize=(10*1.5, 3*1.3), layout="constrained")
    subfigs = fig.subfigures(1, 5, wspace=0.001, hspace=0.05, width_ratios=[1.08, 1, 1, 1.03, 0.97])

    Lx_arr = np.array([40, 80, 10])
    r0 = [80, 40, 200]

    f_profiles = [f0, f1, f2]
    fs ="xx-large"

    ax_rho_x = []
    ax_p_x = []

    ins_types = ["(a) LSI", "(b) LOI", "(c) SOI"]
    for j, subfig in enumerate(subfigs[:3]):
        ax = subfig.subplots(3, 1, sharex=True, height_ratios=[1.2,  1, 1])
        im = mpimg.imread("../data/profile/s%d.png" % j)
        ax[0].imshow(im, extent=[0, Lx_arr[j], 0, Lx_arr[j]/2.5])
        # ax[0].set_xticklabels([])
        ax[0].set_yticks([])
        ax[0].set_title(ins_types[j], fontsize=fs)
        ax_rho_x.append(ax[1])
        ax_p_x.append(ax[2])

        if j == 1:
            ax[0].arrow(0.5, 0.5, -0.15, 0, transform=ax[0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
        elif j == 2:
            ax[0].arrow(0.7, 0.5, -0.1, 0, transform=ax[0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
            ax[0].arrow(0.3, 0.5, 0.1, 0, transform=ax[0].transAxes, width=0.04, color="k", ec="k", head_length=0.04)
        else:
            ax[0].arrow(0.05, 0.1, 0.15, 0, transform=ax[0].transAxes, width=0.03, color="tab:grey", ec="tab:grey", head_length=0.03)
            ax[0].arrow(0.05, 0.1, 0, 0.15*2.5, transform=ax[0].transAxes, width=0.03/2.5, color="tab:grey", ec="tab:grey", head_length=0.03*2.5)
            ax[0].text(0.075, 0.5, r"$y$", transform=ax[0].transAxes, color="k", fontsize="xx-large")
            ax[0].text(0.2, 0.175, r"$x$", transform=ax[0].transAxes, color="k", fontsize="xx-large")
    ylim_rho = (0, 6.5) 
    ylim_p = (-5.5, 5.5)
    for i, f_profile in enumerate(f_profiles):
        with np.load(f_profile, "r") as data:
            x = data["x"]
            frame = data["fields"][-1]
            rho_A = np.mean(frame[0], axis=0)
            rho_B = np.mean(frame[1], axis=0)
            # p_A = np.mean(frame[2], axis=0) / rho_A
            # p_B = np.mean(frame[3], axis=0) / rho_B
            p_A = np.mean(frame[2], axis=0) / r0[i]
            p_B = np.mean(frame[3], axis=0) / r0[i]
            v_A = np.mean(frame[6], axis=0)
            v_B = np.mean(frame[7], axis=0)
            vrho_A = np.mean(frame[6] * frame[0], axis=0)
            vrho_B = np.mean(frame[7] * frame[1], axis=0)
        
        if i > 0:
            # ax_rho_x[i].set_yticklabels([])
            # ax_p_x[i].set_yticklabels([])
            pass
        else:
            ax_rho_x[i].set_ylabel(r"$\langle\rho_S\rangle_{y}/\rho_0$", fontsize=fs)
            ax_p_x[i].set_ylabel(r"$\langle p_{x,S}\rangle_{y}/\rho_0$", fontsize=fs)
        for ax in ax_p_x:
            ax.set_xlabel(r"$x$", fontsize=fs)

        ax_rho_x[i].plot(x, rho_A/r0[i], c="tab:blue", label=r"$S=A$")
        ax_rho_x[i].plot(x, rho_B/r0[i], c="tab:red", label=r"$S=B$")
        

        ax_p_x[i].plot(x, p_A, c="tab:blue", label="S=A")
        ax_p_x[i].plot(x, p_B, c="tab:red", label="S=B")

        ax_rho_x[i].set_ylim(ylim_rho)
        ax_p_x[i].set_ylim(ylim_p)

        ax_rho_x[0].legend(ncols=2, fontsize="large")


   
    ax7 = subfigs[3].subplots()
    ax8 = subfigs[4].subplots()

    with np.load("data/space_time/L80_16_Dr0.100_0.100_k0.70_p40_40_r40_40_40_e-2.000_-2.000_J0.10000_-0.10000_h0.100_2100.npz", "r") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        rhoA = fields[:, 0] / 40
        im1 = ax7.imshow(rhoA, origin="lower", aspect="auto", extent=[0, 80, 0, t[-1]])
    
    with np.load("data/space_time/L10_4_Dr0.100_r200_e0.000_0.000_J1.800000_-1.800000_h0.100_1005.npz", "r") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        rhoA = fields[:, 0] / 200
        im2 = ax8.imshow(rhoA, origin="lower", aspect="auto", extent=[0, 10, 0, t[-1]], vmax=7)
    ax7.set_xlabel(r"$x$", fontsize=fs)
    ax8.set_xlabel(r"$x$", fontsize=fs)
    ax7.set_ylabel(r"$t$", fontsize=fs)
    ax7.set_title("(d) LOI", fontsize=fs)
    ax8.set_title("(e) SOI", fontsize=fs)
    # ax8.set_ylabel(r"$t$", fontsize=fs)

    cb1 = subfigs[3].colorbar(im1, shrink=0.5, ax=ax7, location="right",  extend="max")

    # cb1.set_label(r"$\langle \rho_B(\mathbf{r},t)\rangle_y /\rho_0$", fontsize="x-large")

    subfigs[4].colorbar(im2, shrink=0.5, ax=ax8, location="right",  extend="max")

    subfigs[3].text(0.88, 0.77, r"$\frac{\langle \rho_A\rangle_y}{\rho_0}$", fontsize="xx-large")
    subfigs[4].text(0.83, 0.77, r"$\frac{\langle \rho_A\rangle_y}{\rho_0}$", fontsize="xx-large")

    plt.show()
    # plt.savefig("fig/snap_profile_space_time_quasi_1D.pdf")
    plt.close()

if __name__ == "__main__":
    plot_rho_p_space_time()