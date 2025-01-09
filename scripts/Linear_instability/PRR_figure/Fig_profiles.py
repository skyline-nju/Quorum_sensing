import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
from add_line import add_line
from matplotlib import image as mpimg


def get_q1(q_radial, Sq):
    """ Get the first moment of s(q, t). """
    I0 = np.sum(Sq)
    I1 = np.sum(q_radial * Sq)
    q1 = I1 / I0
    return q1



def plot_profile_R_order_para():
    fig = plt.figure(figsize=(12, 4), layout="constrained")
    subfigs = fig.subfigures(1, 3, wspace=0.0001, hspace=0.0001, width_ratios=[1, 1, 1])
    ax_left = subfigs[0].subplots(3, 1, sharex=True)
    ax1, ax2, ax3 = ax_left
    ax4 = subfigs[1].subplots(1, 1)
    ax5 = subfigs[2].subplots(1, 1)

    with np.load("data/time_ave_profile/LA_CCB_profiles_p18_14.5_L480_160.npz", "rb") as data:
        x, fx = data["x"], data["fx"]
        x /= 480
        ax1.plot(x, fx[0] / 10, label=r"$S=A$", c="tab:blue")
        ax1.plot(x, fx[1] / 10, label=r"$S=B$", c="tab:red")
    
    with np.load("data/instant_profile/L480_160_Dr0.100_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_1001.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax2.plot(x/480, fields[0, 0]/5, c="tab:blue")
        ax2.plot(x/480, fields[0, 1]/5, c="tab:red")

    with np.load("data/instant_profile/L120_40_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.100_1000.npz", "rb") as data:
        t, x, fields = data["t"], data["x"], data["fields"]
        ax3.plot(x/120, fields[-1, 0]/10, c="tab:blue")
        ax3.plot(x/120, fields[-1, 1]/10, c="tab:red")

    files_sq = ["data/Sq/L1280_1280_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p4.5_14.7_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p17_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
                ]
    labels = ["G+LB", "G+CCB", "G+LB+\nCCB", "LB+CCB", "CCB"]
    mks =["o", "s", "v", "p", "D"]
    for j, fin in enumerate(files_sq):
        with np.load(fin, "r") as data:
            t, q, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            # q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        if j == 1:
            mask = t < 2.5e6
            t = t[mask]
            q1_B = q1_B[mask]
        ax4.plot(t, 2*np.pi/q1_B, mks[j], fillstyle="none", label=labels[j], ms=5)
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.legend(loc=(0.03, 0.6), fontsize="large", borderpad=0.05, frameon=True, handletextpad=0.2, handlelength=1)

    L_arr = [40, 80, 160, 320, 640]
    jm1 = [0.02447, 0.01205, 0.00574, 0.00267]
    jm2 = [0.05016, 0.02543, 0.01392, 0.007394, 0.003154]
    jm3 = [0.54276, 0.6023, 0.098, 0.04883, 0.0222]

    ax5.plot(L_arr, jm3, "v", label="G+LB+CCB", fillstyle="none", c="tab:green", ms=5)
    ax5.plot(L_arr, jm2, "p", label="LB+CCB", fillstyle="none", c="tab:red", ms=5)
    ax5.plot(L_arr[:-1], jm1, "D", label="CCB", fillstyle="none", c="tab:purple", ms=5)

    ax5.set_xscale("log")
    ax5.set_yscale("log")

    ax1.set_xlim(0, 1)
    ax1.set_ylabel(r"$\langle\rho_S\rangle_{y, t}$", fontsize="x-large")
    ax2.set_ylabel(r"$\langle\rho_S\rangle_{y}$", fontsize="x-large")
    ax3.set_ylabel(r"$\langle\rho_S\rangle_{y}$", fontsize="x-large")
    ax3.set_xlabel(r"$x/L_x$", fontsize="x-large")
    ax1.legend(loc="lower center", ncol=2, fontsize="large", borderpad=0.2, handlelength=1.5)

    ax4.set_xlabel(r"$t$", fontsize="x-large")
    ax4.set_ylabel(r"$\xi_B$", fontsize="x-large")
    ax5.set_ylabel(r"$\langle|\mathbf{J}|\rangle_t$", fontsize="x-large")
    ax5.set_xlabel(r"$L$", fontsize="x-large")
    ax5.legend(fontsize="large", loc=(0.45, 0.78), borderpad=0.2, handlelength=1)
    # ax5.set_ylim(ymax=0.07)

    label_font_size = "x-large"
    ax1.text(0.91, 0.78, "(a)", fontsize=label_font_size, transform=ax1.transAxes)
    ax2.text(0.91, 0.78, "(b)", fontsize=label_font_size, transform=ax2.transAxes)
    ax3.text(0.91, 0.78, "(c)", fontsize=label_font_size, transform=ax3.transAxes)
    ax4.text(0.91, 0.94, "(d)", fontsize=label_font_size, transform=ax4.transAxes)
    ax5.text(0.91, 0.94, "(e)", fontsize=label_font_size, transform=ax5.transAxes)

    ax2.arrow(0.78, 0.5, 0.08, 0, transform=ax2.transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax3.arrow(0.2, 0.8, 0.08, 0, transform=ax3.transAxes, width=0.04, color="k", ec="k", head_length=0.04)
    ax3.arrow(0.65, 0.8, 0.08, 0, transform=ax3.transAxes, width=0.04, color="k", ec="k", head_length=0.04)

    add_line(ax4, 0.35, 0.54, 0.98, 0.13, label=r"$0.13$", c="tab:orange", xl=0.8, yl=0.65)

    add_line(ax4, 0.3, 0.26, 0.95, 1/3, label=r"$0.33$", c="tab:blue", xl=0.55, yl=0.45)
    # add_line(ax4, 0.35, 0.7, 0.65, 1/3, label=None, c="tab:blue")
    add_line(ax4, 0.25, 0.4, 0.65, 1, label=r"$1$", c="tab:green", xl=0.35, yl=0.88)

    add_line(ax5, 0.01, 0.8, 0.95, -1, label=r"$-1$", xl=0.3, yl=0.5)
    add_line(ax5, 0.01, 0.6, 0.95, -1)
    add_line(ax5, 0.01, 0.4, 0.95, -1)

    plt.show()
    # plt.savefig("fig/profile_R_J.pdf")
    plt.close()
 


def plot_rho_v():
    fins = ["data/time_ave_profile/G_LB_profiles_p5_20_r20_L120_40.npz",
            "data/time_ave_profile/G_CCB_profiles_p7.6525_10.5675_L960_160.npz"]
    
    fig, axes = plt.subplots(3, 2, figsize=(8, 4.5), sharex="col", sharey="row", constrained_layout=True)
    rho0_arr = [20, 10]
    for col, fin in enumerate(fins):
        with np.load(fin, "rb") as data:
            x, fx = data["x"], data["fx"]
        axes[0, col].plot(x, fx[0] / rho0_arr[col], label=r"$S=A$", c="tab:blue")
        axes[0, col].plot(x, fx[1] / rho0_arr[col], label=r"$S=B$", c="tab:red")
        axes[1, col].plot(x, fx[2]/fx[0], label=r"$S=A$", c="tab:blue")
        axes[1, col].plot(x, fx[3]/fx[1], label=r"$S=B$", c="tab:red")
        axes[2, col].plot(x, fx[8] / rho0_arr[col], label=r"$S=A$", c="tab:blue")
        axes[2, col].plot(x, fx[9] / rho0_arr[col], label=r"$S=B$", c="tab:red")

    axes[0, 0].set_xlim(0, 120)
    axes[0, 1].set_xlim(0, 960)
    axes[2, 0].set_ylim(ymax=1.1)


    axes[2, 0].set_xlabel(r"$x$", fontsize="xx-large")
    axes[2, 1].set_xlabel(r"$x$", fontsize="xx-large")

    axes[0, 0].set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    axes[0, 0].set_ylabel(r"$\frac{\langle\rho_S\rangle_{y,t}}{\rho_0}$", fontsize=22)
    # axes[0, 0].set_ylabel(r"$\langle\rho_S\rangle_{y,t}\slash\rho_0$", fontsize="xx-large")

    axes[1, 0].set_ylabel(r"$\frac{\langle p_{x,S}\rangle_{y,t}}{\langle\rho_S\rangle_{y,t}}$", fontsize=22)
    # axes[1, 0].set_ylabel(r"$\langle p_{x,S}\rangle_{y,t}/\langle\rho_S\rangle_{y,t}$", fontsize="xx-large")

    axes[2, 0].set_ylabel(r"$\frac{\langle \rho_S v_S\rangle_{y,t}}{\rho_0}$", fontsize=22)
    # axes[2, 0].set_ylabel(r"$\langle \rho_S v_S\rangle_{y,t}\slash\rho_0$", fontsize="xx-large")

    axes[0, 0].set_title(r"(a) G+LB", fontsize="xx-large")
    axes[0, 1].set_title(r"(b) G+CCB", fontsize="xx-large")

    phi1 = [0.2867309, 0.49119514, 0.1376073, 2.55688951]
    axes[0, 0].axhline(phi1[0], linestyle=":", c="tab:blue")
    axes[0, 0].axhline(phi1[1], linestyle=":", c="tab:red")
    axes[0, 0].axhline(phi1[2], linestyle="--", c="tab:blue")
    axes[0, 0].axhline(phi1[3], linestyle="--", c="tab:red")

    axes[0, 0].legend(loc=(0.68, 0.3), fontsize="large")
    # axes[0, 0].text(0.01, 0.8, "(a)", fontsize="xx-large", transform=axes[0, 0].transAxes)
    # axes[0, 1].text(0.01, 0.8, "(b)", fontsize="xx-large", transform=axes[0, 1].transAxes)
    # axes[1, 0].text(0.01, 0.84, "(c)", fontsize="xx-large", transform=axes[1, 0].transAxes)
    # axes[1, 1].text(0.01, 0.84, "(d)", fontsize="xx-large", transform=axes[1, 1].transAxes)
    # axes[2, 0].text(0.01, 0.84, "(e)", fontsize="xx-large", transform=axes[2, 0].transAxes)
    # axes[2, 1].text(0.01, 0.84, "(f)", fontsize="xx-large", transform=axes[2, 1].transAxes)
    plt.show()
    # plt.savefig("fig/Fig_profile_rho_v.pdf")
    plt.close()
           

def plot_R_J_snaps():
    fig = plt.figure(figsize=(7.5, 6.5), layout="constrained")
    subfigs = fig.subfigures(2, 1, wspace=0.0001, hspace=0.01, height_ratios=[1,0.73])
    (ax1, ax2) = subfigs[0].subplots(1, 2)

    ax_snaps = subfigs[1].subplots(1, 3, gridspec_kw=dict(hspace=0, wspace=0, left=0, right=1, bottom=0., top=1))
    # ax_snaps = subfigs[1].subplots(1, 2, sharex=True)

    # fig, axes = plt.subplots(1, 4, figsize=(13, 3.5), constrained_layout=True, width_ratios=[1, 1, 1, 1])
    # ax1, ax2 = axes[:2]
    # ax_snaps = axes[2:]

    # fig, axes = plt.subplots(2, 2, figsize=(7.5, 7.5), constrained_layout=True)
    # ax1, ax2 = axes[0]
    # ax_snaps = axes[1]

    files_sq = ["data/Sq/L1280_1280_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p4.5_14.7_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
                "data/Sq/L640_640_Dr0.100_k0.70_p17_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
                ]
    labels = ["G+LB", "G+CCB", "G+LB+\nCCB", "LB+\nCCB", "CCB"]
    mks =["o", "s", "v", "p", "D"]
    for j, fin in enumerate(files_sq):
        with np.load(fin, "r") as data:
            t, q, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            # q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        if j == 1:
            mask = t < 2.5e6
            t = t[mask]
            q1_B = q1_B[mask]
        ax1.plot(t, 2*np.pi/q1_B, mks[j], fillstyle="none", label=labels[j], ms=5)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend(loc=(0.007, 0.4), fontsize="x-large", borderpad=0.05, frameon=False, handletextpad=0.2, handlelength=1)
    ax1.set_xlim(xmax=3e6)

    L_arr = [40, 80, 160, 320, 640]
    jm1 = [0.02447, 0.01205, 0.00574, 0.00267]
    jm2 = [0.05016, 0.02543, 0.01392, 0.007394, 0.003154]
    jm3 = [0.54276, 0.6023, 0.098, 0.04883, 0.0222]

    ax2.plot(L_arr, jm3, "v", label="G+LB+CCB", fillstyle="none", c="tab:green", ms=5)
    ax2.plot(L_arr, jm2, "p", label="LB+CCB", fillstyle="none", c="tab:red", ms=5)
    ax2.plot(L_arr[:-1], jm1, "D", label="CCB", fillstyle="none", c="tab:purple", ms=5)

    ax2.set_xscale("log")
    ax2.set_yscale("log")



    # ax1.set_xl(r"$t$", fontsize="x-large")
    ax1.set_title(r"(a) $\zeta_B$", fontsize="xx-large")
    ax2.set_title(r"(b) $\langle|\mathbf{J}|\rangle_t$", fontsize="xx-large")
    ax_snaps[0].set_title("(c) G+LB+CCB", fontsize="xx-large")
    ax_snaps[1].set_title("(d) LB+CCB", fontsize="xx-large")
    ax_snaps[2].set_title("(e) CCB", fontsize="xx-large")
    # ax2.set_xlabel(r"$L$", fontsize="x-large")
    ax2.legend(fontsize="x-large", loc=(0.52, 0.72), borderpad=0.15, handlelength=1)
    ax2.set_ylim(ymin=1e-3)

    label_font_size = "xx-large"
    ax1.text(0.9, 0.02, r"$t$", fontsize=label_font_size, transform=ax1.transAxes)
    ax2.text(0.9, 0.02, r"$L$", fontsize=label_font_size, transform=ax2.transAxes)

    img_names = ["G_LB_CCB_L640.png", "LB_CCB_L640.png", "CCB_L160.jpg"]
    for i, ax in enumerate(ax_snaps):
        fname = f"fig/snap/{img_names[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

    x = np.array([0.68, 0.35, 0.75, 0.1, 0.87, 0.38, 0.65, 0.4])
    y = np.array([0.65, 0.45, 0.25, 0.65, 0.45, 0.85, 0.92, 0.1])
    theta = np.array([-145, -120, 160, -40, 140, -30, 160, 145], float) * np.pi / 180
    for i in range(x.size):
        ax_snaps[1].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=ax_snaps[1].transAxes, width=0.0075, color="k", ec="k", head_length=0.03)
    
    x = np.array([0.4, 0.4, 0.6, 0.8, 0.9, 0.75, 0.9])
    y = np.array([0.9, 0.65, 0.5, 0.45, 0.2, 0.9, 0.75])
    theta = np.array([0, 0, 35, 180, 0, 180, -90], float) * np.pi / 180
    for i in range(x.size):
        ax_snaps[0].arrow(x[i], y[i], 0.05 * np.cos(theta[i]), 0.05 * np.sin(theta[i]), transform=ax_snaps[0].transAxes, width=0.0075, color="k", ec="k", head_length=0.03)
    
    
    
    add_line(ax1, 0.35, 0.54, 0.98, 0.13, label=r"$0.13$", c="tab:orange", xl=0.8, yl=0.65)

    add_line(ax1, 0.3, 0.26, 0.95, 1/3, label=r"$0.33$", c="tab:blue", xl=0.55, yl=0.45)
    # add_line(ax4, 0.35, 0.7, 0.65, 1/3, label=None, c="tab:blue")
    add_line(ax1, 0.25, 0.4, 0.65, 1, label=r"$1$", c="tab:green", xl=0.35, yl=0.88)

    add_line(ax2, 0.01, 0.85, 0.95, -1, label=r"$-1$", xl=0.3, yl=0.5)
    add_line(ax2, 0.01, 0.65, 0.95, -1)
    add_line(ax2, 0.01, 0.45, 0.8, -1)

    plt.show()
    # plt.savefig("fig/R_J_snaps.pdf", dpi=150)
    plt.close()


if __name__ == "__main__":
    # plot_R_J_snaps()
    plot_rho_v()