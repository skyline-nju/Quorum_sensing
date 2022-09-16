import numpy as np
import matplotlib.pyplot as plt
from corr_length import get_Lx, get_para


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_center(x, y, Lx, inversed=False):
    theta = x / Lx * np.pi * 2
    px = np.mean(np.cos(theta) * y)
    py = np.mean(np.sin(theta) * y)
    theta_m = np.arctan2(py, px)
    if inversed:
        theta_m += np.pi
    if theta_m < 0:
        theta_m += np.pi * 2
    elif theta_m >= np.pi * 2:
        theta_m -= np.pi * 2
    xm = theta_m / (np.pi * 2) * Lx
    return xm


def shift_profile(y, xc, Lx):
    shift = int(Lx / 2 - np.around(xc))
    y_new = np.roll(y, shift)
    return y_new


def plot_instant_profile(t, x, frame, normed=False, para=None, r1=None, r2=None):
    Lx = x[0] + x[-1]
    fig, axes = plt.subplots(5, 1, figsize=(6, 6), constrained_layout=True, sharex=True)
    
    fields_x = np.array([np.mean(i, axis=0) for i in frame])
    
    if normed:
        fields_x[2] /= fields_x[0]
        fields_x[4] /= fields_x[0]
        fields_x[3] /= fields_x[1]
        fields_x[5] /= fields_x[1]
    
    for i, ax in enumerate(axes[:4]):
        ax.plot(x, fields_x[2*i], label="A")
        ax.plot(x, fields_x[2*i+1], label="B")

    rv_A = np.mean(frame[0] * frame[6], axis=0)
    rv_B = np.mean(frame[1] * frame[7], axis=0)
    axes[4].plot(x, rv_A)
    axes[4].plot(x, rv_B)
    axes[0].set_xlim(0, Lx)

    axes[-1].set_xlabel(r"$x$")
    ylabels = [r"$\langle \rho \rangle_y$",
               r"$\langle \mathbf{m}_x\rangle_y$",
               r"$\langle \mathbf{m}_y\rangle_y$",
               r"$\langle v\rangle_y$",
               r"$\langle \rho v\rangle_y$"]
    axes[0].legend()
    for i, ax in enumerate(axes):
        ax.set_ylabel(ylabels[i])
    
    if r1 is not None:
        axes[0].axhline(r1[0], linestyle="dotted")
        axes[0].axhline(r1[1], linestyle="dotted")
    if r2 is not None:
        axes[0].axhline(r2[0], linestyle="dashed")
        axes[0].axhline(r2[1], linestyle="dashed")
    if para is None:
        title = r"$t=%g$" % t
    else:
        title = "$L_x=%g,L_y=%g,D_r=%g,\\eta=%g,\\eta_{AB}=-\\eta_{BA}=%g,\\rho_0=\\bar{\\rho}_{A,B}=%g,$\n$\\phi_A=%g,\\phi_B=%g,t=%g$" % (
            para["Lx"], para["Ly"], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para['pA'], para['pB'], t
        )
    fig.suptitle(title)
    plt.show()
    plt.close()


def get_time_ave_density_profile(fin, beg_frame, end_frame=None, dt=0.01, inversed_B=True):
    Lx = get_Lx(fin)
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    if end_frame is None:
        end_frame = t.size
    n_frames = end_frame - beg_frame
    rho_x_A = np.zeros((n_frames, x.size))
    rho_x_B = np.zeros_like(rho_x_A)
    for i, i_frame in enumerate(np.arange(beg_frame, end_frame)):
        rho_A = np.mean(fields[i_frame][0], axis=0)
        rho_B = np.mean(fields[i_frame][1], axis=0)
        xm_A = get_center(x, rho_A, Lx)
        xm_B = get_center(x, rho_B, Lx, inversed=inversed_B)
        rho_x_A[i] = shift_profile(rho_A, xm_A, Lx)
        rho_x_B[i] = shift_profile(rho_B, xm_B, Lx)
    rA = np.mean(rho_x_A, axis=0)
    rB = np.mean(rho_x_B, axis=0)
    return x, rA, rB


def plot_time_ave_density_profile_SP():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    f2 = f"{prefix}/L120_40_Dr0.100_k0.70_p6.4825_27.115_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    f3 = f"{prefix}/L120_40_Dr0.100_k0.70_p7.965_21.23_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    # dt = 0.01
    # Lx = 120
    x, rA_1, rB_1 = get_time_ave_density_profile(f1, 102, 360)
    x, rA_2, rB_2 = get_time_ave_density_profile(f2, 50, 200)
    x, rA_3, rB_3 = get_time_ave_density_profile(f3, 50, 75)

    rA = [rA_1, rA_2, rA_3]
    rB = [rB_1, rB_2, rB_3]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)

    labels = [r"$\phi_A=5,\phi_B=33$", r"$\phi_A=6.4825,\phi_B=27.115$", r"$\phi_A=7.965,\phi_B=21.23$"]
    for i in range(len(rA)):
        ax1.plot(x, rA[i], label=labels[i])
        ax2.plot(x, rB[i], label=labels[i])
    ax1.set_ylabel(r"$\langle\rho_A(\mathbf{r}, t)\rangle_{y,t}$", fontsize="large")
    ax2.set_ylabel(r"$\langle\rho_B(\mathbf{r}, t)\rangle_{y,t}$", fontsize="large")
    ax2.set_xlabel(r"$x$", fontsize="large")
    ax1.legend()
    ax2.legend()
    ax1.set_xlim(0, 120)

    title = r"$L_x=120,L_y=40,D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=2,\rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title)
    plt.show()
    plt.close()


if __name__ == "__main__":
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    # f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    # prefix = f"{root_tahmineh}/QS_MPI/L480_160/coarse_grain_dx4"
    # f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p40_40_r40_e-0.800_J1.000_-1.000_1000.npz"
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L240_80_r10_10_10/coarse_grain_dx1"
    # f1 = f"{prefix}/L240_80_Dr0.100_k0.70_p10_10_r10_10_10_e-2.000_J0.500_-0.500_2000.npz"
    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10/coarse_grain_dx1"
    f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p12.5_10_r10_10_10_e-2.000_J0.500_-0.500_21010.npz"
    Lx = get_Lx(f1)
    para = get_para(f1)
    dt = 0.01
    with np.load(f1, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]

    r1=[3.88, 1.2]
    r2=[10.54, 21.65]
    for i_frame in range(t.size):
        print("frame %d/%d" % (i_frame, t.size))
        plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True, para=para, r1=r1, r2=r2)
    
    beg = 110
    x, rA_1, rB_1 = get_time_ave_density_profile(f1, beg, t.size, inversed_B=False)
    
    lineA, = plt.plot(x, rA_1)
    lineB, = plt.plot(x, rB_1)

    x_edge_1 = [50, 420]
    x_edge_2 = [200, 280]

    plt.axvline(x_edge_1[0], linestyle="dashed", c="tab:grey")
    plt.axvline(x_edge_1[1], linestyle="dashed", c="tab:grey")
    plt.axvline(x_edge_2[0], linestyle="dotted", c="tab:grey")
    plt.axvline(x_edge_2[1], linestyle="dotted", c="tab:grey")

    mask1 = np.logical_or(x <= x_edge_1[0], x > x_edge_1[1])
    r1 = [np.mean(rA_1[mask1]), np.mean(rB_1[mask1])]
    mask2 = np.logical_and(x > x_edge_2[0], x <= x_edge_2[1])
    r2 = [np.mean(rA_1[mask2]), np.mean(rB_1[mask2])]

    plt.axhline(r1[0], c=lineA.get_c(), linestyle="dashed", label="%.4f" % r1[0])
    plt.axhline(r1[1], c=lineB.get_c(), linestyle="dashed", label="%.4f" % r1[1])

    plt.axhline(r2[0], c=lineA.get_c(), linestyle="dotted", label="%.4f" % r2[0])
    plt.axhline(r2[1], c=lineB.get_c(), linestyle="dotted", label="%.4f" % r2[1])

    plt.legend()

    phiA = 7.5
    phiB = 10
    
    l1_A = Lx * (phiA - r2[0]) / (r1[0] - r2[0])
    l1_B = Lx * (phiB - r2[1]) / (r1[1] - r2[1])

    print(l1_A, l1_B)

    plt.axvline(l1_A/2, c="k", linestyle="dashdot")

    drhoA = r2[0] - r1[0]
    drhoB = r2[1] - r1[1]

    dl = -80
    phiA_new = phiA + drhoA * dl / Lx
    phiB_new = phiB + drhoB * dl / Lx
    print(phiA_new, phiB_new)

    plt.xlim(0, Lx)
    plt.show()
    plt.close()


    # beg = 50
    # f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_21010.npz"
    # x, rA_1, rB_1 = get_time_ave_density_profile(f1, beg, inversed_B=False)

    # f2 = f"{prefix}/L480_160_Dr0.100_k0.70_p8.625_13.36_r10_10_10_e-2.000_J0.500_-0.500_27510.npz"
    # x, rA_2, rB_2 = get_time_ave_density_profile(f2, beg, inversed_B=False)

    # fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    # ax1.plot(x, rA_1)
    # ax1.plot(x, rA_2)

    # ax2.plot(x, rB_1)
    # ax2.plot(x, rB_2)

    # plt.show()
    # plt.close()