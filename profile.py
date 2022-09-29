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

    axes[-1].set_xlabel(r"$x$", fontsize="large")
    ylabels = [r"$\langle \rho \rangle_y$",
               r"$\langle \mathbf{m}_x\rangle_y$",
               r"$\langle \mathbf{m}_y\rangle_y$",
               r"$\langle v\rangle_y$",
               r"$\langle \rho v\rangle_y$"]
    if normed:
        ylabels[1] = r"$\langle \mathbf{m}_x\rangle_y / \langle \rho \rangle_y$"
        ylabels[2] = r"$\langle \mathbf{m}_y\rangle_y / \langle \rho \rangle_y$"
    axes[0].legend()
    for i, ax in enumerate(axes):
        ax.set_ylabel(ylabels[i], fontsize="large")
    
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


def get_time_ave_density_profile(fin, beg_frame, end_frame=None, dt=0.01, inversed_A=False, inversed_B=True, shifting_mode="separated"):
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
        xm_A = get_center(x, rho_A, Lx, inversed=inversed_A)
        xm_B = get_center(x, rho_B, Lx, inversed=inversed_B)
        if shifting_mode == "A":
            xm_B = xm_A
        elif shifting_mode == "B":
            xm_A = xm_B
        rho_x_A[i] = shift_profile(rho_A, xm_A, Lx)
        rho_x_B[i] = shift_profile(rho_B, xm_B, Lx)
    rA = np.mean(rho_x_A, axis=0)
    rB = np.mean(rho_x_B, axis=0)
    return x, rA, rB


def get_time_ave_profile(fin, beg_frame, end_frame=None, dt=0.01, inversed_A=False, inversed_B=True, shifting_mode="separated"):
    Lx = get_Lx(fin)
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    if end_frame is None:
        end_frame = t.size
    n_frames = end_frame - beg_frame
    rho_x_A = np.zeros((n_frames, x.size))
    rho_x_B = np.zeros_like(rho_x_A)
    m_x_A = np.zeros_like(rho_x_A)
    m_x_B = np.zeros_like(rho_x_A)
    fields_x = np.zeros((n_frames, 10, x.size))
    for i, frame in enumerate(fields[beg_frame:end_frame]):
        fx = np.zeros((10, x.size))
        fx[:8] = np.mean(frame, axis=1)
        fx[8] = np.mean(frame[0]*frame[6], axis=0)
        fx[9] = np.mean(frame[1]*frame[7], axis=0)

        xm_A = get_center(x, fx[0], Lx, inversed=inversed_A)
        xm_B = get_center(x, fx[1], Lx, inversed=inversed_B)
        if shifting_mode == "A":
            xm_B = xm_A
        elif shifting_mode == "B":
            xm_A = xm_B
        for j in range(fx.shape[0]//2):
            jA = j*2
            jB = jA + 1
            fields_x[i][jA] = shift_profile(fx[jA], xm_A, Lx)
            fields_x[i][jB] = shift_profile(fx[jB], xm_B, Lx)

    fields_x = np.mean(fields_x, axis=0)
    return x, fields_x


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


def SP():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p5_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    f0 = f"{prefix}/L120_40_Dr0.100_k0.70_p0_25_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    # f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    Lx = get_Lx(f1)
    para = get_para(f1)
    dt = 0.01
    with np.load(f1, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    i_frame = -2
    plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True, para=para)

    x, fx = get_time_ave_profile(f1, beg_frame=50, end_frame=None, inversed_A=True, inversed_B=False)
    
    x0, fx0 = get_time_ave_profile(f0, beg_frame=50, end_frame=None, inversed_A=True, inversed_B=False)

    normed = True
    r1, r2 = None, None
    fig, axes = plt.subplots(5, 1, figsize=(6, 6), constrained_layout=True, sharex=True)
    if normed:
        fx[2] /= fx[0]
        fx[4] /= fx[0]
        fx[3] /= fx[1]
        fx[5] /= fx[1]

        fx0[3] /= fx0[1]
        fx0[5] /= fx0[1]
    
    for i, ax in enumerate(axes[:4]):
        ax.plot(x, fx[2*i], label="A")
        ax.plot(x, fx[2*i+1], label="B")
        ax.plot(x0, fx0[2*i+1], label="one species", c="tab:grey")

    axes[4].plot(x, fx[8])
    axes[4].plot(x, fx[9])
    axes[4].plot(x0, fx0[9], c="tab:grey")
    axes[0].set_xlim(0, Lx)

    axes[1].set_ylim(-0.05, 0.05)
    axes[2].set_ylim(-0.05, 0.05)


    axes[-1].set_xlabel(r"$x$", fontsize="large")
    ylabels = [r"$\langle \rho \rangle_{y,t}$",
               r"$\langle \mathbf{m}_x\rangle_{y,t}$",
               r"$\langle \mathbf{m}_y\rangle_{y,t}$",
               r"$\langle v\rangle_{y,t}$",
               r"$\langle \rho v\rangle_{y,t}$"]
    if normed:
        ylabels[1] = r"$\frac{\langle \mathbf{m}_x\rangle_{y,t}}{\langle \rho \rangle_{y,t}}$"
        ylabels[2] = r"$\frac{\langle \mathbf{m}_y\rangle_{y,t}}{\langle \rho \rangle_{y,t}}$"
    axes[0].legend()
    for i, ax in enumerate(axes):
        ax.set_ylabel(ylabels[i], fontsize="large")
    
    if r1 is not None:
        axes[0].axhline(r1[0], linestyle="dotted")
        axes[0].axhline(r1[1], linestyle="dotted")
    if r2 is not None:
        axes[0].axhline(r2[0], linestyle="dashed")
        axes[0].axhline(r2[1], linestyle="dashed")
    if para is None:
        title = r"$t=%g$" % t
    else:
        title = "$L_x=%g,L_y=%g,D_r=%g,\\eta=%g,\\eta_{AB}=-\\eta_{BA}=%g,\\rho_0=\\bar{\\rho}_{A,B}=%g,$\n$\\phi_A=%g,\\phi_B=%g$" % (
            para["Lx"], para["Ly"], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para['pA'], para['pB']
        )
    fig.suptitle(title)
    plt.show()
    plt.close()


def D1S1():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10/coarse_grain_dx1"
    f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_21010.npz"
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    # f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    Lx = get_Lx(f1)
    para = get_para(f1)
    dt = 0.01
    with np.load(f1, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    i_frame = -1
    plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True, para=para)

    x, fx = get_time_ave_profile(f1, beg_frame=50, end_frame=None, inversed_A=False, inversed_B=False)

    normed = True
    r1, r2 = None, None
    fig, axes = plt.subplots(5, 1, figsize=(6, 6), constrained_layout=True, sharex=True)
    if normed:
        fx[2] /= fx[0]
        fx[4] /= fx[0]
        fx[3] /= fx[1]
        fx[5] /= fx[1]
    
    for i, ax in enumerate(axes[:4]):
        ax.plot(x, fx[2*i], label="A")
        ax.plot(x, fx[2*i+1], label="B")

    axes[4].plot(x, fx[8])
    axes[4].plot(x, fx[9])
    axes[0].set_xlim(0, Lx)

    # axes[1].set_ylim(-0.05, 0.05)
    # axes[2].set_ylim(-0.05, 0.05)


    axes[-1].set_xlabel(r"$x$", fontsize="large")
    ylabels = [r"$\langle \rho \rangle_{y,t}$",
               r"$\langle \mathbf{m}_x\rangle_{y,t}$",
               r"$\langle \mathbf{m}_y\rangle_{y,t}$",
               r"$\langle v\rangle_{y,t}$",
               r"$\langle \rho v\rangle_{y,t}$"]
    if normed:
        ylabels[1] = r"$\langle \mathbf{m}_x\rangle_{y,t} / \langle \rho \rangle_{y,t}$"
        ylabels[2] = r"$\langle \mathbf{m}_y\rangle_{y,t} / \langle \rho \rangle_{y,t}$"
    axes[0].legend()
    for i, ax in enumerate(axes):
        ax.set_ylabel(ylabels[i], fontsize="large")
    
    if r1 is not None:
        axes[0].axhline(r1[0], linestyle="dotted")
        axes[0].axhline(r1[1], linestyle="dotted")
    if r2 is not None:
        axes[0].axhline(r2[0], linestyle="dashed")
        axes[0].axhline(r2[1], linestyle="dashed")
    if para is None:
        title = r"$t=%g$" % t
    else:
        title = "$L_x=%g,L_y=%g,D_r=%g,\\eta=%g,\\eta_{AB}=-\\eta_{BA}=%g,\\rho_0=\\bar{\\rho}_{A,B}=%g,$\n$\\phi_A=%g,\\phi_B=%g$" % (
            para["Lx"], para["Ly"], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para['pA'], para['pB']
        )
    fig.suptitle(title)
    plt.show()
    plt.close()


def D1S2():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10/coarse_grain_dx1"
    f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p15_10_r10_10_10_e-2.000_J0.500_-0.500_212510.npz"
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    # f1 = f"{prefix}/L120_40_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"
    Lx = get_Lx(f1)
    para = get_para(f1)
    dt = 0.01
    with np.load(f1, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    i_frame = -1
    plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True, para=para)

    x, fx = get_time_ave_profile(f1, beg_frame=60, end_frame=None, inversed_A=False, inversed_B=False, shifting_mode="B")

    normed = True
    r1, r2 = None, None
    fig, axes = plt.subplots(5, 1, figsize=(6, 6), constrained_layout=True, sharex=True)
    if normed:
        fx[2] /= fx[0]
        fx[4] /= fx[0]
        fx[3] /= fx[1]
        fx[5] /= fx[1]
    
    for i, ax in enumerate(axes[:4]):
        ax.plot(x, fx[2*i], label="A")
        ax.plot(x, fx[2*i+1], label="B")

    axes[4].plot(x, fx[8])
    axes[4].plot(x, fx[9])
    axes[0].set_xlim(0, Lx)

    # axes[1].set_ylim(-0.05, 0.05)
    # axes[2].set_ylim(-0.05, 0.05)


    axes[-1].set_xlabel(r"$x$", fontsize="large")
    ylabels = [r"$\langle \rho \rangle_{y,t}$",
               r"$\langle \mathbf{m}_x\rangle_{y,t}$",
               r"$\langle \mathbf{m}_y\rangle_{y,t}$",
               r"$\langle v\rangle_{y,t}$",
               r"$\langle \rho v\rangle_{y,t}$"]
    if normed:
        ylabels[1] = r"$\langle \mathbf{m}_x\rangle_{y,t} / \langle \rho \rangle_{y,t}$"
        ylabels[2] = r"$\langle \mathbf{m}_y\rangle_{y,t} / \langle \rho \rangle_{y,t}$"
    axes[0].legend()
    for i, ax in enumerate(axes):
        ax.set_ylabel(ylabels[i], fontsize="large")
    
    if r1 is not None:
        axes[0].axhline(r1[0], linestyle="dotted")
        axes[0].axhline(r1[1], linestyle="dotted")
    if r2 is not None:
        axes[0].axhline(r2[0], linestyle="dashed")
        axes[0].axhline(r2[1], linestyle="dashed")
    if para is None:
        title = r"$t=%g$" % t
    else:
        title = "$L_x=%g,L_y=%g,D_r=%g,\\eta=%g,\\eta_{AB}=-\\eta_{BA}=%g,\\rho_0=\\bar{\\rho}_{A,B}=%g,$\n$\\phi_A=%g,\\phi_B=%g$" % (
            para["Lx"], para["Ly"], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para['pA'], para['pB']
        )
    fig.suptitle(title)
    plt.show()
    plt.close()


def phase_transform():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10/coarse_grain_dx1"
    f1 = f"{prefix}/L480_160_Dr0.100_k0.70_p6.5_14.05_r10_10_10_e-2.000_J0.500_-0.500_2000.npz"
    Lx = get_Lx(f1)
    para = get_para(f1)
    dt = 0.01
    with np.load(f1, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]

    # r1=[3.88, 1.2]
    # r2=[10.54, 21.65]
    r1=None
    r2=None
    # for i_frame in range(156, t.size):
    #     print("frame %d/%d" % (i_frame, t.size))
    #     plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True, para=para, r1=r1, r2=r2)
    
    beg = 70
    # shifting_mode = "separated"
    shifting_mode = "B"
    x, rA_1, rB_1 = get_time_ave_density_profile(f1, beg, t.size, inversed_A=False, inversed_B=False, shifting_mode=shifting_mode)
    
    lineA, = plt.plot(x, rA_1)
    lineB, = plt.plot(x, rB_1)

    # x_edge_1 = [32, 456]
    # x_edge_2 = [200, 280]
    x_edge_1 = [35, 440]
    x_edge_2 = [210, 280]

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

    para = get_para(f1)
    phiA, phiB = para["pA"], para["pB"]
    
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



def level_rule():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10/coarse_grain_dx1"
    fins = [f"{prefix}/L480_160_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_21010.npz",
            f"{prefix}/L480_160_Dr0.100_k0.70_p8.625_13.36_r10_10_10_e-2.000_J0.500_-0.500_27510.npz"
    ]

    fig = plt.figure(figsize=(10, 4.2), constrained_layout=True)
    subfigs = fig.subfigures(1, 2, width_ratios=[4, 6])
    ax1 = subfigs[0].subplots(1, 1)
    ax1.plot(7.5, 10, "o")
    ax1.plot(8.625, 13.36, "o")
    ax1.plot(9.73, 16.77, "o")
    ax1.plot([3.8935, 10.5367], [1.2001, 21.6776], linestyle="dotted", c="tab:grey")
    ax1.plot(3.8935, 1.2001, "v", fillstyle="none", c="tab:grey", label="liquid binodal")
    ax1.plot(10.5367, 21.6776, "^", fillstyle="none", c="tab:brown", label="gas binodal")
    ax1.set_xlabel(r"$\phi_A$", fontsize="large")
    ax1.set_ylabel(r"$\phi_B$", fontsize="large")
    ax1.legend()

    (ax2, ax3) = subfigs[1].subplots(2, 1, sharex=True)
    begs = [70, 70]
    for i, fin in enumerate(fins):
        x, rA, rB = get_time_ave_density_profile(fin, begs[i], inversed_A=False, inversed_B=False)
        ax2.plot(x, rA)
        ax3.plot(x, rB)
    
    ax2.axhline(3.8935, linestyle="dotted", c="tab:grey")
    ax2.axhline(10.5367, linestyle="dashed", c="tab:brown")
    ax3.axhline(1.2001, linestyle="dotted", c="tab:grey")
    ax3.axhline(21.6776, linestyle="dashed", c="tab:brown")
    ax2.set_ylabel(r"$\langle \rho_A\rangle_{y,t}$", fontsize="large")
    ax3.set_ylabel(r"$\langle \rho_B\rangle_{y,t}$", fontsize="large")
    ax3.set_xlabel(r"$x$", fontsize="large")
    ax3.set_xlim(0, 480)

    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")

    title = r"$L_x=480,L_y=160,D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=0.5,\rho_0=\bar{\rho}_{A,B}=10$"
    fig.suptitle(title)

    plt.show()
    plt.close()


if __name__ == "__main__":
    # D1S2()
    SP()
    # plot_time_ave_density_profile_SP()
    # level_rule()
    # phase_transform()