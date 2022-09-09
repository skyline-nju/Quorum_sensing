import enum
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from add_line import add_line
from coarse_grain import coarse_grain_all


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_Lx(fname):
    basename = os.path.basename(fname)
    Lx = float(basename.lstrip("L").split("_")[0])
    return Lx


def get_para(fname):
    basename = os.path.basename(fname)
    s = basename.rstrip(".npz").split("_")
    para = {}
    if len(s) == 13:
        para["Lx"] = float(s[0].lstrip("L"))
        para["Ly"] = float(s[1])
        para["Dr"] = float(s[2].lstrip("Dr"))
        para["kappa"] = float(s[3].lstrip("k"))
        para["pA"] = float(s[4].lstrip("p"))
        para["pB"] = float(s[5])
        para["rho0"] = float(s[6].lstrip("r"))
        para["bar_rho_A"] = float(s[7])
        para["bar_rho_B"] = float(s[8])
        para["eta"] = float(s[9].lstrip("e"))
        para["eta_AB"] = float(s[10].lstrip("J"))
        para["eta_BA"] = float(s[11])
        para["seed"] = int(s[12])
    return para


def get_S_tilde(S, q_radial, q_module):
    """ Cal orientation-averaged structure factor"""
    Sq = np.zeros_like(q_radial)
    half_dq = (q_radial[1]-q_radial[0]) / 2

    for i, q in enumerate(q_radial):
        mask = np.logical_and(q_module >= q-half_dq, q_module < q+half_dq)
        Sq[i] = np.mean(S[mask])
    return Sq


def get_corr_length(q_radial, Sq):
    I0 = np.sum(Sq)
    I1 = np.sum(q_radial * Sq)
    if I0 > 0:
        corr_len = np.pi * 2 * I0 / I1
    else:
        corr_len = 0
    return corr_len


def get_R_t(fin, dx=1, dt=0.01):
    """ cal domain length Rt for A, B particles respectively"""

    folder, basename = os.path.split(fin)
    fcache = f"{folder}/cache/Rt/{basename}"
    if os.path.exists(fcache) and os.path.getmtime(fcache) > os.path.getmtime(fin):
        with np.load(fcache) as cache_data:
            t, R_t_A, R_t_B, R_t_tot = cache_data["t"], cache_data["R_t_A"], cache_data["R_t_B"], cache_data["R_t_tot"]
        print("loading data from", fcache)
    else:
        L = get_Lx(fin)
        with np.load(fin, "r") as data:
            t, x, y = data["t"] * dt, data["x"], data["y"]
            fields = data["fields"]
            
        n = int(L / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[1:qx.size//2-1]
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)
        qx_ij, qy_ij = np.meshgrid(qx, qy)
        q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        R_t_A = np.zeros(t.size)
        R_t_B = np.zeros(t.size)
        R_t_tot = np.zeros(t.size)
        for i, t_i in enumerate(t):
            rhoA_m, rhoB_m = np.mean(fields[i][0]), np.mean(fields[i][1])
            rhoA_q = np.fft.fft2(fields[i][0]-rhoA_m, norm="ortho")
            rhoB_q = np.fft.fft2(fields[i][1]-rhoB_m, norm="ortho")

            rhoA_q = np.fft.fftshift(rhoA_q)
            rhoB_q = np.fft.fftshift(rhoB_q)

            S_A = np.abs(rhoA_q) ** 2
            S_B = np.abs(rhoB_q) ** 2
            S_tot = np.abs(rhoA_q + rhoB_q) ** 2
            
            S_A_q = get_S_tilde(S_A, q_radial, q_module)
            S_B_q = get_S_tilde(S_B, q_radial, q_module)
            S_tot_q = get_S_tilde(S_tot, q_radial, q_module)

            R_t_A[i] = get_corr_length(q_radial, S_A_q)
            R_t_B[i] = get_corr_length(q_radial, S_B_q)
            R_t_tot[i] = get_corr_length(q_radial, S_tot_q)
        
        if not os.path.exists(f"{folder}/cache/Rt"):
            os.makedirs(f"{folder}/cache/Rt")
        np.savez_compressed(fcache, t=t, R_t_A=R_t_A, R_t_B=R_t_B, R_t_tot=R_t_tot)
        print("updating cache:", fcache)
    return t, R_t_A, R_t_B, R_t_tot


def plot_Rt(ax, fin, labels, which="both", ms=4):
    t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)
    if which == "both":
        line, = ax.plot(t, Rt_A, "o", fillstyle="none", ms=ms, label=labels[0])
        ax.plot(t, Rt_B, ">", c=line.get_c(), fillstyle="none", ms=ms, label=labels[1])
    elif which == "B":
        ax.plot(t, Rt_B, ">", fillstyle="none", ms=ms, label=labels)
    elif which == "A":
        ax.plot(t, Rt_A, "o", fillstyle="none", ms=ms, label=labels)
    else:
        print("Error, 'which' must be one of 'both', 'A' or 'B'")
        sys.exit()


def show_t(fin, dt=0.01):
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        for i in range(t.size):
            print(i, t[i])

def get_cb_extend(vmin, vmax):
    if vmin is None and vmax is None:
        extend = "neither"
    elif (vmin is None or vmin == 0) and vmax is not None:
        extend = "max"
    elif (vmin is not None and vmin != 0) and vmax is None:
        extend = "min"
    else:
        extend = "both"
    return extend


def plot_rho_A_B(fin, idx_frame, dx=1, dt=0.01, vmin=(None, None), vmax=(None, None)):
    Lx = get_Lx(fin)
    if not isinstance(vmin, tuple):
        vmin = (vmin, vmin)
    if not isinstance(vmax, tuple):
        vmax = (vmax, vmax)
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]

    if not isinstance(idx_frame, list):
        fig, axes = plt.subplots(1, 2, figsize=(6, 3.6), constrained_layout=True, sharey=True)
        fig.suptitle(r"$t=%g$" % t[idx_frame])
        im1 = axes[0].imshow(fields[idx_frame][0], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[0], vmax=vmax[0])
        im2 = axes[1].imshow(fields[idx_frame][1], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[1], vmax=vmax[1])
        cb1 = plt.colorbar(im1, ax=axes[0], orientation="horizontal", extend=get_cb_extend(vmin[0], vmax[0]))
        cb2 = plt.colorbar(im2, ax=axes[1], orientation="horizontal", extend=get_cb_extend(vmin[1], vmax[1]))
        cb1.set_label(r"$\rho_A$")
        cb2.set_label(r"$\rho_B$")
    else:
        figsize = (2.5 * len(idx_frame)+0.4, 5)
        fig, axes = plt.subplots(
            2, len(idx_frame), figsize=figsize, constrained_layout=True, sharex=True, sharey=True)
        for i, i_frame in enumerate(idx_frame):
            im1 = axes[0, i].imshow(fields[i_frame][0], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[0], vmax=vmax[0])
            im2 = axes[1, i].imshow(fields[i_frame][1], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[1], vmax=vmax[1])

            axes[0, i].set_title(r"$t=%g$" % t[i_frame])
            if i == len(idx_frame) - 1:
                cb1 = plt.colorbar(im1, ax=axes[0], extend=get_cb_extend(vmin[0], vmax[0]))
                cb2 = plt.colorbar(im2, ax=axes[1], extend=get_cb_extend(vmin[1], vmax[1]))
                cb1.set_label(r"$\rho_A$", fontsize="large")
                cb2.set_label(r"$\rho_B$", fontsize="large")
    plt.show()
    plt.close()


def plot_rhoA_B_Sq(fin, i_frame, dx=1, dt=0.01, vmin=(None, None), vmax=(None, None)):
    def plot_one_frame(i, axes):
        im1 = axes[0].imshow(fields[i][0], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[0], vmax=vmax[0])
        im2 = axes[1].imshow(fields[i][1], origin="lower", extent=[0, Lx, 0, Lx], vmin=vmin[1], vmax=vmax[1])

        rhoA_q = np.fft.fft2(fields[i][0]-np.mean(fields[i][0]), norm="ortho")
        rhoB_q = np.fft.fft2(fields[i][1]-np.mean(fields[i][1]), norm="ortho")

        rhoA_q = np.fft.fftshift(rhoA_q)
        rhoB_q = np.fft.fftshift(rhoB_q)

        S_A = np.abs(rhoA_q) ** 2
        S_B = np.abs(rhoB_q) ** 2
        
        S_A_q = get_S_tilde(S_A, q_radial, q_module)
        S_B_q = get_S_tilde(S_B, q_radial, q_module)

        axes[2].plot(q_radial, S_A_q, label=r"$\tilde{S}_A$")
        axes[2].plot(q_radial, S_B_q, label=r"$\tilde{S}_B$")
        axes[2].set_xscale("log")
        axes[2].set_yscale("log")
        axes[2].legend()
        return im1, im2

    Lx = get_Lx(fin)
    if not isinstance(vmin, tuple):
        vmin = (vmin, vmin)
    if not isinstance(vmax, tuple):
        vmax = (vmax, vmax)
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    n = int(Lx / dx)
    qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
    qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
    q_radial = qx[1:qx.size//2-1]
    qx = np.fft.fftshift(qx)
    qy = np.fft.fftshift(qy)
    qx_ij, qy_ij = np.meshgrid(qx, qy)
    q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)

    if not isinstance(i_frame, list):
        fig, axes = plt.subplots(1, 3, figsize=(9, 3), constrained_layout=True)
        plot_one_frame(i_frame, axes)
    else:
        figsize = (2.6 * len(i_frame)+0.2, 8)
        fig, axes = plt.subplots(3, len(i_frame), figsize=figsize, constrained_layout=True, sharey='row')
        for i in range(len(i_frame)):
            im1, im2 = plot_one_frame(i_frame[i], axes[:, i])
            axes[0, i].set_title(r"$t=%g$" % t[i_frame[i]])
            axes[0, i].set_xticklabels([])
            axes[2, i].set_xlabel(r"$q$", fontsize="large")
            if i == len(i_frame) - 1:
                cb1 = plt.colorbar(im1, ax=axes[0], extend=get_cb_extend(vmin[0], vmax[0]))
                cb2 = plt.colorbar(im2, ax=axes[1], extend=get_cb_extend(vmin[1], vmax[1]))
                # cb1.set_label(r"$\rho_A$", fontsize="large")
                # cb2.set_label(r"$\rho_B$", fontsize="large")
        axes[0, 0].set_ylabel(r"$\rho_A(\mathbf{r}, t)$", fontsize="large")
        axes[1, 0].set_ylabel(r"$\rho_B(\mathbf{r}, t)$", fontsize="large")
        axes[2, 0].set_ylabel(r"$\tilde{S}(q, t)$", fontsize="large")
    
    para = get_para(fin)        
    if para["rho0"] == para["bar_rho_A"] == para["bar_rho_B"]:
        title = r"$L=%g,D_r=%g,\eta=%g,\eta_{AB}=-\eta_{BA}=%g,\rho_0=\bar{\rho}_{A,B}=%g,\phi_A=%g,\phi_B=%g$" % (
            para["Lx"], para["Dr"], para["eta"], para["eta_AB"], para["rho0"], para["pA"], para["pB"]
        )
    fig.suptitle(title, fontsize="x-large")
    plt.show()
    plt.close()

def plot_Rt_L320():
    f0 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p0_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    f1 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    # f1 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p90_30_r40_e-2.000_J0.500_-0.500_1000.npz"


    t, Rt_A, Rt_B, Rt_tot = get_R_t(f1)
    t0, Rt_A0, Rt_B0, Rt_tot0 = get_R_t(f0)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    line, = ax.plot(t, Rt_A, "o", label=r"$R_A: \phi_A=40,\phi_B=40$", fillstyle="none", ms=5)
    ax.plot(t, Rt_B, "s", label=r"$R_B: \phi_A=40,\phi_B=40$", c=line.get_c(), fillstyle="none", ms=5)
    # ax.plot(t, Rt_tot, "o", label=r"$R_{\rm tot}: \phi_A=40,\phi_B=40$")

    ax.plot(t0, Rt_B0, "s", fillstyle="none", c="tab:grey", label=r"$R_B: \phi_A=0,\phi_B=40$", ms=5)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$R$", fontsize="large")
    ax.legend(fontsize="large")

    title = r"$L=320, \eta=-2, \eta_{\rm AB}=-\eta_{\rm BA}=0.5, \rho_0=\bar{\rho}_A=\bar{\rho}_B=40$"
    ax.set_title(title, fontsize="x-large")

    add_line(ax, 0.25, 0.22, 0.95, 1/3, label="0.33", xl=0.5, yl=0.55, fontsize="large")
    # add_line(ax, 0.2, 0.08, 0.9, 1/3, label="0.33", fontsize="large")

    add_line(ax, 0.33, 0.12, 0.5, 0.18, label="0.18", xl=0.4, yl=0.1, c="tab:pink", fontsize="large")
    add_line(ax, 0.64, 0.29, 0.75, 0.6, label="0.6", xl=0.7, yl=0.3, c="tab:green", fontsize="large")
    # add_line(ax, 0.2, 0.05, 0.6, 0.22, label="0.22", xl=0.5, yl=0.3, c="k", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_bubble_burst():
    prefix1 = f"{root_rudabeh}/QS_MPI_noself/L320_r20_20_20"
    prefix2 = f"{root_rudabeh}/QS_MPI_noself/L160_r20_20_20"
    coarse_grain_all(prefix1, dx=1, pat="*_e-2.000_J2.000_-2.000_*.gsd")
    coarse_grain_all(prefix2, dx=1, pat="*_e-2.000_J2.000_-2.000_*.gsd")
    folder1 = f"{prefix1}/coarse_grain_dx1"
    folder2 = f"{prefix2}/coarse_grain_dx1"

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True)
    ms = 3

    fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p5_32_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder2}/L160_160_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
    ]
    labels = [[r"$R_A:\phi_B=32$", r"$R_B:\phi_B=32$"],
              [r"$R_A:\phi_B=33$", r"$R_B:\phi_B=33$"]]
    for i, fin in enumerate(fins):
        plot_Rt(ax1, fin, labels[i], ms=ms)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$R$", fontsize="large")

    fins = [f"{folder1}/L320_320_Dr0.100_k0.70_p5_32_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p5_34_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p5_35_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
        ]
    labels = [[r"$R_A:\phi_B=32$", r"$R_B:\phi_B=32$"],
              [r"$R_A:\phi_B=33$", r"$R_B:\phi_B=33$"],
              [r"$R_A:\phi_B=34$", r"$R_B:\phi_B=34$"],
              [r"$R_A:\phi_B=35$", r"$R_B:\phi_B=35$"]]
    for i, fin in enumerate(fins):
        plot_Rt(ax2, fin, labels[i], ms=ms)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$t$", fontsize="large")

    fins = [f"{folder1}/L320_320_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder2}/L160_160_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"]
    labels = [[r"$R_A: L=320$", r"$R_B: L=320$"],
              [r"$R_A: L=160$", r"$R_B: L=160$"]]

    for i, fin in enumerate(fins):
        plot_Rt(ax3, fin, labels[i], ms=ms)

    # ax2.set_xscale("log")
    ax3.set_xlabel(r"$t$", fontsize="large")

    titles = [r"(a) $L=160,\phi_A=5$", r"(b) $L=320,\phi_A=5$", r"(c) $\phi_A=5,\phi_B=33$"]
    for i, ax in enumerate([ax1, ax2, ax3]):
        ax.set_title(titles[i])
        ax.legend()

    title = r"$D_r=0.1, \eta=-2,\eta_{AB}=\eta_{BA}=2,\rho_0=\bar{\rho}_A=\bar{\rho}_B=20$"
    fig.suptitle(title, fontsize="x-large")

    add_line(ax2, 0.3, 0.25, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")
    # add_line(ax, 0.3, 0.25, 0.95, 0.3, label="0.33", xl=0.5, yl=0.5, fontsize="large")

    add_line(ax2, 0.2, 0.0, 0.6, 0.2, label="0.2", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()

    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    prefix1 = f"{root_rudabeh}/QS_MPI_noself/L160_r40_40_40"
    coarse_grain_all(prefix1, dx=1, pat="*_p10_*_r40_40_40_e-2.000_J2.000_-2.000_*.gsd")
    folder1 = f"{prefix1}/coarse_grain_dx1"

    fins = [f"{folder1}/L160_160_Dr0.100_k0.70_p10_68_r40_40_40_e-2.000_J2.000_-2.000_21070.npz",
            f"{folder1}/L160_160_Dr0.100_k0.70_p10_70_r40_40_40_e-2.000_J2.000_-2.000_1000.npz",
        ]
    labels = [[r"$R_A:\phi_B=68$", r"$R_B:\phi_B=68$"],
              [r"$R_A:\phi_B=70$", r"$R_B:\phi_B=70$"]]
    for i, fin in enumerate(fins):
        plot_Rt(ax, fin, labels[i], ms=ms)
    ax.legend()
    
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$R$")

    title = r"$L=160,D_r=0.1,\eta=-2,\eta_{\rm AB}=-\eta_{\rm BA}=2,\rho_0=\bar{\rho}_{A,B}=40,\phi_A=10$"
    ax.set_title(title)

    plt.show()
    plt.close()


def plot_Rt_SP_osi_band():
    ms=3.0

    prefix1 = f"{root_rudabeh}/QS_MPI_noself/L320_r10_10_10"
    coarse_grain_all(prefix1, dx=1, pat="*_r10_10_10_e-2.000_J0.500_-0.500_*.gsd")
    folder1 = f"{prefix1}/coarse_grain_dx1"

    prefix2 = f"{root_rudabeh}/QS_MPI_noself/L160_r10_10_10"
    coarse_grain_all(prefix2, dx=1, pat="*_r10_10_10_e-2.000_J0.500_-0.500_*.gsd")
    folder2 = f"{prefix2}/coarse_grain_dx1"

    prefix3 = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10"
    coarse_grain_all(prefix3, dx=1, pat="*_r10_10_10_e-2.000_J0.500_-0.500_*.gsd")
    folder3 = f"{prefix3}/coarse_grain_dx1"

    prefix4 = f"{root_rudabeh}/QS_MPI_noself/L160_r20_20_20"
    coarse_grain_all(prefix4, dx=1, pat="*_r20_20_20_e-2.000_J0.500_-0.500_*.gsd")
    folder4 = f"{prefix4}/coarse_grain_dx1"

    fig, axes = plt.subplots(2, 3, figsize=(10, 8), constrained_layout=True, sharey=True)

    fins = [f"{folder1}/L320_320_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p10_40_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
        ]
    labels = [[r"$R_A:\phi_B=28$", r"$R_B:\phi_B=28$"],
              [r"$R_A:\phi_B=34$", r"$R_B:\phi_B=34$"],
              [r"$R_A:\phi_B=40$", r"$R_B:\phi_B=40$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 0], fin, labels[i], ms=ms)
        plot_Rt(axes[1, 0], fin, labels[i], ms=ms)


    fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder3}/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    # fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
    #         f"{folder2}/L160_160_Dr0.100_k0.70_p20_56_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=160$", r"$R_B:L=160$"],
              [r"$R_A:L=320$", r"$R_B:L=320$"],
              [r"$R_A:L=640$", r"$R_B:L=640$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 1], fin, labels[i], ms=ms)
        plot_Rt(axes[1, 1], fin, labels[i], ms=ms)


    fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder1}/L320_320_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder3}/L640_640_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    # fins = [f"{prefix}/L160_160_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
    #         f"{prefix}/L160_160_Dr0.100_k0.70_p20_68_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=160$", r"$R_B:L=160$"],
              [r"$R_A:L=320$", r"$R_B:L=320$"],
              [r"$R_A:L=640$", r"$R_B:L=640$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 2], fin, labels[i], ms=ms)
        plot_Rt(axes[1, 2], fin, labels[i], ms=ms)
    
    for i in range(3):
        axes[0, i].set_xscale("log")
        axes[0, i].set_yscale("log")
        axes[1, i].set_yscale("log")
        axes[0, i].set_xlabel(r"$t$", fontsize="large")
        axes[1, i].set_xlabel(r"$t$", fontsize="large")
        axes[0, i].legend()
        axes[1, i].legend()

    
    titles = [r"(a) $L=320$", r"(b) $\phi_A=10,\phi_B=28$", r"(c) $\phi_A=10,\phi_B=34$",
              r"(d) $L=320$", r"(e) $\phi_A=10,\phi_B=28$", r"(f) $\phi_A=10,\phi_B=34$"]
    for i, ax in enumerate(axes.flat):
        ax.set_title(titles[i])
    axes[0, 0].set_ylabel(r"$R$", fontsize="large")
    axes[1, 0].set_ylabel(r"$R$", fontsize="large")

    title = r"$D_r=0.1, \eta=-2,\eta_{AB}=-\eta_{BA}=0.5,\rho_0=\bar{\rho}_A=\bar{\rho}_B=10$"
    fig.suptitle(title, fontsize="x-large")

    add_line(axes[0, 0], 0.3, 0.3, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")
    add_line(axes[0, 1], 0.3, 0.3, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")
    add_line(axes[0, 2], 0.3, 0.3, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")


    # add_line(ax, 0.3, 0.25, 0.95, 0.3, label="0.33", xl=0.5, yl=0.5, fontsize="large")

    # add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()

    fig, axes = plt.subplots(2, 2, figsize=(6, 7), constrained_layout=True, sharey=True)

    fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder4}/L160_160_Dr0.100_k0.70_p20_56_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A: \rho_0=10$", r"$R_B: \rho_0=10$"],
              [r"$R_A: \rho_0=20$", r"$R_B: \rho_0=20$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 0], fin, labels[i], ms=ms)
        plot_Rt(axes[1, 0], fin, labels[i], ms=ms)

    fins = [f"{folder2}/L160_160_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder4}/L160_160_Dr0.100_k0.70_p20_68_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A: \rho_0=10$", r"$R_B: \rho_0=10$"],
              [r"$R_A: \rho_0=20$", r"$R_B: \rho_0=20$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 1], fin, labels[i], ms=ms)
        plot_Rt(axes[1, 1], fin, labels[i], ms=ms)

    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_xscale("log")

    axes[0, 1].set_yscale("log")
    titles = [r"(a) $\phi_A/\rho_0=1,\phi_B/\rho_0=2.8$",
              r"(b) $\phi_A/\rho_0=1,\phi_B/\rho_0=3.4$",
              r"(c) $\phi_A/\rho_0=1,\phi_B/\rho_0=2.8$",
              r"(d) $\phi_A/\rho_0=1,\phi_B/\rho_0=3.4$"]
    for i, ax in enumerate(axes.flat):
        ax.legend()
        ax.set_title(titles[i])
        ax.set_xlabel(r"$t$")
    axes[0, 0].set_ylabel(r"$R$")
    axes[1, 0].set_ylabel(r"$R$")
    
    title = r"$L=160, D_r=0.1, \eta=-2,\eta_{AB}=-\eta_{BA}=0.5,\rho_0=\bar{\rho}_A=\bar{\rho}_B$"
    fig.suptitle(title, fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_L320_r16():
    fins = ["/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p36_8_r16_16_16_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p36_9_r16_16_16_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p36_10_r16_16_16_e-2.000_J0.500_-0.500_1000.npz"
        ]

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    labels = [[r"$R_A:\phi_A=36,\phi_B=8$", r"$R_B:\phi_A=36,\phi_B=8$"],
              [r"$R_A:\phi_A=36,\phi_B=9$", r"$R_B:\phi_A=36,\phi_B=9$"],
              [r"$R_A:\phi_A=36,\phi_B=10$", r"$R_B:\phi_A=36,\phi_B=10$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)

        line, = ax.plot(t, Rt_A, "o", fillstyle="none", ms=5, label=labels[i][0])
        ax.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=5, label=labels[i][1])
        # ax.plot(t, Rt_tot, ">", c=line.get_c(), fillstyle="none")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$R$", fontsize="large")
    ax.legend(fontsize="large")

    title = r"$L=320, \eta=-2, \eta_{\rm AB}=-\eta_{\rm BA}=2, \rho_0=\bar{\rho}_A=\bar{\rho}_B=20$"
    ax.set_title(title, fontsize="x-large")

    add_line(ax, 0.3, 0.25, 0.95, 1/3, label="0.33", xl=0.5, yl=0.5, fontsize="large")
    # add_line(ax, 0.3, 0.25, 0.95, 0.3, label="0.33", xl=0.5, yl=0.5, fontsize="large")

    # add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_SP_spiral_cluster():
    prefix1 = f"{root_rudabeh}/QS_MPI_noself/L160_r20_20_20"
    coarse_grain_all(prefix1, dx=1, pat="*_J2.000_-2.000_*.gsd")
    folder1 = f"{prefix1}/coarse_grain_dx1"

    prefix2 = f"{root_rudabeh}/QS_MPI_noself/L320_r20_20_20"
    coarse_grain_all(prefix2, dx=1, pat="*._J2.000_-2.000_*.gsd")
    folder2 = f"{prefix2}/coarse_grain_dx1"
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True, sharey=True)

    fins = [f"{folder1}/L160_160_Dr0.100_k0.70_p15_3_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder1}/L160_160_Dr0.100_k0.70_p15_4_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder1}/L160_160_Dr0.100_k0.70_p15_5_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"]
    labels = [[r"$R_A:\phi_B=3$", r"$R_B:\phi_B=3$"],
              [r"$R_A:\phi_B=4$", r"$R_B:\phi_B=4$"],
              [r"$R_A:\phi_B=5$", r"$R_B:\phi_B=5$"]]

    ms = 4
    for i, fin in enumerate(fins):
        plot_Rt(ax1, fin, labels[i], ms=ms)
        plot_Rt(ax2, fin, labels[i], ms=ms)

    ax1.set_xlabel(r"$t$")
    ax2.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$R$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    ax1.set_title(r"(a) $L=160,\phi_A=15$")
    ax2.set_title(r"(b) $L=160,\phi_A=15$")
    ax1.legend()
    # ax2.legend(title=r"$L=160,\phi_A=15$")

    title = r"$D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=2, \rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title, fontsize="x-large")

    plt.show()
    plt.close()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True, sharey=True)
    fins = [f"{folder1}/L160_160_Dr0.100_k0.70_p15_5_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            f"{folder2}/L320_320_Dr0.100_k0.70_p15_5_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"]
    labels = [[r"$R_A:L=160$", r"$R_B:L=160$"],
              [r"$R_A:L=320$", r"$R_B:L=320$"]]

    ms = 4
    for i, fin in enumerate(fins):
        plot_Rt(ax1, fin, labels[i], ms=ms)
        plot_Rt(ax2, fin, labels[i], ms=ms)

    ax1.set_xlabel(r"$t$")
    ax2.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$R$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    ax1.set_title(r"(a) $\phi_A=15,\phi_B=5$")
    ax2.set_title(r"(b) $\phi_A=15,\phi_B=5$")
    ax1.legend()
    # ax2.legend(title=r"$L=160,\phi_A=15$")

    title = r"$D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=2, \rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title, fontsize="x-large")

    plt.show()
    plt.close()


def plot_Rt_L160():
    # fins = ["/home/yduan/data/tmp/L160_160_Dr0.100_k0.70_p10_68_r40_40_40_e-2.000_J2.000_-2.000_21070.npz",
    #         "/home/yduan/data/tmp/L160_160_Dr0.100_k0.70_p10_70_r40_40_40_e-2.000_J2.000_-2.000_1000.npz"
    #     ]
    fins = ["/home/yduan/data/tmp/L160_160_Dr0.100_k0.70_p10_70_r40_40_40_e-2.000_J2.000_-2.000_1000.npz"]
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True, sharey=True)

    labels = [[r"$R_A:\phi_A=10,\phi_B=70$", r"$R_B:\phi_A=10,\phi_B=70$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)
        for ax in axes:
            line, = ax.plot(t, Rt_A, "o", fillstyle="none", ms=4, label=labels[i][0])
            ax.plot(t, Rt_B, "s", fillstyle="none", ms=4, label=labels[i][1])

    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")


    axes[0].set_xlabel(r"$t$", fontsize="large")
    axes[1].set_xlabel(r"$t$", fontsize="large")
    axes[0].set_ylabel(r"$R$", fontsize="large")
    axes[1].set_ylabel(r"$R$", fontsize="large")
    axes[0].legend(fontsize="large")

    title = r"$L=160, \eta=-2, \eta_{\rm AB}=-\eta_{\rm BA}=2, \rho_0=\bar{\rho}_A=\bar{\rho}_B=40$"
    fig.suptitle(title, fontsize="x-large")

    add_line(axes[0], 0.01, 0.1, 0.5, 1/3, label="0.33", xl=0.15, yl=0.5, fontsize="large")
    # add_line(ax, 0.3, 0.25, 0.95, 0.3, label="0.33", xl=0.5, yl=0.5, fontsize="large")

    # add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_varied_composition():
    prefix1 = f"{root_tahmineh}/QS_MPI_noself/L640_r20_20_20"
    coarse_grain_all(prefix1, dx=1, pat="*_e-2.000_J0.500_-0.500_*.gsd")
    folder1 = f"{prefix1}/coarse_grain_dx1"

    fig, axes = plt.subplots(2, 2, figsize=(8, 7), sharex=True, sharey=True, constrained_layout=True)

    f0 = f"{folder1}/L640_640_Dr0.100_k0.70_p0_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    plot_Rt(axes[0,0], f0, labels=r"$R_B:\phi_A=0,\phi_B=20$", which="B")

    f1 = [f"{folder1}/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
          f"{folder1}/L640_640_Dr0.100_k0.70_p12_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:\phi_A=20,\phi_B=20$", r"$R_B:\phi_A=20,\phi_B=20$"],
              [r"$R_A:\phi_A=12,\phi_B=20$", r"$R_B:\phi_A=12,\phi_B=20$"]]
    for i, fin in enumerate(f1):
        plot_Rt(axes[0, 1], fin, labels[i])

    f2 = [f"{folder1}/L640_640_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
          f"{folder1}/L640_640_Dr0.100_k0.70_p14_10_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]

    labels = [[r"$R_A:\phi_A=12,\phi_B=6$", r"$R_B:\phi_A=12,\phi_B=6$"],
              [r"$R_A:\phi_A=14,\phi_B=10$", r"$R_B:\phi_A=14,\phi_B=10$"]]
    for i, fin in enumerate(f2):
        plot_Rt(axes[1, 0], fin, labels[i], ms=4)
    
    prefix2 = f"{root_rudabeh}/QS_MPI_noself/L320_r20_20_20"
    coarse_grain_all(prefix2, dx=1, pat="*_e-2.000_J0.500_-0.500_*.gsd")
    folder2 = f"{prefix2}/coarse_grain_dx1"

    f3 = f"{folder2}/L320_320_Dr0.100_k0.70_p40_40_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    labels = [r"$R_A:\phi_A=40,\phi_B=40$", r"$R_B:\phi_A=40,\phi_B=40$"]
    plot_Rt(axes[1, 1], f3, labels=labels)

    titles = ["(a)", "(b)", "(c)", "(d)"]
    for i, ax in enumerate(axes.flat):
        ax.set_xscale("log")
        ax.set_yscale("log")
        if i == 3:
            ax.legend(title=r"$L=320$")
        else:
            ax.legend(title=r"$L=640$")
        ax.set_title(titles[i])
    
    axes[0, 0].set_ylabel(r"$R(t)$")
    axes[1, 0].set_ylabel(r"$R(t)$")
    axes[1, 0].set_xlabel(r"$t$")
    axes[1, 1].set_xlabel(r"$t$")

    title = r"$\eta=-2,\eta_{AB}=-\eta_{BA}=0.5, D_r=0.1, \rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title)

    add_line(axes[0, 0], 0.2, 0.25, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")
    add_line(axes[0, 1], 0.2, 0.25, 0.95, 1/3, label="0.33", xl=0.2, yl=0.4, fontsize="large")
    add_line(axes[1, 0], 0.2, 0.25, 0.95, 1/3, label="0.33", xl=0.2, yl=0.4, fontsize="large")
    add_line(axes[1, 0], 0.55, 0.2, 0.67, 1, label="1", xl=0.65, yl=0.35, fontsize="large")
    add_line(axes[1, 0], 0.7, 0.53, 0.95, 0.2, label="0.2", xl=0.93, yl=0.62, fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_one_5th():    
    pat = "*_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.gsd"
    prefix1 = f"{root_rudabeh}/QS_MPI_noself/L320_r20_20_20"
    coarse_grain_all(prefix1, dx=1, pat=pat)
    folder1 = f"{prefix1}/coarse_grain_dx1"

    prefix2 = f"{root_tahmineh}/QS_MPI_noself/L640_r20_20_20"
    coarse_grain_all(prefix2, dx=1, pat=pat)
    folder2 = f"{prefix2}/coarse_grain_dx1"

    prefix3 = f"{root_rudabeh}/QS_MPI_noself/L1280_r20_20_20"
    coarse_grain_all(prefix3, dx=1, pat=pat)
    folder3 = f"{prefix3}/coarse_grain_dx1"

    prefix4 = f"{root_rudabeh}/QS_MPI_noself/L320_r10_10_10"
    coarse_grain_all(prefix4, dx=1, pat="*_r10_10_10_e-2.000_J0.500_-0.500_1000.gsd")
    folder4 = f"{prefix4}/coarse_grain_dx1"

    prefix5 = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10"
    coarse_grain_all(prefix5, dx=1, pat="*_r10_10_10_e-2.000_J0.500_-0.500_1000.gsd")
    folder5 = f"{prefix5}/coarse_grain_dx1"

    prefix6 = f"{root_tahmineh}/QS_MPI_noself/L640_r5_5_5"
    coarse_grain_all(prefix6, dx=1, pat="*_r5_5_5_e-2.000_J0.500_-0.500_1000.gsd")
    folder6 = f"{prefix6}/coarse_grain_dx1"

    prefix7 = f"{root_rudabeh}/QS_MPI_noself/L1280_r10_10_10"
    coarse_grain_all(prefix7, dx=1)
    folder7 = f"{prefix7}/coarse_grain_dx1"

    prefix8 = f"{root_rudabeh}/QS_MPI_noself/L1280_r5_5_5"
    coarse_grain_all(prefix8, dx=1)
    folder8 = f"{prefix8}/coarse_grain_dx1"

    prefix9 = f"{root_rudabeh}/QS_MPI_noself/L320_r5_5_5"
    coarse_grain_all(prefix9, dx=1)
    folder9 = f"{prefix9}/coarse_grain_dx1"

    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    fins = [f"{folder1}/L320_320_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder2}/L640_640_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder3}/L1280_1280_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=320$", r"$R_B:L=320$"],
              [r"$R_A:L=640$", r"$R_B:L=640$"],
              [r"$R_A:L=1280$", r"$R_B:L=1280$"]]

    for i, fin in enumerate(fins):
        plot_Rt(ax, fin, labels[i], ms=3.5)
    
    ax.legend(title=r"$\phi_A=12,\phi_B=6$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$R$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    title = r"$\eta=-2,\eta_{AB}=-\eta_{BA}=0.5, D_r=0.1, \rho_0=\bar{\rho}_{A,B}=20$"
    fig.suptitle(title)

    add_line(ax, 0.75, 0.65, 0.99, 0.2, label="0.2", xl=0.93, yl=0.7, fontsize="large")
    add_line(ax, 0.45, 0.35, 0.65, 1.0, label="1", xl=0.5, yl=0.7, fontsize="large")
    # add_line(ax, 0.4, 0.05, 0.65, 0.5, label="0.5", xl=0.5, yl=0.1, fontsize="large")

    plt.show()
    plt.close()

    fig, axes = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True, sharex=True, sharey=True)

    fins = [f"{folder1}/L320_320_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder4}/L320_320_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:\rho_0=20$", r"$R_B:\rho_0=20$"],
              [r"$R_A:\rho_0=10$", r"$R_B:\rho_0=10$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 0], fin, labels[i], ms=3.5)
    axes[0, 0].set_ylabel(r"$R$")

    fins = [f"{folder2}/L640_640_Dr0.100_k0.70_p12_6_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder5}/L640_640_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder6}/L640_640_Dr0.100_k0.70_p3_1.5_r5_5_5_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:\rho_0=20$", r"$R_B:\rho_0=20$"],
              [r"$R_A:\rho_0=10$", r"$R_B:\rho_0=10$"],
              [r"$R_A:\rho_0=5$", r"$R_B:\rho_0=5$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[0, 1], fin, labels[i], ms=3.5)

    fins = [f"{folder4}/L320_320_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder5}/L640_640_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder7}/L1280_1280_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=320$", r"$R_B:L=320$"],
              [r"$R_A:L=640$", r"$R_B:L=640$"],
              [r"$R_A:L=1280$", r"$R_B:L=1280$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[1, 0], fin, labels[i], ms=3.5)
    axes[1, 0].set_xlabel(r"$t$")
    axes[1, 0].set_ylabel(r"$R$")

    fins = [f"{folder9}/L320_320_Dr0.100_k0.70_p3.25_2.25_r5_5_5_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder6}/L640_640_Dr0.100_k0.70_p3.25_2.25_r5_5_5_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder8}/L1280_1280_Dr0.100_k0.70_p3.25_2.25_r5_5_5_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=320$", r"$R_B:L=320$"],
              [r"$R_A:L=640$", r"$R_B:L=640$"],
              [r"$R_A:L=1280$", r"$R_B:L=1280$"]]
    for i, fin in enumerate(fins):
        plot_Rt(axes[1, 1], fin, labels[i], ms=3.5)
    axes[1, 1].set_xlabel(r"$t$")

    titles = [r"(a) $L=320, \phi_A/\rho_0=0.6,\phi_B/\rho_0=0.3$",
              r"(b) $L=640, \phi_A/\rho_0=0.6,\phi_B/\rho_0=0.3$",
              r"(c) $\phi_A=6,\phi_B=3,\rho_0=10$",
              r"(d) $\phi_A=3.25,\phi_B=2.25,\rho_0=5$"]
    for i, ax in enumerate(axes.flat):
        ax.set_title(titles[i])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()

    title = r"$\eta=-2,\eta_{AB}=-\eta_{BA}=0.5, D_r=0.1, \rho_0=\bar{\rho}_{A,B}$"
    fig.suptitle(title)

    add_line(axes[1, 0], 0.75, 0.65, 0.99, 0.2, label="0.2", xl=0.93, yl=0.7, fontsize="large")
    add_line(axes[0, 0], 0.75, 0.65, 0.99, 0.2, label="0.2", xl=0.93, yl=0.7, fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_one_3rd():
    prefix0 = f"{root_tahmineh}/QS_MPI_noself/L640_r5_5_5"
    coarse_grain_all(prefix0, dx=1)
    folder0 = f"{prefix0}/coarse_grain_dx1"

    prefix1 = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10"
    coarse_grain_all(prefix1, dx=1)
    folder1 = f"{prefix1}/coarse_grain_dx1"

    prefix2 = f"{root_tahmineh}/QS_MPI_noself/L640_r20_20_20"
    coarse_grain_all(prefix2, dx=1)
    folder2 = f"{prefix2}/coarse_grain_dx1"

    prefix3 = f"{root_rudabeh}/QS_MPI_noself/L1280_r5_5_5"
    coarse_grain_all(prefix3, dx=1)
    folder3 = f"{prefix3}/coarse_grain_dx1"


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5), sharex=True, sharey=True, constrained_layout=True)
    fins = [f"{folder0}/L640_640_Dr0.100_k0.70_p5_5_r5_5_5_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder1}/L640_640_Dr0.100_k0.70_p10_10_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder2}/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:\phi_A=\phi_B=\rho_0=5$", r"$R_B:\phi_A=\phi_B=\rho_0=5$"],
              [r"$R_A:\phi_A=\phi_B=\rho_0=10$", r"$R_B:\phi_A=\phi_B=\rho_0=10$"],
              [r"$R_A:\phi_A=\phi_B=\rho_0=20$", r"$R_B:\phi_A=\phi_B=\rho_0=20$"]]

    for i, fin in enumerate(fins):
        plot_Rt(ax1, fin, labels[i], ms=3.5)
    
    ax1.legend(title=r"$L=640$")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$R$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    fins = [f"{folder0}/L640_640_Dr0.100_k0.70_p5_5_r5_5_5_e-2.000_J0.500_-0.500_1000.npz",
            f"{folder3}/L1280_1280_Dr0.100_k0.70_p5_5_r5_5_5_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=640$", r"$R_B:L=640$"],
              [r"$R_A:L=1280$", r"$R_B:L=1280$"]]
    for i, fin in enumerate(fins):
        plot_Rt(ax2, fin, labels[i], ms=3.5)
    ax2.legend(title=r"$\phi_A=\phi_B=\rho_0=5$")
    ax2.set_xlabel(r"$t$")

    title = r"$\eta=-2,\eta_{AB}=-\eta_{BA}=0.5, D_r=0.1, \rho_0=\bar{\rho}_{A,B}$"
    fig.suptitle(title)

    add_line(ax1, 0.75, 0.65, 0.99, 0.33, label="0.33", xl=0.93, yl=0.7, fontsize="large")
    add_line(ax2, 0.75, 0.65, 0.99, 0.25, label="0.25", xl=0.85, yl=0.65, fontsize="large")


    plt.show()
    plt.close()


if __name__ == "__main__":
    # f0 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
    # show_t(f0)

    # plot_rhoA_B_Sq(f0, [0, 8, 21, 33, 54, 74], vmin=None, vmax=(40, 100))

    # f0 = "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    # show_t(f0)
    # plot_rhoA_B_Sq(f0, [0, 8, 21, 33, 49], vmin=None, vmax=(80, 80))

    # f0 = "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p14_10_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    # show_t(f0)
    # plot_rhoA_B_Sq(f0, [0, 8, 21, 33, 54], vmin=None, vmax=(80, 80))

    # f0 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
    # show_t(f0)
    # plot_rhoA_B_Sq(f0, [0, 8, 21, 33, 50], vmin=None, vmax=(100, 100))


    # plot_Rt_L320_r16()

    # plot_Rt_SP_osi_band()

    # plot_Rt_SP_spiral_cluster()

    # plot_Rt_L160()

    plot_Rt_bubble_burst()
    # plot_Rt_varied_composition()

    # plot_Rt_one_5th()

    # plot_Rt_one_3rd()

