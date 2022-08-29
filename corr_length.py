import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from add_line import add_line


def get_Lx(fname):
    basename = os.path.basename(fname)
    Lx = float(basename.lstrip("L").split("_")[0])
    return Lx


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
            axes[2, i].set_xlabel(r"$q$")
            if i == len(i_frame) - 1:
                cb1 = plt.colorbar(im1, ax=axes[0], extend=get_cb_extend(vmin[0], vmax[0]))
                cb2 = plt.colorbar(im2, ax=axes[1], extend=get_cb_extend(vmin[1], vmax[1]))
                # cb1.set_label(r"$\rho_A$", fontsize="large")
                # cb2.set_label(r"$\rho_B$", fontsize="large")
        axes[0, 0].set_ylabel(r"$\rho_A(\mathbf{r}, t)$")
        axes[1, 0].set_ylabel(r"$\rho_B(\mathbf{r}, t)$")
        axes[2, 0].set_ylabel(r"$\tilde{S}(q, t)$")
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


def plot_Rt_L640():
    fins = ["/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p12_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p14_10_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p0_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
        ]

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    labels = [[r"$R_A:\phi_A=20,\phi_B=20$", r"$R_B:\phi_A=20,\phi_B=20$"],
              [r"$R_A:\phi_A=12,\phi_B=20$", r"$R_B:\phi_A=12,\phi_B=20$"],
              [r"$R_A:\phi_A=14,\phi_B=10$", r"$R_B:\phi_A=14,\phi_B=10$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)
        if i < len(fins)-1:
            line, = ax.plot(t, Rt_A, "o", fillstyle="none", ms=5, label=labels[i][0])
            ax.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=5, label=labels[i][1])
            # ax.plot(t, Rt_tot, ">", c=line.get_c(), fillstyle="none")
        else:
            label=r"$R_B:\phi_A=0,\phi_B=20$"
            ax.plot(t, Rt_B, "s", c="tab:grey", fillstyle="none", ms=5, label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$R$", fontsize="large")
    ax.legend(fontsize="large")

    title = r"$L=640, \eta=-2, \eta_{\rm AB}=-\eta_{\rm BA}=0.5, \rho_0=\bar{\rho}_A=\bar{\rho}_B=20$"
    ax.set_title(title, fontsize="x-large")

    add_line(ax, 0.3, 0.25, 0.95, 1/3, label="0.33", xl=0.5, yl=0.5, fontsize="large")
    add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_L320_pA5():
    fins = ["/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p5_34_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p5_35_r20_20_20_e-2.000_J2.000_-2.000_1000.npz",
            # "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            # "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p20_20_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
        ]

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    labels = [[r"$R_A:\phi_A=5,\phi_B=33$", r"$R_B:\phi_A=5,\phi_B=33$"],
              [r"$R_A:\phi_A=5,\phi_B=34$", r"$R_B:\phi_A=5,\phi_B=34$"],
              [r"$R_A:\phi_A=5,\phi_B=35$", r"$R_B:\phi_A=5,\phi_B=35$"]]
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

    add_line(ax, 0.2, 0.0, 0.6, 0.2, label="0.2", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_L320_r10():
    ms=4
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5), constrained_layout=True, sharey=True, sharex=True)
    prefix = "/home/yduan/data/tmp"
    fins = [f"{prefix}/L320_320_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{prefix}/L320_320_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{prefix}/L320_320_Dr0.100_k0.70_p10_40_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
        ]
    labels = [[r"$R_A:\phi_B=28$", r"$R_B:\phi_B=28$"],
              [r"$R_A:\phi_B=34$", r"$R_B:\phi_B=34$"],
              [r"$R_A:\phi_B=40$", r"$R_B:\phi_B=40$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)

        line, = ax1.plot(t, Rt_A, "o", fillstyle="none", ms=ms, label=labels[i][0])
        ax1.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=ms, label=labels[i][1])
        # ax.plot(t, Rt_tot, ">", c=line.get_c(), fillstyle="none")

    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$R$", fontsize="large")
    ax1.legend(fontsize="large")

    fins = [f"{prefix}/L160_160_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{prefix}/L320_320_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=160$", r"$R_B:L=160$"],
              [r"$R_A:L=320$", r"$R_B:L=320$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)
        line, = ax2.plot(t, Rt_A, "o", fillstyle="none", ms=ms, label=labels[i][0])
        ax2.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=ms, label=labels[i][1])
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax2.legend(fontsize="large")


    fins = [f"{prefix}/L160_160_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz",
            f"{prefix}/L320_320_Dr0.100_k0.70_p10_34_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"]
    labels = [[r"$R_A:L=160$", r"$R_B:L=160$"],
              [r"$R_A:L=320$", r"$R_B:L=320$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)
        line, = ax3.plot(t, Rt_A, "o", fillstyle="none", ms=ms, label=labels[i][0])
        ax3.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=ms, label=labels[i][1])
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel(r"$t$", fontsize="large")
    ax3.legend(fontsize="large")

    ax1.set_title(r"(a) $L=320$", fontsize="x-large")
    ax2.set_title(r"(b) $\phi_B=28$", fontsize="x-large")
    ax3.set_title(r"(c) $\phi_B=34$", fontsize="x-large")

    add_line(ax1, 0.3, 0.3, 0.95, 1/3, label="0.33", xl=0.5, yl=0.6, fontsize="large")

    # add_line(ax, 0.3, 0.25, 0.95, 0.3, label="0.33", xl=0.5, yl=0.5, fontsize="large")

    # add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    # add_line(ax, 0.68, 0.35, 0.8, 0.6, label="0.6", xl=0.8, yl=0.4, c="tab:green", fontsize="large")

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


def plot_Rt_cluster_spiral():
    fins = ["/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p30_10_r40_e-2.000_J2.000_-2.000_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p30_20_r40_e-2.000_J2.000_-2.000_1000.npz",
            "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p15_10_r20_e-2.000_J2.000_-2.000_1000.npz"
        ]

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    labels = [[r"$R_A:\phi_A=30,\phi_B=10$", r"$R_B:\phi_A=30,\phi_B=10$"],
              [r"$R_A:\phi_A=30,\phi_B=20$", r"$R_B:\phi_A=30,\phi_B=20$"],
              [r"$R_A:\phi_A=15,\phi_B=10$", r"$R_B:\phi_A=15,\phi_B=10$"]]
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B, Rt_tot = get_R_t(fin)

        line, = ax.plot(t, Rt_A, "o", fillstyle="none", ms=5, label=labels[i][0])
        ax.plot(t, Rt_B, "s", c=line.get_c(), fillstyle="none", ms=5, label=labels[i][1])
        # ax.plot(t, Rt_tot, ">", c=line.get_c(), fillstyle="none")

    # ax.set_xscale("log")
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


if __name__ == "__main__":
    f0 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p0_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    f1 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    f2 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.000_-0.000_1000.npz"

    f3 = "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"

    # print(os.path.split(f3))
    # show_t(f3)

    # plot_rho_A_B(f0, [0, 9, 18, 36], vmin=0, vmax=160)

    # plot_rhoA_B_Sq(f0, [0, 9, 18, 36], vmin=0, vmax=160)


    # plot_Rt_L640()

    # plot_Rt_L320_pA5()

    plot_Rt_L320_r16()

    # plot_Rt_L320_r10()

    # plot_Rt_cluster_spiral()

    # plot_Rt_L160()
