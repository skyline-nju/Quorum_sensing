import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from add_line import add_line


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
            t, R_t_A, R_t_B = cache_data["t"], cache_data["R_t_A"], cache_data["R_t_B"]
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
        for i, t_i in enumerate(t):
            rhoA_q = np.fft.fft2(fields[i][0]-np.mean(fields[i][0]), norm="ortho")
            rhoB_q = np.fft.fft2(fields[i][1]-np.mean(fields[i][1]), norm="ortho")

            rhoA_q = np.fft.fftshift(rhoA_q)
            rhoB_q = np.fft.fftshift(rhoB_q)

            S_A = np.abs(rhoA_q) ** 2
            S_B = np.abs(rhoB_q) ** 2
            
            S_A_q = get_S_tilde(S_A, q_radial, q_module)
            S_B_q = get_S_tilde(S_B, q_radial, q_module)

            R_t_A[i] = get_corr_length(q_radial, S_A_q)
            R_t_B[i] = get_corr_length(q_radial, S_B_q)
        
        if not os.path.exists(f"{folder}/cache/Rt"):
            os.makedirs(f"{folder}/cache/Rt")
        np.savez_compressed(fcache, t=t, R_t_A=R_t_A, R_t_B=R_t_B)
        print("updating cache:", fcache)
    return t, R_t_A, R_t_B


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


    t, Rt_A, Rt_B = get_R_t(f1)
    t0, Rt_A0, Rt_B0 = get_R_t(f0)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(t, Rt_A, "o", label=r"$R_A: \phi_A=40,\phi_B=40$")
    ax.plot(t, Rt_B, "o", label=r"$R_B: \phi_A=40,\phi_B=40$")
    ax.plot(t0, Rt_B0, "o", fillstyle="none", c="tab:grey", label=r"$R_B: \phi_A=0,\phi_B=40$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$R$", fontsize="large")
    ax.legend()

    add_line(ax, 0.2, 0.2, 0.9, 1/3, label="0.33", fontsize="large")
    add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    add_line(ax, 0.65, 0.25, 0.8, 0.6, label="0.6", xl=0.7, yl=0.3, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def plot_Rt_L640():
    fins = ["/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p12_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz",
            "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p0_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"
        ]

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    # label_pat_A = r"$$"
    for i, fin in enumerate(fins):
        t, Rt_A, Rt_B = get_R_t(fin)
        if i < len(fins)-1:
            line, = ax.plot(t, Rt_A, "o")
            ax.plot(t, Rt_B, "s", c=line.get_c())
        else:
            ax.plot(t, Rt_B, "s", c="tab:grey", fillstyle="none")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$R$", fontsize="large")
    # ax.legend()

    add_line(ax, 0.2, 0.2, 0.9, 1/3, label="0.33", fontsize="large")
    add_line(ax, 0.2, 0.0, 0.6, 0.18, label="0.18", xl=0.35, yl=0.01, c="tab:pink", fontsize="large")
    add_line(ax, 0.65, 0.25, 0.85, 0.5, label="0.5", xl=0.7, yl=0.3, c="tab:green", fontsize="large")

    plt.show()
    plt.close()


def get_Lx(fname):
    basename = os.path.basename(fname)
    Lx = float(basename.lstrip("L").split("_")[0])
    return Lx

if __name__ == "__main__":
    f0 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p0_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    f1 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    f2 = "/home/yduan/data/tmp/L320_320_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.000_-0.000_1000.npz"

    f3 = "/home/yduan/data/tmp/L640_640_Dr0.100_k0.70_p20_20_r20_20_20_e-2.000_J0.500_-0.500_1000.npz"

    # print(os.path.split(f3))
    # show_t(f3)

    # plot_rho_A_B(f0, [0, 9, 18, 36], vmin=0, vmax=160)

    # plot_rhoA_B_Sq(f0, [0, 9, 18, 36], vmin=0, vmax=160)


    plot_Rt_L320()