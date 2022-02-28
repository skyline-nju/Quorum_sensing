import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_nframes, read_frames, read_one_frame
from gsd import hoomd


def get_rho_x(snap: hoomd.Snapshot, bins: int):
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    bin_area = Lx * Ly / bins

    n_types = len(snap.particles.types)
    rho_x = np.zeros((n_types, bins))
    for i in range(n_types):
        mask = snap.particles.typeid == i
        pos_x = snap.particles.position[:, 0][mask] + 0.5 * Lx
        hist, bin_edges = np.histogram(pos_x,
                                       bins=bins,
                                       range=(0, Lx),
                                       density=False)
        rho_x[i] = hist / bin_area
        x = (bin_edges[:-1] + bin_edges[1:]) * 0.5
    return x, rho_x


def get_rhox_px(snap: hoomd.Snapshot, bins: int):
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    dx = Lx / bins
    bin_area = dx * Ly
    x = np.arange(bins) * dx + dx / 2

    n_types = len(snap.particles.types)
    rho_x = np.zeros((n_types, bins))
    p_x = np.zeros_like(rho_x)

    for i in range(n_types):
        mask = snap.particles.typeid == i
        pos_x = snap.particles.position[:, 0][mask] + 0.5 * Lx
        theta = snap.particles.position[:, 2][mask]
        ux = np.cos(theta)
        for j in range(pos_x.size):
            idx = int(pos_x[j] / dx)
            if idx == bins:
                idx = 0
            rho_x[i, idx] += 1
            p_x[i, idx] += ux[j]
    p_x /= rho_x
    rho_x /= bin_area
    return x, rho_x, p_x


def get_rho_x_t(fname: str, bins: int):
    n_frames = get_nframes(fname)
    frames = read_frames(fname)
    bins = 40
    rho_x = np.zeros((n_frames, 2, bins))
    for i, snap in enumerate(frames):
        x, rho_x[i] = get_rho_x(snap, bins)
    fout = os.path.basename(fname).replace(".gsd", ".npz")
    np.savez_compressed(f"data/{fout}", x=x, rho_x=rho_x)


def plot_instant_profile(fname, rho0, i_frame=-1, bins=40):
    snap = read_one_frame(fname, i_frame)
    x, rho_x, p_x = get_rhox_px(snap, bins=bins)

    fig, (ax1, ax2) = plt.subplots(2,
                                   1,
                                   sharex=True,
                                   constrained_layout=True,
                                   figsize=(8, 5))
    ax1.plot(x, rho_x[0], label=r"$A$", c="tab:blue")
    ax1.plot(x, rho_x[1], label=r"$B$", c="tab:red")
    ax1.axhline(rho0, linestyle="dashed", c="k")
    ax1.set_xlim(0, x[-1] + x[0])
    ax1.set_ylim(0)
    ax1.legend()
    ax1.set_title(r"$t=%g$" % (i_frame * 10))
    # plt.savefig("data/%s_%g.png" % (basename.rstrip(".gsd"), i_frame))

    ax2.plot(x, p_x[0], label=r"$A$", c="tab:blue")
    ax2.plot(x, p_x[1], label=r"$B$", c="tab:red")
    ax2.axhline(0, linestyle="dashed", c="k")

    ax2.set_ylim(-0.75, 0.75)

    ax2.set_xlabel(r"$x$", fontsize="x-large")
    ax2.set_ylabel(r"$\langle {\bf p}_x\rangle_y$", fontsize="x-large")
    ax1.set_ylabel(r"$\langle\rho\rangle_y$", fontsize="x-large")
    ax2.legend()

    plt.show()
    plt.close()


def plot_space_time(bins=40):
    # folder = '/scratch03.local/yduan/QS2_varied_Dr/L20_5_Dr0.1_h0.01_e-4'
    # basename = "20_5_200_100_-4.00_3.60_1.0_0.1_0.001.gsd"
    # folder = "/scratch03.local/yduan/QS2_varied_Dr/four_bands"
    # basename = "40_5_200_100_0.00_2.00_1.0_0.02_0.gsd"

    # folder = '/scratch03.local/yduan/QS2_varied_Dr/four_bands'
    # basename = "20_5_200_100_0.00_2.00_1.0_0.02_0.gsd"

    folder = "/scratch03.local/yduan/QS2_varied_Dr/restart_h0.01"
    basename = "L40_5_pA100_pB80_r100_e0.000_a0.600_Dr0.02_Dt0.gsd"
    fname = f"{folder}/{basename}"

    # get_rho_x_t(fname, bins)

    fnpz = "data/" + os.path.basename(fname).replace(".gsd", ".npz")
    with np.load(fnpz) as data:
        # x = data["x"]
        rho_x = data["rho_x"]
        print(rho_x.shape)

    plt.imshow(rho_x[:, 1, :], origin="lower", aspect="auto", vmax=200)
    plt.show()
    plt.close()


if __name__ == "__main__":
    # fname = "D:/data/QS/20_5_200_100_0.00_2.00_1.0_0.02_0.gsd"
    fname = "D:/data/QS/L20_5_Dr1_Dt0_k0.70_pA80_pB80_r80_e-2.000_J0.200_-0.200.gsd"
    plot_instant_profile(fname, rho0=100, bins=40, i_frame=400)
