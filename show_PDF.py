import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_one_frame
from gsd import hoomd


def get_PDF(snap: hoomd.Snapshot, dx: float):
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    bin_area = dx ** 2
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)

    par_num = np.zeros((2, nrows, ncols), int)
    mx = np.zeros((2, nrows, ncols))
    my = np.zeros((2, nrows, ncols))
    v = np.zeros((2, nrows, ncols))

    x = snap.particles.position[:,0] + Lx / 2
    y = snap.particles.position[:,1] + Ly / 2

    for i in range(snap.particles.N):
        col = int(x[i]/dx)
        row = int(y[i]/dx)
        if col < 0:
            col += ncols
        elif col >= ncols:
            col -= ncols
        if row < 0:
            row += nrows
        elif row >= nrows:
            row -= nrows

        id_type = snap.particles.typeid[i]
        par_num[id_type, row, col] += 1
        theta = snap.particles.position[i, 2]
        mx[id_type, row, col] += np.cos(theta)
        my[id_type, row, col] += np.sin(theta)
        v[id_type, row, col] += snap.particles.charge[i]

    mask = par_num > 0
    v[mask] /= par_num[mask]

    rho = par_num / bin_area
    m = np.sqrt(mx ** 2 + my ** 2) / bin_area
    return rho, m, v


def tmp():
    data = np.load("data/PDF/DP3.npz")
    rho_PDF = data["rho_PDF"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    hist, bin_edges = np.histogram(rho_PDF[0], bins=20, density=True)
    rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    ax1.plot(rho, hist)
    hist, bin_edges = np.histogram(rho_PDF[1], bins=20, density=True)
    rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    ax1.plot(rho, hist)

    H, xedges, yedges = np.histogram2d(rho_PDF[1].flatten(), rho_PDF[0].flatten(), bins=20, density=True)
    im2 = ax2.imshow(H, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plt.colorbar(im2, ax=ax2)
    plt.show()
    plt.close()


def plot_profile(i_frame, npzdata=None, fin=None):
    if npzdata is None:
        npzdata = np.load(fin, "r")
    t, x, y = npzdata["t"], npzdata["x"], npzdata["y"]
    fields = npzdata["fields"]

    fig, axes = plt.subplots(4, 1, figsize=(8, 6), constrained_layout=True, sharex=True)
    ylabels = [r"$\langle \rho\rangle_y$",
        r"$\langle \mathbf{m}_x \rangle_y$",
        r"$\langle \mathbf{m}_y\rangle_y$",
        r"$\langle v\rangle_y$"]
    for (i, ax) in enumerate(axes):
        fx_A = np.mean(fields[i_frame][i * 2], axis=0)
        fx_B = np.mean(fields[i_frame][i * 2 + 1], axis=0)
        ax.plot(x, fx_A, label="A")
        ax.plot(x, fx_B, label="B")
        ax.set_ylabel(ylabels[i])
        ax.legend()
    
    ax.set_xlabel(r"$x$")
    plt.show()
    plt.close()


def plot_density_fields(i_frame, npzdata=None, fin=None):
    if npzdata is None:
        npzdata = np.load(fin, "r")
    t, x, y = npzdata["t"], npzdata["x"], npzdata["y"]
    fields = npzdata["fields"]
    Lx = x[0] + x[-1]
    Ly = y[0] + y[-1]
    extent = [0, Lx, 0, Ly]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True, sharex=True)
    im1 = ax1.imshow(fields[i_frame][0], origin="lower", extent=extent)
    im2 = ax2.imshow(fields[i_frame][1], origin="lower", extent=extent)
    ax2.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax2.set_ylabel(r"$y$")
    cb1 = plt.colorbar(im1, ax=ax1)
    cb2 = plt.colorbar(im2, ax=ax2)
    cb1.set_label(r"$\rho_A$", fontsize="x-large")
    cb2.set_label(r"$\rho_B$", fontsize="x-large")

    plt.show()
    plt.close()


def plot_density_PDF(i_frame, npzdata=None, fin=None, ncols=2):
    if npzdata is None:
        npzdata = np.load(fin, "r")
    t, x, y = npzdata["t"], npzdata["x"], npzdata["y"]
    fields = npzdata["fields"]

    if ncols == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
        ax1.plot(fields[i_frame][0], fields[i_frame][1], ".", c="tab:blue", ms=0.2)
        # hist_range = [[0, 150], [0, 120]]
        hist_range = None
        H, xedges, yedges = np.histogram2d(
            fields[i_frame][1].flatten(), fields[i_frame][0].flatten(), bins=40, range=hist_range, density=True)
        print(H.shape)
        extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]]
        im2 = ax2.imshow(H, origin='lower', extent=extent, vmax=0.002)
        plt.colorbar(im2, ax=ax2)
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
        hist, bin_edges = np.histogram(fields[i_frame][0], bins=20, density=True)
        rho = (bin_edges[1:] + bin_edges[:-1]) / 2
        ax1.plot(rho, hist)
        hist, bin_edges = np.histogram(fields[i_frame][1], bins=20, density=True)
        rho = (bin_edges[1:] + bin_edges[:-1]) / 2
        ax1.plot(rho, hist)

        ax2.plot(fields[i_frame][0], fields[i_frame][1], ".", c="tab:blue", ms=0.2)

        # hist_range = [[0, 150], [0, 120]]
        hist_range = None
        H, xedges, yedges = np.histogram2d(
            fields[i_frame][1].flatten(), fields[i_frame][0].flatten(), bins=40, range=hist_range, density=True)
        print(H.shape)
        extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]]
        im3 = ax3.imshow(H, origin='lower', extent=extent, vmax=0.002)
        plt.colorbar(im3, ax=ax3)

    plt.show()
    plt.close()


if __name__ == "__main__":
    from datetime import datetime
    prefix = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan/QS_MPI/L160_e-2.0_J0.5/coarse_grain_dx4"
    # prefix = "/home/yduan/data/tmp"
    basename = "L160_160_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_1000.npz"
    fin = f"{prefix}/{basename}"

    data = np.load(fin, "r")
    t, x, y = data["t"], data["x"], data["y"]
    fields = data["fields"]

    i_frame = -2

    plot_density_fields(i_frame, npzdata=data)
    plot_profile(i_frame, npzdata=data)

    plot_density_PDF(i_frame, npzdata=data)

    

    