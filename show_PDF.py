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

    rho = np.zeros((2, nrows, ncols))
    x = snap.particles.position[:,0] + Lx / 2
    y = snap.particles.position[:,1] + Ly / 2

    # print(y.min(), y.max())

    for i in range(snap.particles.N):
        col = int(x[i]/dx)
        row = int(y[i]/dx)
        id_type = snap.particles.typeid[i]
        rho[id_type, row, col] += 1
    
    rho /= (dx * dx)
    return rho


if __name__ == "__main__":
    # fin = "/scratch03.local/yduan/QS_MPI/L160_e-2.0_J0.5/L160_160_Dr0.100_k0.70_p40_40_r40_e-2.000_J0.500_-0.500_1000.gsd"
    # snap = read_one_frame(fin, -1)
    # rho_PDF = get_PDF(snap, dx=4)

    # np.savez_compressed("data/PDF/DP3.npz", rho_PDF=rho_PDF)


    data = np.load("data/PDF/DP3.npz")
    rho_PDF = data["rho_PDF"]

    hist, bin_edges = np.histogram(rho_PDF[0], bins=20)
    plt.plot(hist)
    hist, bin_edges = np.histogram(rho_PDF[1], bins=20)
    plt.plot(hist)
    plt.show()
    plt.close()

    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6), sharex=True, constrained_layout=True)

    # rho_PDF = np.load("data/PDF/DP2.npz")["rho_PDF"]
    # hist, bin_edges = np.histogram(rho_PDF[0], bins=20, density=True)
    # rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    # ax1.plot(rho, hist, label=r"$\rho_A$")
    # hist, bin_edges = np.histogram(rho_PDF[1], bins=20, density=True)
    # rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    # ax1.plot(rho, hist, label=r"$\rho_B$")
    # # ax1.set_title("(a) DP II")
    # ax1.set_ylabel("PDF")
    # ax1.axvline(40, linestyle="dashed", c="k", label=r"$\rho_0=40$")
    # ax2.axvline(40, linestyle="dashed", c="k")
    # rho_PDF = np.load("data/PDF/DP1.npz")["rho_PDF"]
    # hist, bin_edges = np.histogram(rho_PDF[0], bins=20, density=True)
    # rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    # ax2.plot(rho, hist, label=r"$\rho_A$")
    # hist, bin_edges = np.histogram(rho_PDF[1], bins=20, density=True)
    # rho = (bin_edges[1:] + bin_edges[:-1]) / 2
    # ax2.plot(rho, hist, label=r"$\rho_B$")
    # # ax2.set_title("(b) DP I")
    # ax2.set_ylabel("PDF")
    # ax2.set_xlabel(r"$\rho$")

    # ax1.legend(title="DP II")
    # ax2.legend(title="DP I")


    # plt.show()
    # plt.close()


    