import os
from posixpath import basename
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_nframes, read_frames, read_one_frame
from gsd import hoomd


def get_rho_x(snap: hoomd.Snapshot, bins: int) -> tuple[np.ndarray, np.ndarray]:
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    bin_area = Lx * Ly / bins

    n_types = len(snap.particles.types)
    rho_x = np.zeros((n_types, bins))
    for i in range(n_types):
        mask = snap.particles.typeid == i
        pos_x = snap.particles.position[:, 0][mask] + 0.5 * Lx
        hist, bin_edges = np.histogram(
            pos_x, bins=bins, range=(0, Lx), density=False)
        rho_x[i] = hist / bin_area
        x = (bin_edges[:-1] + bin_edges[1:]) * 0.5
    return x, rho_x


def get_rho_x_t(fname: str, bins:int):
    n_frames = get_nframes(fname)
    frames = read_frames(fname)
    bins = 40
    rho_x = np.zeros((n_frames, 2, bins))
    for i, snap in enumerate(frames):
        x, rho_x[i] = get_rho_x(snap, bins)
    fout = os.path.basename(fname).replace(".gsd", ".npz")
    np.savez_compressed(f"data/{fout}", x=x, rho_x=rho_x)


def plot_instant_profile(i_frame=-1, bins=40):
    # folder = '/scratch03.local/yduan/QS2_varied_Dr/L20_5_pA100_r100_Dr0.02_Dt0.001'
    # basename = "pB120_e0.000_a3.000.gsd"
    # folder = '/scratch03.local/yduan/QS2_varied_Dr/four_bands'
    # basename = "20_5_200_100_0.00_2.00_1.0_0.02_0.gsd"

    folder = "/scratch03.local/yduan/QS2_varied_Dr/restart_h0.01"
    basename = "L40_5_pA100_pB40_r100_e0.000_a0.800_Dr0.02_Dt0.gsd"
    fname = f"{folder}/{basename}"
    snap = read_one_frame(fname, i_frame)
    x, rho_x = get_rho_x(snap, bins=bins)
    plt.plot(x, rho_x[0], label=r"$\rho_A$")
    plt.plot(x, rho_x[1], label=r"$\rho_B$")
    plt.axhline(100, linestyle="dashed", c="tab:red")
    plt.xlim(0)
    plt.ylim(0)
    plt.legend()
    plt.title(r"$t=%g$" % (i_frame))
    plt.savefig("data/%s_%g.png" % (basename.rstrip(".gsd"), i_frame))
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
        x = data["x"]
        rho_x = data["rho_x"]
        print(rho_x.shape)

    plt.imshow(rho_x[:,1,:], origin="lower", aspect="auto", vmax=200)
    plt.show()
    plt.close()


if __name__ == "__main__":
    # plot_space_time(80)
    
    for i in range(0, 1900, 100):
        plot_instant_profile(i, 80)