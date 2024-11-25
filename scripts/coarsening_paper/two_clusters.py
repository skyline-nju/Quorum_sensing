import numpy as np
import matplotlib.pyplot as plt


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def two_clusters(seed=2026, axes=None):
    if axes is None:
        fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True)
        ax1, ax2 = axes
        flag_show = True
    else:
        ax1, ax2 = axes
        flag_show = False

    dx = 4
    Lx = 1280
    Ly = 320
    folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx:d}_{Ly:d}_r10_10_10_h0.333/coarse_grain_dx{dx:d}"

    fname = f"{folder}/L1280_320_Dr0.100_k0.70_p5_2.41_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"

    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]

    rhoA = np.mean(fields[1:, 0], axis=0)
    rhoB = np.mean(fields[1:, 1], axis=0)

    rho_v_A = np.mean(fields[1:, 0] * fields[1:, 6], axis=0)
    rho_v_B = np.mean(fields[1:, 1] * fields[1:, 7], axis=0)
    # rho_v_A = np.mean(fields[1:2, 1] * fields[1:2, 7], axis=0)


    vx_A = np.mean(fields[1:, 8], axis=0)
    vx_B = np.mean(fields[1:, 9], axis=0)


    # ax1.plot(x, np.mean(rhoA, axis=0))
    # ax1.plot(x, np.mean(rhoB, axis=0))
    ax1.plot(x, np.mean(rho_v_A, axis=0))
    ax1.plot(x, np.mean(rho_v_B, axis=0))
    ax2.plot(x, np.mean(vx_A, axis=0))
    ax2.plot(x, np.mean(vx_B, axis=0))

    if flag_show:
        plt.show()
        plt.close()


def two_droplets():
    dx = 4
    Lx = 960
    Ly = 480
    seed = 2085
    folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx:d}_{Ly:d}_r10_10_10_h0.333/coarse_grain_dx{dx:d}"

    fname = f"{folder}/L960_480_Dr0.100_k0.70_p9.6_16.15_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"

    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]

    rhoA = np.mean(fields[1:, 0], axis=0)
    rhoB = np.mean(fields[1:, 1], axis=0)

    # rho_v_A = np.mean(fields[1:, 0] * fields[1:, 6], axis=0)
    rho_v_A = np.mean(fields[1:, 0] * fields[1:, 6], axis=0)
    rho_v_B = np.mean(fields[1:, 1] * fields[1:, 7], axis=0)


    vx_A = np.mean(fields[1:, 8], axis=0)
    vx_B = np.mean(fields[1:, 9], axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)

    # ax1.plot(x, np.mean(rhoA, axis=0))
    # ax1.plot(x, np.mean(rhoB, axis=0))
    ax1.plot(x, np.mean(rho_v_A, axis=0))
    ax1.plot(x, np.mean(rho_v_B, axis=0))
    ax2.plot(x, np.mean(vx_A, axis=0))
    ax2.plot(x, np.mean(vx_B, axis=0))
    plt.show()
    plt.close()


if __name__ == "__main__":
    # two_droplets()
    fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True)
    two_clusters(seed=2030, axes=axes)
    # two_clusters(seed=2116, axes=axes)

    plt.show()
    plt.close()

    # folder = f"{root_tahmineh}/QS_MPI_noself/MC_L640_r10_10_10_h0.333/coarse_grain_dx4"

    # fname = f"{folder}/L640_640_Dr0.100_k0.70_p9.6_16.15_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    # with np.load(fname, "r") as data:
    #     t = data["t"]
    #     x = data["x"]
    #     fields = data["fields"]
    #     fields_name = data["fields_name"]

    # rhoA = np.mean(fields[40, 0], axis=0)

    # plt.plot(x, rhoA)
    # plt.show()
    # plt.close()