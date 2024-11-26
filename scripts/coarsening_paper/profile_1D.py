import numpy as np
import matplotlib.pyplot as plt


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def time_ave_profile():
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1"
    fname = f"{folder}/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_21000160.npz"

    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]

    rho_A = np.mean(fields[1:, 0], axis=0)
    rho_B = np.mean(fields[1:, 1], axis=0)
    
    plt.plot(x, rho_A)
    plt.plot(x, rho_B)
    plt.show()
    plt.close()


def show_space_time():
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1"
    fname = f"{folder}/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"

    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]

    rhoA = fields[160:, 0]
    rhoB = fields[160:, 1]

    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 6), constrained_layout=True)
    ax1.imshow(rhoA, origin="lower", aspect="auto")
    ax2.imshow(rhoB, origin="lower", aspect="auto")
    ax1.axvline(150)
    ax1.axvline(250)
    ax1.axvline(750)
    ax1.axvline(850)

    plt.show()
    plt.close()

    rhoA_left = np.mean(rhoA[:, 150:250], axis=1)
    rhoA_right = np.mean(rhoA[:, 750:850], axis=1)

    plt.plot(t[160:], rhoA_left)
    plt.plot(t[160:], rhoA_right)
    plt.show()
    plt.close()
 

if __name__ == "__main__":
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx4"
    fname = f"{folder}/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_21000220.npz"

    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]
        fields_name = data["fields_name"]

    rhoA = np.mean(fields[:, 0], axis=0)
    rhoB = np.mean(fields[:, 1], axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 4), constrained_layout=True)
    ax1.plot(x, rhoA)
    ax2.plot(x, rhoB)
    ax1.set_ylabel(r"$\rho_A$", fontsize="x-large")
    ax2.set_ylabel(r"$\rho_B$", fontsize="x-large")
    ax2.set_xlabel(r"$x$", fontsize="x-large")
    ax1.set_xlim(0, 1280)
    ax1.set_ylim(0, 4.5)
    ax2.set_ylim(0, 4.5)

    plt.show()
    plt.close()

