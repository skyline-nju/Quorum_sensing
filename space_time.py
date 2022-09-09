import numpy as np
import matplotlib.pyplot as plt
import os


def band_osi_static(Lx):
    if Lx == 120:
        prefix = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan/QS_MPI_noself/L120_40_r10_10_10/coarse_grain_dx1"
    elif Lx == 240:
        prefix = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan/QS_MPI_noself/L240_40_r10_10_10/coarse_grain_dx1"
    Ly = 40
    Dr=0.1
    k=0.7
    pA = 10
    rho0 = bar_rho_A = bar_rho_B = 10
    eta=-2
    J=0.5
    seed=1000

    dt = 0.01
    pB_arr = np.array([28, 34, 40])

    v0 = [5, 23]
    v1 = [25, 70]
    # v0 = [None, None]
    # v1 = [None, None]

    fig, axes = plt.subplots(2, 3, figsize=(8, 6), constrained_layout=True, sharex=True, sharey=True)
    for i, pB in enumerate(pB_arr):
        fin = f"{prefix}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
            f"r{rho0:g}_{bar_rho_A:g}_{bar_rho_B:g}_" \
            f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_{seed:d}.npz"
        with np.load(fin, "r") as data:
            t, x, y = data["t"] * dt, data["x"], data["y"]
            fields = data["fields"]
            rho_A = fields[:, 0]
            rho_B = fields[:, 1]

            rho_A_x = np.mean(rho_A, axis=1)
            rho_B_x = np.mean(rho_B, axis=1)

        extent = [0, Lx, 0, t[-1]]
        im0 = axes[0, i].imshow(rho_A_x, origin="lower", extent=extent, aspect="auto", vmin=v0[0], vmax=v0[1])
        im1 = axes[1, i].imshow(rho_B_x, origin="lower", extent=extent, aspect="auto", vmin=v1[0], vmax=v1[1])
        axes[0, i].set_title(r"$\phi_B=%g$" % pB)
        axes[1, i].set_xlabel(r"$x$", fontsize="large")

    axes[0, 0].set_ylabel(r"$t$", fontsize="large")
    axes[1, 0].set_ylabel(r"$t$", fontsize="large")

    cb0 = plt.colorbar(im0, ax=axes[0, -1], extend="both")
    cb1 = plt.colorbar(im1, ax=axes[1, -1], extend="both")
    cb0.set_label(r"$\langle \rho_{A}(\mathbf{r})\rangle_y$", fontsize="large")
    cb1.set_label(r"$\langle \rho_{B}(\mathbf{r})\rangle_y$", fontsize="large")

    title = r"$L_x=%g, L_y=%g, D_r=%g, \eta=%g,\eta_{AB}=-\eta_{BA}=%g, \rho_0=\bar{\rho}_A=\bar{\rho}_B=%g, \phi_A=%g$" % (
        Lx, Ly, Dr, eta, J, rho0, pA)
    fig.suptitle(title)
    plt.show()
    plt.close()



def show_osi_snaps():
    prefix = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan/QS_MPI_noself/L120_40_r10_10_10/coarse_grain_dx1"

    Lx = 120
    Ly = 40
    Dr=0.1
    k=0.7
    pA = 10
    rho0 = bar_rho_A = bar_rho_B = 10
    eta=-2
    J=0.5
    seed=1000
    dt = 0.01
    frames = np.arange(9) * 8 + 250

    pB = 34
    fig, axes = plt.subplots(9, 2, figsize=(6, 9),
                             constrained_layout=True,
                             sharex=True, sharey=True)
    vA = [0, 35]
    vB = [0, 80]

    fin = f"{prefix}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
        f"r{rho0:g}_{bar_rho_A:g}_{bar_rho_B:g}_" \
        f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_{seed:d}.npz"
    with np.load(fin, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
        rho_A = fields[:, 0]
        rho_B = fields[:, 1]
    extent = [0, Lx, 0, Ly]
    for i, i_frame in enumerate(frames):
        im0 = axes[i, 0].imshow(rho_A[i_frame], origin="lower", extent=extent, vmin=vA[0], vmax=vA[1])
        im1 = axes[i, 1].imshow(rho_B[i_frame], origin="lower", extent=extent, vmin=vB[0], vmax=vB[1])
        axes[i, 0].set_xticks([])
        axes[i, 0].set_yticks([])
        axes[i, 1].set_xticks([])
        axes[i, 1].set_yticks([])
        axes[i, 0].set_ylabel(r"$t=%g$" % t[i_frame], rotation="horizontal", ha="right", fontsize="large")
    cb0 = plt.colorbar(im0, ax=axes[-1, 0], orientation="horizontal", extend="max")
    cb1 = plt.colorbar(im1, ax=axes[-1, 1], orientation="horizontal", extend="max")
    cb0.set_label(r"$\rho_A$", fontsize="large")
    cb1.set_label(r"$\rho_B$", fontsize="large")

    title = r"$L_x=%g, L_y=%g, \phi_A=%g, \phi_B=%g$" % (
        Lx, Ly, pA, pB)
    fig.suptitle(title, fontsize="large")
    plt.show()
    plt.close()


if __name__ == "__main__":
    band_osi_static(Lx=120)
    show_osi_snaps()