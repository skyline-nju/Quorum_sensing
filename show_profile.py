import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_one_frame
from gsd import hoomd


def get_profile(snap: hoomd.Snapshot, bins: int, dt=0.01):
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    dx = Lx / bins
    bin_area = dx * Ly
    x = np.arange(bins) * dx + dx / 2

    n_types = len(snap.particles.types)
    rho_x = np.zeros((n_types, bins))
    p_x = np.zeros_like(rho_x)
    speed = np.zeros_like(rho_x)
    m_x = np.zeros_like(rho_x)
    t = snap.configuration.step * dt

    for i in range(n_types):
        mask = snap.particles.typeid == i
        pos_x = snap.particles.position[:, 0][mask] + 0.5 * Lx
        theta = snap.particles.position[:, 2][mask]
        vel = snap.particles.charge[mask]
        ux = np.cos(theta)
        for j in range(pos_x.size):
            idx = int(pos_x[j] / dx)
            if idx == bins:
                idx = 0
            rho_x[i, idx] += 1
            p_x[i, idx] += ux[j]
            speed[i, idx] += vel[j]
            m_x[i, idx] += ux[j] * vel[j]
    p_x /= rho_x
    speed /= rho_x
    rho_x /= bin_area
    m_x /= bin_area
    return t, x, rho_x, p_x, speed, m_x


def plot_instant_profile(fname, i_frame=-1, bins=40):
    snap = read_one_frame(fname, i_frame)
    t, x, rho_x, p_x, speed, m_x = get_profile(snap, bins=bins)

    para = get_para(fname)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,
                                             1,
                                             sharex=True,
                                             constrained_layout=True,
                                             figsize=(8, 6))
    ax1.plot(x, rho_x[0], label=r"$A$", c="tab:blue")
    ax1.plot(x, rho_x[1], label=r"$B$", c="tab:red")
    ax1.axhline(para["rho0"], linestyle="dashed", c="k")
    ax1.set_xlim(0, x[-1] + x[0])
    ax1.set_ylim(0)
    ax1.legend()
    ax1.set_title(r"$t=%g$" % t)


    ax2.plot(x, speed[0], label=r"$A$", c="tab:blue")
    ax2.plot(x, speed[1], label=r"$B$", c="tab:red")
    ax2.axhline(1, linestyle="dashed", c="k")

    ax3.plot(x, p_x[0], label=r"$A$", c="tab:blue")
    ax3.plot(x, p_x[1], label=r"$B$", c="tab:red")
    ax3.axhline(0, linestyle="dashed", c="k")
    ax3.set_ylim(-0.75, 0.75)

    ax4.plot(x, m_x[0], label=r"$A$", c="tab:blue")
    ax4.plot(x, m_x[1], label=r"$B$", c="tab:red")

    ax4.set_xlabel(r"$x$", fontsize="x-large")
    ax1.set_ylabel(r"$\langle\rho\rangle_y$", fontsize="x-large")
    ax2.set_ylabel(r"$\langle v_m \rangle_y$", fontsize="x-large")
    ax3.set_ylabel(r"$\langle {\bf p}_x\rangle_y$", fontsize="x-large")
    ax4.set_ylabel(r"$\langle {\bf m}_x\rangle_y$", fontsize="x-large")

    ax2.legend()

    plt.savefig("fig/snap/%s_%g.png" % (os.path.basename(fname).rstrip(".gsd"), t))
    # plt.show()
    plt.close()

def profiles_varied_etaAB():
    prefix = "/scratch03.local/yduan/QS5/L40_5_k0.7_p40_beta"
    J_arr = np.array([0.1, 0.3, 0.5, 0.8])
    seed = 1200
    Lx = 40
    Ly = 5
    k = 0.7
    eta = -2.
    rho0 = 40
    bins = 40
    Dr = 0.1

    fig, axes = plt.subplots(3, 4, sharex=True, sharey="row",
        figsize=(9, 4.5), constrained_layout=True)
    for i, J in enumerate(J_arr):
        fin = f"{prefix}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{rho0}_{rho0}_r{rho0}_" \
            f"e{eta:.3f}_J{J:.3f}_{-J:.3f}_{seed}.gsd"
        print(fin)
        if J == 0.8:
            i_frame = -2
        else:
            i_frame = -1
        snap = read_one_frame(fin, i_frame)
        t, x, rho_x, p_x, speed, m_x = get_profile(snap, bins=bins)

        axes[0, i].plot(x, rho_x[0], label=r"$A$", c="tab:blue")
        axes[0, i].plot(x, rho_x[1], label=r"$B$", c="tab:red")
        axes[0, i].axhline(rho0, linestyle="dashed", c="k")
    
        axes[1, i].plot(x, speed[0], label=r"$A$", c="tab:blue")
        axes[1, i].plot(x, speed[1], label=r"$B$", c="tab:red")
        axes[1, i].axhline(1, linestyle="dashed", c="k")

        axes[2, i].plot(x, p_x[0], label=r"$A$", c="tab:blue")
        axes[2, i].plot(x, p_x[1], label=r"$B$", c="tab:red")
        axes[2, i].axhline(0, linestyle="dashed", c="k")
    
    for ax in axes[-1]:
        ax.set_xlabel(r"$x$", fontsize="x-large")
        ax.set_xlim(0, 40)
    for i, ax in enumerate(axes[0]):
        ax.set_title(r"$\eta_{AB}=%g$" % J_arr[i])

    ylabels = [r"$\langle\rho \rangle_y$", r"$\langle v_m\rangle_y$", r"$\langle \mathbf{p}_y/\rho\rangle_y$"]
    for i, ax in enumerate(axes[:, 0]):
        ax.set_ylabel(ylabels[i], fontsize="x-large")

            
    fout = f"fig/snap/40_5_e{eta:g}_Dr{Dr:g}.pdf"
    plt.savefig(fout)
    plt.close()


if __name__ == "__main__":
    # fname = "D:/data/QS/20_5_200_100_0.00_2.00_1.0_0.02_0.gsd"
    # fname = "/scratch03.local/yduan/QS5/L80_5_k0.7_p40_beta/L80_5_Dr0.020_k0.70_p40_40_r40_e-2.000_J0.300_-0.300_1200.gsd"
    # plot_instant_profile(fname, bins=160, i_frame=-1)

    profiles_varied_etaAB()
