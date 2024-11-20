import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import sys
sys.path.append("../../")

from structure_factor import get_h, get_Sq_t_rho, get_q1, get_npz_file_name
from add_line import add_line


def plot_Rt(ax):
    para = {}
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["Dt"] = 0
    para["eta"] = -2
    para["eta_AB"] = 0.5

    para["h"] = 1 / 3.

    para["Lx"] = 1280
    para["Ly"] = para["Lx"]
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["pA"] = 5.5
    para["pB"] = 3.85

    folder = "/home/yduan/code/Quorum_sensing"
    f_Sqt = get_npz_file_name(para)
    f_Sqt = f"{folder}/{f_Sqt}"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])

    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
    mask = t < 3e6
    ax.plot(t[mask], 2 * np.pi/q1_A[mask], "o", fillstyle="none", c="tab:pink", label=r"$(0.55, 0.385)$")

    para["pA"] = 7.5
    para["pB"] = 10
    f_Sqt = get_npz_file_name(para)
    f_Sqt = f"{folder}/{f_Sqt}"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
    mask = t<1.2e6
    ax.plot(t[mask], 2 * np.pi/q1_A[mask], "o", fillstyle="none", c="tab:cyan", label=r"$(0.76, 1.06)$")


    para["pA"] = 10
    para["pB"] = 17.4

    f_Sqt = get_npz_file_name(para)
    f_Sqt = f"{folder}/{f_Sqt}"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
    ax.plot(t, 2 * np.pi/q1_B, "o", fillstyle="none", c="tab:brown", label=r"$(1, 1.74)$")
    for i in range(t.size):
        print(i, t[i])


    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.65, 0.7, 0.95, 0.33, label="0.33", xl=0.75, yl=0.75, fontsize="large")
    # add_line(ax, 0.75, 0.6, 0.98, 0.2, label="0.2", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax, 0.1, 0.3, 0.95, 0.13, label="0.13", xl=0.8, yl=0.5, fontsize="large")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$\zeta$", fontsize="x-large")
    ax.legend(title=r"$(\bar{\rho}_A/\rho_0, \bar{\rho}_B/\rho_0)$", fontsize="large", title_fontsize="x-large")


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    plot_Rt(axes[0])


    snaps = ["fig/snaps/CCB_G_2.png", "fig/snaps/CCB_G_3.png"]

    for i, ax in enumerate(axes[1:]):
        im = mpimg.imread(snaps[i])
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
    # plt.show()
    plt.savefig("fig/Fig2.pdf")
    plt.close()
