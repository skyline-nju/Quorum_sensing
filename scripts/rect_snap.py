import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from plot_snap import plot_one_panel, add_xlabel, add_ylabel


def eta_vs_J_L160_5():
    fig = plt.figure(figsize=(10, 8))
    subfigs = fig.subfigures(3, 2, wspace=0.01, hspace=0.01)
    
    eta_arr = np.array([0, -0.2, -0.4, -0.6, -1.0, -2.0])
    J_arr = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0])
    
    Lx = 160
    Ly = 5
    Dr = 0.1
    rho0 = pA = pB = 40
    k = 0.7
    seed = 1200
    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/L{Lx}_{Ly}_k{k:g}_p{rho0}_beta"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
        f"r{rho0:g}_e%.3f_J%.3f_%.3f_{seed:d}.gsd"
    fout_name = f"fig/PD/travelling_bands/L{Lx}_{Ly}_Dr{Dr:g}_r{rho0:g}_s{seed}.png"

    axes_list = []
    for i, subfig in enumerate(subfigs.flat):
        if i < eta_arr.size:
            eta = eta_arr[i]
            axes = subfig.subplots(J_arr.size, 1, sharex=True, sharey=True)
            for j, J in enumerate(J_arr):
                fin = file_pat % (eta, J, -J)
                print(fin)
                plot_one_panel(fin, ax=axes[j], ms=0.05)
                axes[j].axis("off")
            subfig.subplots_adjust(hspace=0.0001)
            # subfig.suptitle(r"$\eta=%g$" % eta, fontsize="x-large")
            axes_list.append(axes)

    rect = [0.1, -0.03, 1.005, 0.98] # left, bottom, right, top
    plt.tight_layout(h_pad=0.15, w_pad=0.05, rect=rect)
    
    for i, axes in enumerate(axes_list):
        subfig = subfigs.flat[i]
        for j, ax in enumerate(axes):
            # bbox = [xmin, ymin, xmax, ymax]
            bbox = ax.get_position().get_points().flatten()
            y = (bbox[1] + bbox[3]) * 0.5
            x = bbox[0]
            subfig.text(x, y, r"$%g$" % J_arr[j], fontsize=15, ha="right", va="center")
            if i % subfigs.shape[1] == 0 and j == J_arr.size // 2:
                subfig.text(0, y, r"$\eta_{AB}$", fontsize=15, ha="left", va="center",
                    rotation="vertical")
            if j == 0:
                x = (bbox[0] + bbox[2]) * 0.5
                y = bbox[1] + 0.02
                subfig.text(x, y, r"$\eta=%g$" % eta_arr[i], fontsize=15,
                    ha="center", va="bottom")

    plt.savefig(fout_name)
    plt.close()


def etaAB_vs_replicas(Lx=320, Ly=5, eta=-0.4, Dr=0.1, rho0=40, pA=None, pB=None, k=0.7):
    if pA is None:
        pA = rho0
    if pB is None:
        pB = rho0
    prefix = "/scratch03.local/yduan/QS5"
    if Lx == 320 and Ly == 5:
        seed_arr = np.array([1200, 1201, 1202, 1203, 1204])
        if eta == -0.4 or eta == -0.5:
            J_arr = np.array([0.5])
        elif eta == -0.6:
            J_arr = np.array([0.4, 0.5])
        folder = f"{prefix}/L{Lx}_{Ly}_k{k:g}_JAB"
        file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
            f"r{rho0:g}_e{eta:.3f}_J%.3f_%.3f_%d.gsd"
        fout_name = f"fig/PD/travelling_bands/L{Lx}_{Ly}_Dr{Dr:g}_r{rho0:g}_e{eta:.3f}.png"

    fig = plt.figure(figsize=(12, 4))
    subfigs = fig.subfigures(2, 1, wspace=0.01, hspace=0.01)

    axes_list = []
    for i, subfig in enumerate(subfigs.flat): 
        J = J_arr[i]
        axes = subfig.subplots(seed_arr.size, 1, sharex=True, sharey=True)
        for j, seed in enumerate(seed_arr):
            fin = file_pat % (J, -J, seed)
            print(fin)
            plot_one_panel(fin, ax=axes[j], ms=0.05)
            axes[j].axis("off")
        subfig.subplots_adjust(hspace=0.0001)
        # subfig.suptitle(r"$\eta=%g$" % eta, fontsize="x-large")
        axes_list.append(axes)

    rect = [0.03, -0.05, 1.01, 0.96] # left, bottom, right, top
    plt.tight_layout(h_pad=0.15, w_pad=0.05, rect=rect)

    for i, axes in enumerate(axes_list):
        subfig = subfigs.flat[i]
        for j, ax in enumerate(axes):
            # bbox = [xmin, ymin, xmax, ymax]
            bbox = ax.get_position().get_points().flatten()
            y = (bbox[1] + bbox[3]) * 0.5
            x = bbox[0]
            subfig.text(x, y, r"$%g$" % j, fontsize=15, ha="right", va="center")
            if j == seed_arr.size // 2:
                subfig.text(0, y, r"replica", fontsize=15, ha="left", va="center",
                    rotation="vertical")
            if j == 0:
                x = (bbox[0] + bbox[2]) * 0.5
                y = bbox[1] + 0.1
                subfig.text(x, y, r"$\eta_{AB}=%g$" % J_arr[i], fontsize=15,
                    ha="center", va="bottom")

    plt.savefig(fout_name)
    plt.close()


def varied_Dr(Lx=20, Ly=5, eta=0, alpha=10, rho0=80, phiA=None, phiB=None, k=0.7, seed=1001):
    if phiA is None:
        phiA = rho0
    if phiB is None:
        phiB = rho0
    
    prefix = "/home/yduan/data/QS5/last/L20_5_k0.7_alpha"
    pat = f"{prefix}/L{Lx}_{Ly}_Dr%.3f_k{k:.2f}_p{phiA:g}_{phiB:g}_r{rho0:g}_" \
            f"e{eta:.3f}_a{alpha:.3f}_{seed}.gsd"
    Dr_arr = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.2, 0.4, 0.6,
        0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0])
    
    fig, axes = plt.subplots(4, 4, figsize=(10, 4.0), constrained_layout=True)

    for i, ax in enumerate(axes.flat):
        fin = pat % (Dr_arr[i])
        print(fin)
        plot_one_panel(fin, ax=ax, ms=0.08)
        ax.axis("off")
        ax.set_title(r"$D_r=%g$" % Dr_arr[i])
    title = f"$L_x={Lx}, L_y={Ly}, \\alpha={alpha}, \\eta={eta}," \
        f"\\phi_A=\\phi_B=\\rho_0={rho0},\\kappa={k}$"
    plt.suptitle(title, fontsize="x-large")
    # plt.show()
    plt.savefig("fig/snap/L20_5_a%g_e0_varied_Dr.png" % alpha)
    plt.close()


def varied_alpha(Lx=20, Ly=5, eta=0, Dr=0.02, rho0=80, phiA=None, phiB=None, k=0.7, seed=1001):
    if phiA is None:
        phiA = rho0
    if phiB is None:
        phiB = rho0

    prefix = "/home/yduan/data/QS5/last/L20_5_k0.7_alpha"

    pat = f"{prefix}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{phiA:g}_{phiB:g}_r{rho0:g}_" \
            f"e{eta:.3f}_a%.3f_{seed}.gsd"
    Dr_arr = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    
    fig, axes = plt.subplots(4, 4, figsize=(10, 4.0), constrained_layout=True)

    for i, ax in enumerate(axes.flat):
        fin = pat % (Dr_arr[i])
        print(fin)
        plot_one_panel(fin, ax=ax, ms=0.08)
        ax.axis("off")
        ax.set_title(r"$\alpha=%g$" % Dr_arr[i])
    title = f"$L_x={Lx}, L_y={Ly}, D_r={Dr}, \\eta={eta}," \
        f"\\phi_A=\\phi_B=\\rho_0={rho0},\\kappa={k}$"
    plt.suptitle(title, fontsize="x-large")
    # plt.show()
    plt.savefig("fig/snap/L20_5_Dr%g_e0_varied_Dr.png" % Dr)
    plt.close()

     
if __name__ == "__main__":
    # eta_vs_J_L160_5()
    # etaAB_vs_replicas(eta=-0.4)

    # varied_Dr(alpha=1)
    varied_alpha()
