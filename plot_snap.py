import os
import sys
from read_gsd import read_one_frame
import matplotlib.pyplot as plt
import numpy as np


def plot_one_panel(fname, i_frame=-1, ax=None, ms=0.05):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    
    snap = read_one_frame(fname, i_frame)
    Lx = snap.configuration.box[0]
    Ly = snap.configuration.box[1]
    x = snap.particles.position[:, 0] + Lx / 2
    y = snap.particles.position[:, 1] + Ly / 2
    type_id = snap.particles.typeid

    mask0 = type_id == 0
    mask1 = type_id == 1
    ax.plot(x[mask0], y[mask0], ".", ms=ms, c="b", alpha=1)
    ax.plot(x[mask1], y[mask1], ".", ms=ms, c="r", alpha=0.7)

    ax.axis("scaled")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)

    if flag_show:
        plt.show()
        plt.close()


def add_xlabel(fig, axes, loc, para_name, para_arr, fontsize):
    if loc == "bottom":
        j = axes.shape[0] - 1
        va = "top"
    elif loc == "top":
        j = 0
        va = "bottom"
    else:
        print("Error, loc must be 'top' or 'bottom'")
        sys.exit(1)

    for i, para in enumerate(para_arr):
        if i == 0:
            if para_name == "replica":
                text = r"replica $%d$" % (i)
            else:
                text = r"$%s=%g$" % (para_name, para)
        else:
            if para_name == "replica":
                text = r"%d" % i
            else:
                text = r"$%g$" % para
        # bbox = [xmin, ymin, xmax, ymax]
        bbox = axes[j, i].get_position().get_points().flatten()
        x = (bbox[0] + bbox[2]) * 0.5
        if loc == "top":
            y = bbox[3]
        else:
            y = bbox[1]
        fig.text(x, y, text, fontsize=fontsize, ha="center", va=va)


def add_ylabel(fig, axes, loc, para_name, para_arr, fontsize, vertical=True, reverse=False):
    if loc == "left":
        col = 0
        ha = "right"
    elif loc == "right":
        col = axes.shape[1] - 1
        ha = "left"
    else:
        print("Error, loc must be 'left' or 'right'")
        sys.exit(1)
    if vertical:
        rot = "vertical"
    else:
        rot = "horizontal"

    for j, para in enumerate(para_arr):
        if j == 0:
            text = r"$%s=%g$" % (para_name, para)
        else:
            text = "$%g$" % para
        
        if reverse:
            row = axes.shape[0] - 1 - j
        else:
            row = j
        
        bbox = axes[row, col].get_position().get_points().flatten()
        if loc == "left":
            x = bbox[0]
        else:
            x = bbox[2]
        y = (bbox[1] + bbox[3]) * 0.5
        fig.text(x, y, text, fontsize=fontsize, ha=ha, va="center",
            rotation=rot)
        

def PD_pA_pB(Lx, J, rho0, Dr=0.1, eta=-2, Ly=None):
    if Ly is None:
        Ly = Lx
    if Lx == Ly == 40:
        if rho0 == 40:
            pA_arr = np.arange(1, 8) * 10
            pB_arr = np.arange(1, 8) * 10
        elif rho0 == 20:
            pA_arr = np.arange(1, 8) * 5
            pB_arr = np.arange(1, 8) * 5
    elif Lx == Ly == 80:
        if J == 0.3:
            if rho0 == 20:
                pA_arr = np.arange(2, 7) * 5
                pB_arr = np.arange(2, 7) * 5
        else:
            if rho0 == 40:
                pA_arr = np.array([20, 25, 30])
                pB_arr = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60])
            elif rho0 == 20:
                pA_arr = np.array([20, 25, 30]) / 2
                pB_arr = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60]) / 2
    
    prefix = "/scratch03.local/yduan/QS5"
    if Lx == Ly:
        folder = f"{prefix}/L{Lx}_k0.7_r{rho0:g}_e{eta:g}_J{J:g}"
    if not os.path.exists(folder):
        folder = f"{prefix}/L{Lx}_k0.7_r{rho0:g}_e{eta:g}_J{J:.1f}"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_" \
        f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_1002.gsd"
    fout_name = f"fig/PD/pA_pB/L{Lx}_{Ly}_Dr{Dr:g}_r{rho0:g}_e{eta:g}_J{J:.2f}.png"
    nrows = pB_arr.size
    ncols = pA_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, pB in enumerate(pB_arr[::-1]):
        for i, pA in enumerate(pA_arr):
            fin = file_pat % (pA, pB)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i])
            axes[j, i].axis("off")

    fontsize = 30
    if nrows == ncols == 7:
        rect = [0.02, 0.02, 1.005, 1.005]  # left, bottom, right, top
    elif nrows == 9 and ncols == 3:
        rect = [0.015, 0.015, 1.025, 1.008]
    else:
        rect = [0.015, 0.015, 1.01, 1.01]
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    add_xlabel(fig, axes, "bottom", r"\phi_A", pA_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\phi_B", pB_arr, fontsize,
        vertical=True, reverse=True)

    if fout_name is None:
        plt.show()
    else:
        plt.savefig(fout_name)
    plt.close()


def PD_JAB_JBA(Lx, rho0, Dr, phiA=None, phiB=None, eta=-2, Ly=None):
    if phiA is None:
        phiA = rho0
    if phiB is None:
        phiB = rho0
    if Ly is None:
        Ly = Lx
    if Lx == Ly == 40:
        JAB_arr = np.array([0, 0.04, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.2, 2.8])
        JBA_arr = np.array([0, -0.04, -0.1, -0.2, -0.4, -0.8, -1.2, -1.6, -2.2, -2.8])
        if rho0 == 20:
            seed = 1234
    elif Lx == 20 and Ly == 5:
        # JAB_arr = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
        JAB_arr = np.array([1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
        JBA_arr = np.array([0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6,
            -1.8, -2.0, -2.2, -2.4, -2.6, -2.8, -3.0])
        seed = 1315
    elif Lx == 40 and Ly == 5:
        # JAB_arr = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
        JAB_arr = np.array([1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
        JBA_arr = np.array([0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6,
            -1.8, -2.0, -2.2, -2.4, -2.6, -2.8, -3.0])
        seed = 1314
    prefix = "/scratch03.local/yduan/QS5"
    if Lx == Ly:
        folder = f"{prefix}/L{Lx}_k0.7_p{rho0:g}_e{eta:g}"
    else:
        folder = f"{prefix}/L{Lx}_{Ly}_k0.7_JAB"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k0.70_p{phiA:g}_{phiB:g}_" \
        f"r{rho0:g}_e{eta:.3f}_J%.3f_%.3f_{seed}.gsd"
    fout_name = f"fig/PD/JAB_JBA/L{Lx}_{Ly}_Dr{Dr:g}_r{rho0:g}_e{eta:g}.png"
    ncols = JAB_arr.size
    nrows = JBA_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 * Ly / Lx + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, JBA in enumerate(JBA_arr):
        for i, JAB in enumerate(JAB_arr):
            fin = file_pat % (JAB, JBA)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i])
            axes[j, i].axis("off")
    
    if Lx == Ly:
        fontsize = 30
        rect = [0.02, -0.005, 1.005, 0.98]  # left, bottom, right, top
    else:
        fontsize = 15
        rect = [0.04, -0.015, 1.01, 0.98]
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    add_xlabel(fig, axes, "top", r"\eta_{AB}", JAB_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\eta_{BA}", JBA_arr, fontsize, Lx == Ly, False)

    if fout_name is None:
        plt.show()
    else:
        plt.savefig(fout_name)
    plt.close()


def PD_eta_J(Lx=40, Ly=None, k=0.7, rho0=40, Dr=3, pA=None, pB=None, seed=1200):
    if Ly is None:
        Ly = Lx
    if pA is None:
        pA = rho0
    if pB is None:
        pB = rho0
    
    eta_arr = np.array([-0.9, -1.0, -1.1, -1.4, -2.0, -3.0])
    J_arr = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/L{Lx}_{Ly}_k{k:g}_p{rho0:g}_beta"
    print(folder)
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
        f"r{rho0:g}_e%.3f_J%.3f_%.3f_{seed}.gsd"
    fout_name = f"fig/PD/eta_etaAB/L{Lx}_{Ly}_Dr{Dr:g}_r{rho0:g}_s{seed}.png"

    nrows = J_arr.size
    ncols = eta_arr.size
    if Lx == Ly:
        figsize = (ncols * 2 + 0.25, nrows * 2 * Ly / Lx + 0.25)
    else:
        figsize = (ncols * 2 + 0.5, nrows * 2 * Ly / Lx + 0.25)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, J in enumerate(J_arr):
        for i, eta in enumerate(eta_arr):
            fin = file_pat % (eta, J, -J)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i], ms=0.1)
            axes[j, i].axis("off")
    
    if Lx == 40 and Ly == 5:
        fontsize = 15
        rect = [0.06, -0.02, 1.01, 0.99] # left, bottom, right, top
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    add_xlabel(fig, axes, "top", r"\eta", eta_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\eta_{AB}", J_arr, fontsize, Lx == Ly, False)

    plt.savefig(fout_name)
    plt.close()


def PD_J_seed(Lx=20, Ly=None, k=0.7, rho0=40, pA=None, pB=None, eta=0, Dr=0.1):
    if pA is None:
        pA = rho0
    if pB is None:
        pB = rho0
    if Ly is None:
        Ly = Lx
    
    if Lx == Ly ==20 and k == 0.7 and Dr == 0.1:
        if rho0 == pA == pB == 40 or rho0 == pA == pB == 80:
            J_arr = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
            seed_arr = np.arange(1001, 1011)

    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/CS_L{Lx}_k0.7_beta"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{pA:g}_{pB:g}_" \
        f"r{rho0:g}_e{eta:.3f}_J%.3f_%.3f_%d.gsd"
    fout_name = f"fig/PD/CS/L{Lx}_Dr{Dr:g}_r{rho0:g}_e{eta:g}.png"

    nrows = J_arr.size
    ncols = seed_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, J in enumerate(J_arr):
        for i, seed in enumerate(seed_arr):
            fin = file_pat % (J, -J, seed)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i], ms=0.1)
            axes[j, i].axis("off")

    fontsize = 30
    rect = [0.015, -0.005, 1.01, 0.985] # left, bottom, right, top
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    add_xlabel(fig, axes, "top", r"replica", seed_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\eta_{AB}", J_arr, fontsize, Lx == Ly, False)


    plt.savefig(fout_name)
    plt.close()


def PD_rho0_seed(Lx=20, Ly=None, k=0.7, eta=0, alpha=0.8, Dr=0.1):
    if Ly is None:
        Ly = Lx

    if Lx == Ly == 20 and k == 0.7 and Dr == 0.1:
        rho0_arr = np.array([20, 40, 60, 80, 100])
        seed_arr = np.arange(1001, 1011)
    
    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/CS_L{Lx}_k0.7_beta"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p%g_%g_" \
        f"r%g_e{eta:.3f}_J{alpha:.3f}_{-alpha:.3f}_%d.gsd"
    fout_name = f"fig/PD/CS/L{Lx}_Dr{Dr:g}_J{alpha:.2f}_e{eta:g}.png"

    nrows = rho0_arr.size
    ncols = seed_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    ms_dict = {20: 0.25, 40: 0.2, 60: 0.15, 80: 0.1, 100: 0.1}
    for j, rho0 in enumerate(rho0_arr):
        for i, seed in enumerate(seed_arr):
            fin = file_pat % (rho0, rho0, rho0, seed)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i], ms=ms_dict[rho0])
            axes[j, i].axis("off")
    
    fontsize = 25
    rect = [0.015, -0.003, 1.01, 0.98] # left, bottom, right, top
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    add_xlabel(fig, axes, "top", r"replica", seed_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\rho_0", rho0_arr, fontsize, Lx == Ly, False)

    plt.savefig(fout_name)
    plt.close()


def PD_eta_seed(Lx=20, Ly=None, k=0.7, alpha=0.8, Dr=0.1, rho0=80):
    if Ly is None:
        Ly = Lx
    if Lx == Ly == 20 and k == 0.7 and Dr == 0.1 and alpha == 0.8:
        eta_arr = np.array([0, -0.1, -0.2, -0.4, -0.6, -0.8, -1.0,
            -1.2, -1.4, -1.6, -1.8, -2.0])
        seed_arr = np.arange(1001, 1011)

    prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/CS_L{Lx}_k0.7_beta"
    file_pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_p{rho0:g}_{rho0:g}" \
        f"_r{rho0:g}_e%.3f_J{alpha:.3f}_{-alpha:.3f}_%d.gsd"
    fout_name = f"fig/PD/CS/L{Lx}_{Ly}_Dr{Dr:g}_J{alpha:.2f}_r{rho0:g}.png"

    nrows = eta_arr.size
    ncols = seed_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eta in enumerate(eta_arr):
        for i, seed in enumerate(seed_arr):
            fin = file_pat % (eta, seed)
            print(fin)
            plot_one_panel(fin, ax=axes[j, i], ms=0.1)
            axes[j, i].axis("off")

    rect = [0.015, -0.003, 1.01, 0.98] # left, bottom, right, top
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)

    fontsize = 25
    add_xlabel(fig, axes, "top", r"replica", seed_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\eta", eta_arr, fontsize, Lx == Ly, False)

    plt.savefig(fout_name)
    plt.close()


if __name__ == "__main__":
    # fin = r"D:\tmp\L40_k0.7_r40_e-2_J0.5\L40_40_Dr0.100_k0.70_p30_40_r40_e-2.000_J0.500_-0.500_1002.gsd"
    # plot_one_panel(fin)
    
    # PD_pA_pB(Lx=80, J=0.5, rho0=20, eta=-2, Dr=0.1)

    # PD_eta_J(Ly=5)

    # PD_JAB_JBA(Lx=40, rho0=20, Dr=0.1, eta=-2, Ly=40)

    ### CS
    PD_J_seed(rho0=80)
    # PD_rho0_seed()
    # PD_eta_seed()
