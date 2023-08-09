# plot phase diagram in (eta_AA, eta_BB) plane

import os
import numpy as np
import matplotlib.pyplot as plt
from plot_snap import plot_one_panel, add_xlabel, add_ylabel


def L20_v8_1(eta_AA, eta_BB, eta_AB=-0.1, eta_BA=0.1, rho0=60, Dr_A=3, Dr_B=6, ms=0.05, show_labels=True, seed=1000, h=0.025):
    ncols = eta_AA.size
    nrows = eta_BB.size

    data_dir = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan/QS_not_inc_self/L20_v8_1_J0.1"
    fig_dir = f"fig/PD/etaAA_etaBB"
    fout_name = f"{fig_dir}/L20_v8_1_Dr%g_%g_r%g_J%g_%g.jpg" % (Dr_A, Dr_B, rho0, eta_AB, eta_BA)
    file_pat = f"{data_dir}/L20_20_Dr{Dr_A:g}_{Dr_B:g}_v8_1_r{rho0:g}_e%.4f_%.4f_J{eta_AB:.6f}_{eta_BA:.6f}_h{h:.3f}_{seed}.gsd"

    figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)
    if show_labels:
        rect = [0.015, 0.015, 1.01, 1.01]
    else:
        rect = [-0.01, -0.01, 1.01, 1.01]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eBB in enumerate(eta_BB):
        for i, eAA in enumerate(eta_AA):
            fin = file_pat % (eAA, eBB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    if show_labels:
        add_xlabel(fig, axes, "bottom", r"\eta_{AA}", eta_AA, fontsize)
        add_ylabel(fig, axes, "left", r"\eta_{BB}", eta_BB, fontsize,
            vertical=True, reverse=False)

    plt.savefig(fout_name)
    plt.close()


def L20_v3_1(eta_AA, eta_BB, eta_AB=-0.1, eta_BA=0.1, rho0=40, Dr_A=6, Dr_B=6, ms=0.05, show_labels=True, seed=1000, h=0.1):
    ncols = eta_AA.size
    nrows = eta_BB.size

    data_dir = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan/QS_not_inc_self/L20_v3_1_J0.1"
    fig_dir = f"fig/PD/etaAA_etaBB"
    fout_name = f"{fig_dir}/L20_v3_1_Dr%g_%g_r%g_J%g_%g.jpg" % (Dr_A, Dr_B, rho0, eta_AB, eta_BA)
    file_pat = f"{data_dir}/L20_20_Dr{Dr_A:.3f}_r{rho0:g}_e%.3f_%.3f_J{eta_AB:.6f}_{eta_BA:.6f}_h{h:.3f}_{seed}.gsd"

    figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)
    if show_labels:
        rect = [0.015, 0.015, 1.01, 1.01]
    else:
        rect = [-0.01, -0.01, 1.01, 1.01]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eBB in enumerate(eta_BB):
        for i, eAA in enumerate(eta_AA):
            fin = file_pat % (eAA, eBB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    if show_labels:
        add_xlabel(fig, axes, "bottom", r"\eta_{AA}", eta_AA, fontsize)
        add_ylabel(fig, axes, "left", r"\eta_{BB}", eta_BB, fontsize,
            vertical=True, reverse=False)

    plt.savefig(fout_name)
    plt.close()


def L20_v1_1(eta_AA, eta_BB, eta_AB=-0.1, eta_BA=0.1, rho0=40, Dr_A=6, Dr_B=6, ms=0.05, show_labels=True, seed=1000, h=0.1):
    ncols = eta_AA.size
    nrows = eta_BB.size

    data_dir = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan/QS_not_inc_self/L20_v1_1_J0.1"
    fig_dir = f"fig/PD/etaAA_etaBB"
    fout_name = f"{fig_dir}/L20_v1_1_Dr%g_%g_r%g_J%g_%g.jpg" % (Dr_A, Dr_B, rho0, eta_AB, eta_BA)
    file_pat = f"{data_dir}/L20_20_Dr{Dr_A:.3f}_r{rho0:g}_e%.3f_%.3f_J{eta_AB:.6f}_{eta_BA:.6f}_h{h:.3f}_{seed}.gsd"

    figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)
    if show_labels:
        rect = [0.015, 0.015, 1.01, 1.01]
    else:
        rect = [-0.01, -0.01, 1.01, 1.01]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eBB in enumerate(eta_BB):
        for i, eAA in enumerate(eta_AA):
            fin = file_pat % (eAA, eBB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    if show_labels:
        add_xlabel(fig, axes, "bottom", r"\eta_{AA}", eta_AA, fontsize)
        add_ylabel(fig, axes, "left", r"\eta_{BB}", eta_BB, fontsize,
            vertical=True, reverse=False)

    plt.savefig(fout_name)
    plt.close()


def L40_v1_1(eta_AA, eta_BB, eta_AB=0.1, eta_BA=-0.1, rho0=40, Dr_A=3, Dr_B=3, ms=0.05, show_labels=True, seed=1000, h=0.1):
    ncols = eta_AA.size
    nrows = eta_BB.size

    data_dir = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan/QS_not_inc_self/L40_v1_1_J0.1"
    fig_dir = f"fig/PD/etaAA_etaBB"
    fout_name = f"{fig_dir}/L40_v1_1_Dr%g_%g_r%g_J%g_%g.jpg" % (Dr_A, Dr_B, rho0, eta_AB, eta_BA)
    file_pat = f"{data_dir}/L40_40_Dr{Dr_A:.3f}_r{rho0:g}_e%.3f_%.3f_J{eta_AB:.6f}_{eta_BA:.6f}_h{h:.3f}_{seed}.gsd"

    figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)
    if show_labels:
        rect = [0.015, 0.015, 1.01, 1.01]
    else:
        rect = [-0.01, -0.01, 1.01, 1.01]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eBB in enumerate(eta_BB):
        for i, eAA in enumerate(eta_AA):
            fin = file_pat % (eAA, eBB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    if show_labels:
        add_xlabel(fig, axes, "bottom", r"\eta_{AA}", eta_AA, fontsize)
        add_ylabel(fig, axes, "left", r"\eta_{BB}", eta_BB, fontsize,
            vertical=True, reverse=False)

    plt.savefig(fout_name)
    plt.close()



def L30_10_v1_1(eta_AA, eta_BB, eta_AB=0.1, eta_BA=-0.1, rho0=60, Dr_A=3, Dr_B=6, ms=0.05, show_labels=True, seed=1000, h=0.025):
    ncols = eta_AA.size
    nrows = eta_BB.size

    data_dir = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan/QS_not_inc_self/L30_10_v1_1_J0.1"
    fig_dir = f"fig/PD/etaAA_etaBB"
    fout_name = f"{fig_dir}/L30_10_v1_1_Dr%g_%g_r%g_J%g_%g.jpg" % (Dr_A, Dr_B, rho0, eta_AB, eta_BA)
    file_pat = f"{data_dir}/L30_10_Dr{Dr_A:.3f}_r{rho0:g}_e%.3f_%.3f_J{eta_AB:.6f}_{eta_BA:.6f}_h{h:.3f}_{seed}.gsd"

    figsize = (ncols * 2 + 2, nrows * 2/3 + 0.5)
    if show_labels:
        rect = [0.05, 0.03, 1.0, 1.01]
    else:
        rect = [-0.01, -0.01, 1.01, 1.01]

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, eBB in enumerate(eta_BB):
        for i, eAA in enumerate(eta_AA):
            fin = file_pat % (eAA, eBB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    if show_labels:
        add_xlabel(fig, axes, "bottom", r"\eta_{AA}", eta_AA, fontsize)
        add_ylabel(fig, axes, "left", r"\eta_{BB}", eta_BB, fontsize,
            vertical=False, reverse=False)

    plt.savefig(fout_name)
    plt.close()


if __name__ == "__main__":
    eta_AA = np.array([-4, -3, -2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0])
    eta_BB = np.array([0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6, -1.8, -2, -3, -4])
    L20_v8_1(eta_AA, eta_BB, rho0=40, h=0.04)

    # eta_AA = np.array([-2, -1.2, -1.1, -1, -0.989, -0.7, -0.6, -0.5, 0])
    # eta_BB = np.array([0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6, -1.8, -2, -3, -4, -6, -8])
    # L20_v8_1(eta_AA, eta_BB)


    # eta_AA = np.array([-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0])
    # eta_BB = np.array([0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -1.1, -1.2])
    # L30_10_v1_1(eta_AA, eta_BB)

    # eta_AA = np.array([-2, -1.5, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0])
    # eta_BB = np.array([0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -1.1, -1.2, -1.5, -2])
    # L40_v1_1(eta_AA, eta_BB)

    # eta_AA = np.array([-2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0])
    # eta_BB = np.array([0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, -1.2, -1.4, -1.6, -1.8, -2])
    # L20_v3_1(eta_AA, eta_BB, Dr_A=0.1, Dr_B=0.1)
    # L20_v3_1(eta_AA, eta_BB)
    # L20_v1_1(eta_AA, eta_BB)