# plot phase diagram in the composition plane

import os
import numpy as np
import matplotlib.pyplot as plt
from plot_snap import plot_one_panel, add_xlabel, add_ylabel


def square(ncols=None, nrows=None, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2.0, zero_drho=False, inc_self=True, pA_arr=None, pB_arr=None, subfix=None, bar_rhoA=None, bar_rhoB=None, seed=None, ms=0.05):
    if ncols is not None and nrows is not None:
        pA_arr = np.arange(1, ncols+1) * 10
        pB_arr = np.arange(1, nrows+1) * 10
    elif pA_arr is not None and pB_arr is not None:
        ncols = pA_arr.size
        nrows = pB_arr.size

    if inc_self:
        if not zero_drho:
            prefix = "/scratch03.local/yduan/QS5"
            fig_dir = f"fig/PD/pA_pB"
            if seed is None:
                seed = 1002
        else:
            prefix = "/scratch03.local/yduan/QS_zero_drho"
            fig_dir = f"fig/PD/pA_pB_zero_drho"
            if seed is None:
                seed = 1010
    else:
        if zero_drho:
            prefix = "/scratch03.local/yduan/QS_zero_drho_noself"
            fig_dir = f"fig/PD/pA_pB_zero_drho_noself"
            if seed is None:
                seed = 1011
        else:
            prefix = "/scratch03.local/yduan/QS_not_inc_self"
            fig_dir = f"fig/PD/pA_pB_noself"
            if seed is None:
                seed = 1000

    if bar_rhoA is None and bar_rhoB is None:
        fout_name = f"{fig_dir}/L{L}_{L}_Dr{Dr:g}_r{rho0:g}_e{eta:g}_J{J:.2f}"
    else:
        fout_name = f"{fig_dir}/L{L}_{L}_Dr{Dr:g}_r{rho0:g}_{bar_rhoA:g}_{bar_rhoB:g}_e{eta:g}_J{J:.2f}"

    if subfix is None:
        fout_name += ".jpeg"
    else:
        fout_name += f"_{subfix}.jpeg"

    if bar_rhoA is None and bar_rhoB is None:
        data_dir = f"{prefix}/L{L}_k0.7_r{rho0:g}_e{eta:g}_J{J:g}"
        file_pat = f"{data_dir}/L{L}_{L}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_" \
            f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_{seed}.gsd"
    else:
        data_dir = f"{prefix}/L{L}_k0.7_r{bar_rhoA:g}_{bar_rhoB:g}_e{eta:g}_J{J:g}"
        file_pat = f"{data_dir}/L{L}_{L}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_" \
            f"{bar_rhoA}_{bar_rhoB}_e{eta:.3f}_J{J:.3f}_-{J:.3f}_{seed}.gsd"
    nrows = pB_arr.size
    ncols = pA_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    rect = [0.015, 0.015, 1.01, 1.01]

    if bar_rhoA == 38 and bar_rhoB == 8:
        figsize = (ncols * 2 + 1, nrows * 2 + 1)
        rect = [0.005, 0.005, 1., 1.]
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, pB in enumerate(pB_arr[::-1]):
        for i, pA in enumerate(pA_arr):
            fin = file_pat % (pA, pB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
            axes[j, i].axis("off")
    
    fontsize = 30
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    add_xlabel(fig, axes, "bottom", r"\phi_A", pA_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\phi_B", pB_arr, fontsize,
        vertical=True, reverse=True)

    plt.savefig(fout_name)
    plt.close()


def square_lower_left(L, J, rho0, Dr=0.1, eta=-2.0, inc_self=True, zero_drho=False):
    # Close up for the lower-left corner of the phase diagram on the composition plane
    if inc_self:
        if not zero_drho:
            prefix = "/scratch03.local/yduan/QS5"
            fig_dir = f"fig/PD/pA_pB"
            seed = 1002
            pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 72, 74, 76, 78, 80])
            pB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 72, 74, 76, 78, 80])
        else:
            prefix = "/scratch03.local/yduan/QS_zero_drho"
            fig_dir = f"fig/PD/pA_pB_zero_drho"
            pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
            pB_arr = np.array([0, 4, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
            seed = 1010
    else:
        if zero_drho:
            prefix = "/scratch03.local/yduan/QS_zero_drho_noself"
            fig_dir = f"fig/PD/pA_pB_zero_drho_noself"
            pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 50, 60, 70, 80])
            pB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 50, 60, 70, 80])
            seed = 1011

    # prefix = "/scratch03.local/yduan/QS5"
    folder = f"{prefix}/L{L}_k0.7_r{rho0:g}_e{eta:g}_J{J:g}"
    file_pat = f"{folder}/L{L}_{L}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_" \
        f"e{eta:.3f}_J{J:.3f}_-{J:.3f}_{seed}.gsd"
    fout_name = f"{fig_dir}/L{L}_{L}_Dr{Dr:g}_r{rho0:g}_e{eta:g}_J{J:.2f}_fine.jpeg"
    nrows = pB_arr.size
    ncols = pA_arr.size
    figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    for j, pB in enumerate(pB_arr[::-1]):
        for i, pA in enumerate(pA_arr):
            fin = file_pat % (pA, pB)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i])
            axes[j, i].axis("off")
    
    fontsize = 30
    rect = [0.015, 0.015, 1.005, 1.005]
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    add_xlabel(fig, axes, "bottom", r"\phi_A", pA_arr, fontsize)
    add_ylabel(fig, axes, "left", r"\phi_B", pB_arr, fontsize,
        vertical=True, reverse=True)

    plt.savefig(fout_name)
    plt.close()

def update():
    square(10, 10, L=40, J=2, rho0=40, Dr=0.1, eta=-0.5)
    square(10, 10, L=40, J=0.9, rho0=40, Dr=0.1, eta=-0.5)
    square(10, 10, L=40, J=0.5, rho0=40, Dr=0.1, eta=-0.5)
    square(10, 10, L=40, J=0.1, rho0=40, Dr=0.1, eta=-0.5)
    pA_arr = np.arange(14) * 10
    pB_arr = np.arange(14) * 10
    square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2)
    square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2)



if __name__ == "__main__":
    # square(15, 15, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2, zero_drho=True)
    # update()

    # square(10, 10, L=40, J=0.5, rho0=40, Dr=0.1, eta=-0.5, zero_drho=True, inc_self=False)
    # square(13, 13, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2)
    # square_lower_left(L=40, J=0.5, rho0=40, Dr=0.1, eta=-2.0, inc_self=False, zero_drho=True)

    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=1, rho0=40, Dr=0.1, eta=-1, inc_self=False, zero_drho=True)


    # # rudabeh002

    # pA_arr = np.arange(14) * 10
    # pB_arr = np.arange(14) * 10
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2)
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2)

    # pA_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 72, 74, 76, 78, 80])
    # pB_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 72, 74, 76, 78, 80])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2, subfix="fine")
    # square(10, 10, L=40, J=2, rho0=40, Dr=0.1, eta=-0.5)
    # square(10, 10, L=40, J=0.9, rho0=40, Dr=0.1, eta=-0.5)
    # square(10, 10, L=40, J=0.5, rho0=40, Dr=0.1, eta=-0.5)
    # square(10, 10, L=40, J=0.1, rho0=40, Dr=0.1, eta=-0.5)
    # square_lower_left(L=40, J=2, rho0=40, Dr=0.1, eta=-2.0, inc_self=True, zero_drho=False)

    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    # pB_arr = np.array([20, 22, 24, 26, 28, 30, 40, 50, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, subfix="CL")

    # pA_arr = np.array([20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 50, 60, 70, 72, 74, 76, 78, 80])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, subfix="LC")

 
    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]) / 2
    # pB_arr = np.array([4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 30, 40, 50, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88
    #     ]) / 2
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=20, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1004)

    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
    # pB_arr = np.array([16, 18, 20, 22, 24, 30, 40, 50, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84
    #     #  ,86, 88, 90
    #     ])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1004)

    # pA_arr = np.array([9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
    # pB_arr = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=20, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1000, ms=0.1, bar_rhoA=19, bar_rhoB=4)
    # pA_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]) / 4
    # pB_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]) / 4
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=10, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1004, ms=0.1)
    # pA_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80]) / 4
    # pB_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80]) / 4
    # for Dr in [0.2, 0.3, 0.5, 1]:
    #     square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=10, Dr=Dr, eta=-2, inc_self=False, zero_drho=False, seed=1004, ms=0.1)


    # pA_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]) / 4
    # pB_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]) / 4
    # for eta in [-0.5, -0.75, -1]:
    #     square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=10, Dr=0.1, eta=eta, inc_self=False, zero_drho=False, seed=1004, ms=0.1)

    # pA_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80]) / 2
    # pB_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 70, 80]) / 2

    # for Dr in [0.1, 0.2, 0.3, 1, 3]:
    #     square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=20, Dr=Dr, eta=-2, inc_self=False, zero_drho=False, seed=1004, ms=0.08)
 
    pA_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]) / 2
    pB_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]) / 2
    square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=20, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1004, ms=0.08, subfix="full")
    # sohrab003

    # # eta = -2, eta_AB = 2
    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=True)

    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    # pB_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=True, subfix="CL")

    # pA_arr = np.array([10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=True, subfix="LC")

    # eta = -2, eta_AB = 0.5, rho_0=40, bar_rho_A=pA, bar_rho_B=pB
    # square_lower_left(40, 0.5, 40, Dr=0.1, eta=-2, inc_self=False, zero_drho=True)
    # square_lower_left(40, 0.5, 40, Dr=3.0, eta=-2, inc_self=False, zero_drho=True)
    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=4, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=True)

    # eta = -2, eta_AB = 1, rho_0=bar_rho_A=bar_rho_B=40
    # pA_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # pB_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=1, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, seed=1004)
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=1, rho0=40, Dr=3, eta=-2, inc_self=False, zero_drho=False, seed=1004)

    # # eta = -1, eta_AB = 1, rho_0=bar_rho_A=bar_rho_B=40
    # pA_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # pB_arr = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=1, rho0=40, Dr=0.1, eta=-1, inc_self=False, zero_drho=False, seed=1004)
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=1, rho0=40, Dr=3, eta=-1, inc_self=False, zero_drho=False, seed=1004)


    # # rudabeh003
    # pA_arr = np.array([0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80])
    # pB_arr = np.array([0, 4, 8, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, bar_rhoA=40, bar_rhoB=10)

    # pA_arr = np.array([20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80])
    # pB_arr = np.array([2, 6, 10, 14, 16, 18, 20, 22, 24, 26, 30, 34, 38, 42, 46, 50])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=3.0, eta=-2, inc_self=False, zero_drho=False, bar_rhoA=38, bar_rhoB=8)


    # pA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 50, 60])
    # pB_arr = np.array([16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, bar_rhoA=8, bar_rhoB=38)


    # pA_arr = np.array([18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82])
    # pB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 50, 60])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=0.5, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, bar_rhoA=38, bar_rhoB=8)


    # pA_arr = np.array([0, 2, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
    # pB_arr = np.array([14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70])
    # square(pA_arr=pA_arr, pB_arr=pB_arr, L=40, J=2, rho0=40, Dr=0.1, eta=-2, inc_self=False, zero_drho=False, bar_rhoA=6, bar_rhoB=30)
