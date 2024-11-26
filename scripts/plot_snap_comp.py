import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from read_gsd import read_one_frame
from gsd import hoomd
from plot_snap import add_xlabel, add_ylabel, plot_one_panel


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_last_frame(prefix, out_folder):
    pat = f"{prefix}/*.gsd"
    fins = glob.glob(pat)
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    for fin in fins:
        basename = os.path.basename(fin)
        fout = f"{out_folder}/{basename}"
        if os.path.exists(fout) and os.path.getmtime(fout) > os.path.getmtime(fin):
            print("skipping", fout)
        else:
            snap = read_one_frame(fin, -1)
            if snap is not None:
                with hoomd.open(name=fout, mode="wb") as f:
                    f.append(snap)


def PD_comp(phiA_arr, phiB_arr, ms=0.05, L=40, Dr=1., rho0=20., eta=-2., eta_AB=0.5, eta_BA=None, h=0.1):
    if eta_BA is None:
        eta_BA = -eta_AB
    nrows = phiB_arr.size
    ncols = phiA_arr.size
    if ncols == nrows:
        figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    elif ncols > nrows:
        figsize = (ncols * 2 + 0.25, nrows * 2 + 0.5)
    else:
        figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    # prefix = f"{root_sohrab}/IQS/L40"
    prefix = f"data/snap_L{L:d}_phiA_phiB"

    if eta_BA != 0.:
        file_pat = f"{prefix}/L{L:d}_{L:d}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_e{eta:.3f}_J{eta_AB:.3f}_{eta_BA:.3f}_h{h:.3f}_1004.gsd"
    else:
        file_pat = f"{prefix}/L{L:d}_{L:d}_Dr{Dr:.3f}_k0.70_p%g_%g_r{rho0:g}_e{eta:.3f}_J{eta_AB:.3f}_-{eta_BA:.3f}_h{h:.3f}_1000.gsd"
    print(file_pat)

    rect = [0.015, 0.015, 1.005, 1.01]
    for j, phi_B in enumerate(phiB_arr[::-1]):
        for i, phi_A in enumerate(phiA_arr):
            fin = file_pat % (phi_A, phi_B)
            # print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms)
                # print("find", fin)
            axes[j, i].axis("off")
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    fontsize = 30
    add_xlabel(fig, axes, "bottom", r"\bar{\rho}_A/\rho_0", phiA_arr/rho0, fontsize)
    add_ylabel(fig, axes, "left", r"\bar{\rho}_B/\rho_0", phiB_arr/rho0, fontsize, vertical=True, reverse=True)

    # fout_name = f"fig/PD_L{L:d}_phiA_phiB/PD_comp_Dr{Dr:.3f}_r{rho0:.3f}_e{eta:.3f}_J{eta_AB:.3f}_{eta_BA:.3f}.jpg"
    # plt.savefig(fout_name)
    plt.show()
    plt.close()



if __name__ == "__main__":
    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0.5_h0.1"
    # out_folder = "data/snap_L40_phiA_phiB"
    # get_last_frame(prefix, out_folder)

    # phiA_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70])
    # phiB_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70])

    # PD_comp(phiA_arr, phiB_arr, Dr=3)

    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0_-1_h0.1"
    # out_folder = "data/snap_L40_phiA_phiB"
    # get_last_frame(prefix, out_folder)

    # phiA_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
    # phiB_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])

    # PD_comp(phiA_arr, phiB_arr, Dr=0.1, eta=-2, eta_AB=0, eta_BA=-1)

    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0.25_-0.75_h0.1"
    # out_folder = "data/snap_L40_phiA_phiB"
    # get_last_frame(prefix, out_folder)

    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0.1_-0.9_h0.1"
    # out_folder = "data/snap_L40_phiA_phiB"
    # get_last_frame(prefix, out_folder)

    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0.1_-0.9_h0.1"

    eta_AB = 0.5
    eta_BA = -0.5
    Dr = 0.1
    rho0 = 20
    prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r20_e-2_J0.5"
    # prefix = f"{root_sohrab}/QS_not_inc_self/L40_k0.7_r{rho0:g}_e-2_J{eta_AB:g}_{eta_BA:g}_h0.1"


    out_folder = "data/snap_L40_phiA_phiB"
    get_last_frame(prefix, out_folder)

    # phiA_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
    # phiB_arr = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
    phiA_arr = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3]) * 20
    phiB_arr = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3]) * 20


    PD_comp(phiA_arr, phiB_arr, Dr=Dr, eta=-2, eta_AB=eta_AB, eta_BA=eta_BA, rho0=rho0)
