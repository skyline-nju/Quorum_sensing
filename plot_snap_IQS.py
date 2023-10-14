import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from read_gsd import read_one_frame
from gsd import hoomd
from plot_snap import add_xlabel, add_ylabel


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def plot_one_panel(fname, i_frame=-1, ax=None, ms=0.05, color_coding="none", frac=1.):
    snap = read_one_frame(fname, i_frame)
    if snap is not None:
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
            flag_show = True
        else:
            flag_show = False
        
        Lx = snap.configuration.box[0]
        Ly = snap.configuration.box[1]
        x = snap.particles.position[:, 0] + Lx / 2
        y = snap.particles.position[:, 1] + Ly / 2
        theta = -snap.particles.position[:,2]
        if frac < 1.:
            n = int(x.size * frac)
            x, y, theta = x[:n], y[:n], theta[:n]
        if color_coding == "none":
            ax.plot(x, y, ".", ms=ms, c="b", alpha=1)
        elif color_coding == "ori":
            ax.scatter(x, y, s=ms, c=theta, cmap="hsv")

        ax.axis("scaled")
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)

        if flag_show:
            plt.show()
            plt.close()
        

def PD_gamma_Dr(gamma_arr, Dr_arr, ms=0.05, frac=1., L=40, eta=6.):
    nrows = Dr_arr.size
    ncols = gamma_arr.size
    if ncols == nrows:
        figsize = (ncols * 2 + 0.25, nrows * 2 + 0.25)
    elif ncols > nrows:
        figsize = (ncols * 2 + 0.25, nrows * 2 + 0.5)
    else:
        figsize = (ncols * 2 + 0.5, nrows * 2 + 0.25)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)

    # prefix = f"{root_sohrab}/IQS/L40"
    prefix = f"data/snap_IQS/L{L:d}"
    file_pat = f"{prefix}/L{L:d}_{L:d}_Dr%.3f_Dt0_r20_p20_e{eta:.3f}_g%.3f_h0.020_1000.gsd"

    rect = [0.015, 0.015, 1.005, 1.01]
    for j, Dr in enumerate(Dr_arr[::-1]):
        for i, gamma in enumerate(gamma_arr):
            fin = file_pat % (Dr, gamma)
            print(fin)
            if os.path.exists(fin):
                plot_one_panel(fin, ax=axes[j, i], ms=ms, color_coding="ori", frac=frac)
                # print("find", fin)
            axes[j, i].axis("off")
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    fontsize = 30
    add_xlabel(fig, axes, "bottom", r"\gamma_v", gamma_arr, fontsize)
    add_ylabel(fig, axes, "left", r"D_r", Dr_arr, fontsize, vertical=True, reverse=True)

    fout_name = f"fig/IQS/L{L:d}_e{eta:g}.jpg"
    plt.savefig(fout_name)
    # plt.show()
    plt.close()



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
            with hoomd.open(name=fout, mode="wb") as f:
                f.append(snap)


if __name__ == "__main__":
    prefix = f"{root_sohrab}/IQS/L40"
    out_folder = "data/snap_IQS/L40"
    get_last_frame(prefix, out_folder)


    Dr_arr = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1])
    gamma_arr = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20])

    PD_gamma_Dr(gamma_arr, Dr_arr, ms=0.025, frac=0.5, eta=4)
    PD_gamma_Dr(gamma_arr, Dr_arr, ms=0.025, frac=0.5, eta=2)



