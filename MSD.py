import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from gsd import fl
from add_line import add_line


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def find_center(i_frame: int, fgsd: fl.GSDFile, Lx: float, Ly: float):
    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    x = pos[:,0] / Lx * np.pi * 2
    y = pos[:,1] / Ly * np.pi * 2

    x_c = np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))) / (2 * np.pi) * Lx
    y_c = np.arctan2(np.mean(np.sin(y)), np.mean(np.cos(y))) / (2 * np.pi) * Ly
    return x_c, y_c


def find_center_series(fin: str, fout: str, beg_frame=0):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            y0 = data["y"]
            existed_frames = t0.size - beg_frame
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            print(nframes)
            print(nframes, "frames found for", fin)
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
            
            t = np.zeros(nframes-beg_frame, int)
            x = np.zeros(nframes-beg_frame)
            y = np.zeros(nframes-beg_frame)

            if existed_frames > 0:
                t[:existed_frames] = t0
                x[:existed_frames] = x0
                y[:existed_frames] = y0

            for i_frame in range(beg_frame + existed_frames, nframes):
                try:
                    t[i_frame-beg_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                except KeyError:
                    t[i_frame-beg_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes-beg_frame, t[i_frame-beg_frame]))
                x[i_frame-beg_frame], y[i_frame-beg_frame] = find_center(i_frame, fgsd, Lx, Ly)
            np.savez_compressed(fout, t=t, x=x, y=y)


def untangle(x, Lx):
    x_new = np.copy(x)
    phase_shift = 0 
    for i in range(1, x.size):
        dx = x[i] - x[i-1]
        if dx < -Lx/2:
            phase_shift += Lx
        elif dx > Lx/2:
            phase_shift -= Lx
        x_new[i] += phase_shift
    return x_new


def update_center_traj():
    for L in [160, 320, 640]:
        if L == 160:
            beg_frame = 40
            seed = 1000
        else:
            beg_frame = 10
            seed = 2000
        prefix = f"{root_tahmineh}/QS_MPI_noself/L{L:d}_r10_10_10_h0.333"
        out_folder = f"{prefix}/CM"
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
        fin = f"{prefix}/L{L:d}_{L:d}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.gsd"
        fout = f"{out_folder}/L{L:d}_{L:d}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"
        find_center_series(fin, fout, beg_frame=beg_frame)

        with np.load(fout, "r") as data:

            fig, ax= plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
            t, x, y = data["t"], data["x"], data["y"]
            # print(x[0], y[0])
            x, y = untangle(x, L), untangle(y, L)
            plt.plot(x, y)
            ax.plot(x[0], y[0], "ro")
            ax.plot(x[-1], y[-1], "ko")

            plt.show()
            plt.close()
    
    L = 320
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{L:d}_r10_10_10_h0.333"
    out_folder = f"{prefix}/CM"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    fin = f"{prefix}/L{L:d}_{L:d}_Dr0.100_k0.70_p4.5_1.53_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    fout = f"{out_folder}/L{L:d}_{L:d}_Dr0.100_k0.70_p4.5_1.53_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    find_center_series(fin, fout, beg_frame=10)

    with np.load(fout, "r") as data:

        fig, ax= plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
        t, x, y = data["t"], data["x"], data["y"]
        # x += L/2
        # x[x>L/2] -= L
        print(x[0], y[0])

        # ax.plot(x, y, "o")
        ax.plot(untangle(x, L), untangle(y, L))
        ax.plot(x[0], y[0], "ro")
        plt.show()
        plt.close()


def show_traj():
    fig, ax= plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    L_arr = [160, 320, 320]
    f_arr = [f"{root_tahmineh}/QS_MPI_noself/L{L_arr[0]}_r10_10_10_h0.333/CM/L{L_arr[0]}_{L_arr[0]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p4.5_1.53_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             ]
    
    labels = [
        r"$L=160, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=4.5,\bar{\rho}_B=1.53$",
    ]
    for i, L in enumerate(L_arr):
        with np.load(f_arr[i], "r") as data:
            t, x, y = data["t"], data["x"], data["y"]
            line, = ax.plot(untangle(x[:], L), untangle(y[:], L), label=labels[i])
            ax.plot(x[0], y[0], "o", c=line.get_c())
    ax.set_xlabel(r"$x$", fontsize="x-large")
    ax.set_ylabel(r"$y$", fontsize='x-large')
    ax.legend()
    plt.show()
    plt.close()


def show_MSD():
    fig, ax= plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    L_arr = [160, 320, 320]
    f_arr = [f"{root_tahmineh}/QS_MPI_noself/L{L_arr[0]}_r10_10_10_h0.333/CM/L{L_arr[0]}_{L_arr[0]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p4.5_1.53_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             ]
    labels = [
        r"$L=160, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=4.5,\bar{\rho}_B=1.53$",
    ]
    for i, L in enumerate(L_arr):
        with np.load(f_arr[i], "r") as data:
            t, x, y = data["t"] / 3, data["x"], data["y"]
            x, y = untangle(x[:], L), untangle(y[:], L)
            rr = (x-x[0])**2 + (y-y[0])**2 
            ax.plot(t-t[0], rr, label=labels[i])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$|\mathbf{r}(t)-\mathbf{r}(0)|^2$", fontsize="x-large")
    ax.legend()

    add_line(ax, 0.1, 0.6, 0.95, 1, label="slope=1", xl=0.60, yl=0.9, fontsize="large")
    add_line(ax, 0.1, 0.05, 0.95, 2, label="slope=2", xl=0.10, yl=0.1, fontsize="large")


    plt.show()
    plt.close()


    fig, ax= plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    L_arr = [160, 320, 640]
    f_arr = [f"{root_tahmineh}/QS_MPI_noself/L{L_arr[0]}_r10_10_10_h0.333/CM/L{L_arr[0]}_{L_arr[0]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[2]}_r10_10_10_h0.333/CM/L{L_arr[2]}_{L_arr[2]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             ]
    
    for i, L in enumerate(L_arr):
        with np.load(f_arr[i], "r") as data:
            t, x, y = data["t"] / 3, data["x"], data["y"]
            x, y = untangle(x[:], L), untangle(y[:], L)
            rr = (x-x[0])**2 + (y-y[0])**2 
            ax.plot(t-t[0], rr, label=labels[i])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$|\mathbf{r}(t)-\mathbf{r}(0)|^2$", fontsize="x-large")
    ax.legend()

    add_line(ax, 0.1, 0.6, 0.95, 1, label="slope=1", xl=0.60, yl=0.9, fontsize="large")
    add_line(ax, 0.1, 0.05, 0.95, 2, label="slope=2", xl=0.10, yl=0.1, fontsize="large")


    plt.show()
    plt.close()


def show_MSD_averaged():
    L_arr = [160, 320, 320]
    f_arr = [f"{root_tahmineh}/QS_MPI_noself/L{L_arr[0]}_r10_10_10_h0.333/CM/L{L_arr[0]}_{L_arr[0]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p4.5_1.53_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             ]

    i_fin = 0
    with np.load(f_arr[i_fin], "r") as data:
        t, x, y = data["t"] / 3, data["x"], data["y"]
        x, y = untangle(x[:], L_arr[i_fin]), untangle(y[:], L_arr[i_fin])
    
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    for n in [100, 10]:
        t_new = np.arange(1, n+1) * 10000
        rr_new = np.zeros(n)

        count = 0
        for j in range(t.size - n):
            x0, y0 = x[j], y[j]
            rr_new += (x[j+1:j+n+1] - x0)**2 + (y[j+1: j+n+1] - y0)**2
            count += 1
        rr_new /= count
        ax.plot(t_new, rr_new, lw=2, label=r"$\Delta t=%g$" % (n * 1e4))
    rr = (x-x[0])**2 + (y-y[0])**2 
    ax.plot(t-t[0], rr, label="raw", lw=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$|\mathbf{r}(t)-\mathbf{r}(0)|^2$", fontsize="x-large")
    ax.legend(fontsize="x-large")

    add_line(ax, 0.1, 0.2, 0.95, 1, label="slope=1", xl=0.7, yl=0.55, fontsize="large")

    plt.show()
    plt.close()


    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    labels = [
        r"$L=160, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=5,\bar{\rho}_B=2.96$",
        r"$L=320, \bar{\rho}_A=4.5,\bar{\rho}_B=1.53$",
    ]
    n = 10
    for i_fin in range(3):
        with np.load(f_arr[i_fin], "r") as data:
            t, x, y = data["t"] / 3, data["x"], data["y"]
            x, y = untangle(x[:], L_arr[i_fin]), untangle(y[:], L_arr[i_fin])

        t_new = np.arange(1, n+1) * 10000
        rr_new = np.zeros(n)

        count = 0
        for j in range(t.size - n):
            x0, y0 = x[j], y[j]
            rr_new += (x[j+1:j+n+1] - x0)**2 + (y[j+1: j+n+1] - y0)**2
            count += 1
        rr_new /= count

        
        ax.plot(t_new, rr_new, lw=2, label=labels[i_fin])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$|\mathbf{r}(t)-\mathbf{r}(0)|^2$", fontsize="x-large")
    ax.legend(fontsize="large")
    add_line(ax, 0.1, 0.2, 0.95, 1, label="slope=1", xl=0.60, yl=0.6, fontsize="large")

    plt.show()
    plt.close()

    L_arr = [160, 320, 640]
    f_arr = [f"{root_tahmineh}/QS_MPI_noself/L{L_arr[0]}_r10_10_10_h0.333/CM/L{L_arr[0]}_{L_arr[0]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[1]}_r10_10_10_h0.333/CM/L{L_arr[1]}_{L_arr[1]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
             f"{root_tahmineh}/QS_MPI_noself/L{L_arr[2]}_r10_10_10_h0.333/CM/L{L_arr[2]}_{L_arr[2]}_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz",
            ]

    labels = [
        r"$L=160$",
        r"$L=320$",
        r"$L=640$",
    ]
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    n = 10
    for i_fin in range(3):
        with np.load(f_arr[i_fin], "r") as data:
            t, x, y = data["t"] / 3, data["x"], data["y"]
            x, y = untangle(x[:], L_arr[i_fin]), untangle(y[:], L_arr[i_fin])

        t_new = np.arange(1, n+1) * 10000
        rr_new = np.zeros(n)

        count = 0
        for j in range(t.size - n):
            x0, y0 = x[j], y[j]
            rr_new += (x[j+1:j+n+1] - x0)**2 + (y[j+1: j+n+1] - y0)**2
            count += 1
        rr_new /= count
        ax.plot(t_new, rr_new, lw=2, label=labels[i_fin])
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel(r"$|\mathbf{r}(t)-\mathbf{r}(0)|^2$", fontsize="x-large")
    ax.legend(fontsize="large", title=r"$\bar{\rho}_A=5,\bar{\rho}_B=2.96$")

    add_line(ax, 0.05, 0.1, 0.95, 1, label="slope=1", xl=0.60, yl=0.6, fontsize="large")
    add_line(ax, 0.05, 0.3, 0.95, 1, label="slope=1", xl=0.60, yl=0.6, fontsize="large")
    add_line(ax, 0.05, 0.5, 0.95, 1, label="slope=1", xl=0.60, yl=0.6, fontsize="large")


    plt.show()
    plt.close()


if __name__ == "__main__":
    # update_center_traj()
    show_traj()
    show_MSD()

    