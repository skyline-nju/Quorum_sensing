import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("../")

from add_line import add_line
from center_traj import read_traj, untangle
from cluster_tracking import read_one_file, detect_cluster, cal_dis2


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    f_single = "../MSD/data/center_traj/mc_L1280_1280_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.dat"

    ncut = 500
    t, x, y = read_traj(f_single, 1280, 1280, ncut)
    n = 1000
    t_new = np.arange(1, n+1) * 100

    rr_new = np.zeros(n)

    count = 0
    for j in range(t.size - n):
        x0, y0 = x[j], y[j]
        rr_new += (x[j+1:j+n+1] - x0)**2 + (y[j+1: j+n+1] - y0)**2
        count += 1
    rr_new /= count
    ax.plot(t_new, rr_new/rr_new[0], lw=2, label="single cluster, time-averaged")
    ax.set_xscale("log")
    ax.set_yscale("log")


    folder = "D:/tmp/data"
    f_mult = f"{folder}/L1280_1280_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2070_t300_0.bin"

    beg_frame = 0
    # beg_frame = 9849
    end_frame = 1092

    n_frames = end_frame - beg_frame
    nc = 9
    xc_arr = np.zeros((n_frames, nc))
    yc_arr = np.zeros((n_frames, nc))
    r_arr = np.zeros((n_frames, nc))
    frames = read_one_file(f_mult, 1280, 1280, 4, 300, beg_frame, end_frame)
    nx = 1280 // 4
    ny = 1280 // 4
    n_pad = 60
    for i_frame, (t, rhoA, rhoB) in enumerate(frames):
        flag_show = False
        xc, yc, r = detect_cluster(rhoB, n_pad, sigma=4, flag_show=flag_show, rhoB_thresh=13)
        if i_frame == 0:
            xc_arr[i_frame] = np.array(xc)
            yc_arr[i_frame] = np.array(yc)
        else:
            ic_left_list = [i for i in range(nc)]
            for j in range(nc):
                xc_left = np.array([xc[i] for i in ic_left_list])
                yc_left = np.array([yc[i] for i in ic_left_list])
                r_left = np.array([r[i] for i in ic_left_list])
                dis2 = cal_dis2(xc_arr[i_frame-1, j], yc_arr[i_frame-1, j], xc_left, yc_left, nx, ny)
                ic_min = np.argmin(dis2)
                xc_arr[i_frame, j] = xc_left[ic_min]
                yc_arr[i_frame, j] = yc_left[ic_min]
                r_arr[i_frame, j] = r_left[ic_min]
                del ic_left_list[ic_min]

    t = np.arange(1, xc_arr.shape[0]) * 100
    MSD = np.zeros(t.size)
    for j in range(MSD.size):
        dx2 = (xc_arr[j+1] - xc_arr[0]) ** 2 
        dy2 = (yc_arr[j+1] - yc_arr[0]) ** 2
        dis2 = dx2 + dy2
        MSD[j] = np.mean(dis2) * 16
    ax.plot(t, MSD/MSD[0], label="multi clusters")

    ax.legend(fontsize="large")
    ax.set_xlabel(r"$t$", fontsize="x-large")
    ax.set_ylabel("MSD (shifted)", fontsize="x-large")
    add_line(ax, 0, 0, 1, slope=0.5, label="0.5", xl=0.75)
    add_line(ax, 0, 0, 1, slope=1, label="1", xl=0.5)

    plt.show()
    plt.close()
