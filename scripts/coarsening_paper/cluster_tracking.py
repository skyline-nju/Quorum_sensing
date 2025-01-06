import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import measure
from read_fields import read_one_file
sys.path.append("../")

from add_line import add_line

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def pad(rho, dn_pad):
    nrows, ncols = rho.shape
    nrows_new, ncols_new = nrows + 2 * dn_pad, ncols + 2 * dn_pad
    rho_new = np.zeros((nrows_new, ncols_new))
    rho_new[dn_pad: dn_pad+nrows, : dn_pad] = rho[:, ncols-dn_pad:]
    rho_new[dn_pad: dn_pad+nrows, dn_pad: dn_pad + ncols] = rho
    rho_new[dn_pad: dn_pad+nrows, dn_pad + ncols:] = rho[:, :dn_pad]
    rho_new[:dn_pad] = rho_new[nrows: dn_pad+nrows]
    rho_new[dn_pad+nrows:] = rho_new[dn_pad: dn_pad * 2]
    return rho_new


def find_center(x, y):
    xm = np.mean(x)
    ym = np.mean(y)
    return xm, ym


def shoelace_formula(x_arr, y_arr):
    n = x_arr.size
    area = 0.
    for i in range(n):
        x1, y1 = x_arr[i], y_arr[i]
        x2, y2 = x_arr[(i+1) % n], y_arr[(i+1) % n]
        area += x1 * y2 - x2 * y1
    return abs(area) / 2


def plot_circle(xc, yc, radii, ax):
    theta = np.linspace(0, np.pi * 2, 100)
    x = xc + radii * np.cos(theta)
    y = yc + radii * np.sin(theta)
    ax.plot(x, y, ":r")


def detect_cluster(rhoB, n_pad, sigma=4, flag_show=False, rhoB_thresh=11):
    ny, nx = rhoB.shape
    rhoB_smooth = ndi.gaussian_filter(rhoB, sigma, mode="wrap")
    rhoB_pad = pad(rhoB_smooth, n_pad)
    
    if flag_show:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
        ax.imshow(rhoB_pad, origin="lower")
        ax.axhline(n_pad, color="w")
        ax.axhline(n_pad + ny, color="w")
        ax.axvline(n_pad, color="w")
        ax.axvline(n_pad + nx, color="w")

    contours = measure.find_contours(rhoB_pad, rhoB_thresh)
    count = 0
    xc_arr, yc_arr, r_arr = [], [], []
    for contour in contours:
        x, y = contour[:, 1], contour[:, 0]
        xm, ym = find_center(x, y)
        xm_shift = xm - n_pad
        ym_shift = ym - n_pad

        if 0 <= xm_shift < nx and 0 <= ym_shift < ny:
            area = shoelace_formula(x, y)
            radius = np.sqrt(area/np.pi)
            if flag_show:
                ax.plot(x, y, linewidth=2)
                plot_circle(xm, ym, radius, ax)
            count += 1
            xc_arr.append(xm_shift)
            yc_arr.append(ym_shift)
            r_arr.append(radius)
    if flag_show:
        plt.show()
        plt.close()
    return xc_arr, yc_arr, r_arr


def cal_dis2(x0, y0, x_arr, y_arr, nx, ny):
    half_nx = nx / 2
    half_ny = ny / 2
    dx_arr = x_arr - x0
    dy_arr = y_arr - y0
    dx_arr[dx_arr < -half_nx] += nx
    dx_arr[dx_arr > half_nx] -= nx
    dy_arr[dy_arr < -half_ny] += ny
    dy_arr[dy_arr > half_ny] -= ny
    dis2 = dx_arr**2 + dy_arr**2
    return dis2


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


def untangle_multi_traj(xc_arr, yc_arr, Lx, Ly):
    xc_arr_new = np.zeros_like(xc_arr)
    yc_arr_new = np.zeros_like(yc_arr)
    for i in range(xc_arr.shape[1]):
        xc_arr_new[:, i] = untangle(xc_arr[:, i], Lx)
        yc_arr_new[:, i] = untangle(yc_arr[:, i], Ly)
    return xc_arr_new, yc_arr_new


def L1280():
    # folder = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333/coarse_grain_dx4"
    folder = "D:/tmp/data"
    # fin = f"{folder}/L640_640_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    fin = f"{folder}/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2040_t300_0.bin"

    # beg_frame = 8638
    # end_frame = 9849
    # frames = read_one_file(fin, 1280, 1280, 4, 300, beg_frame, end_frame)
    # n_pad = 30
    # flag_nc_increase = False
    # nc_pre = None
    # for (t, rhoA, rhoB) in frames:
    #     flag_show = flag_nc_increase
    #     xc, yc, r = detect_cluster(rhoB, n_pad, sigma=4, flag_show=flag_show, rhoB_thresh=13)
    #     nc_cur = len(xc)
    #     if nc_pre is not None and nc_cur > nc_pre:
    #         flag_nc_increase = True
    #     else:
    #         flag_nc_increase = False
    #         if nc_pre is None or nc_cur < nc_pre:
    #             print(t//300, nc_cur)
    #     nc_pre = nc_cur
    
    beg_frame = 8639
    # beg_frame = 9849
    end_frame = 9847

    n_frames = end_frame - beg_frame
    nc = 32
    xc_arr = np.zeros((n_frames, nc))
    yc_arr = np.zeros((n_frames, nc))
    r_arr = np.zeros((n_frames, nc))
    frames = read_one_file(fin, 1280, 1280, 4, 300, beg_frame, end_frame)
    nx = 1280 // 4
    ny = 1280 // 4
    n_pad = 30
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

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    xc_arr, yc_arr = untangle_multi_traj(xc_arr, yc_arr, nx, ny)
    for j in range(nc):
        xc = xc_arr[:, j]
        yc = yc_arr[:, j]
        # xc = untangle(xc, nx)
        # yc = untangle(yc, ny)
        ax.plot(xc, yc)
        ax.plot([xc[0], xc[-1]], [yc[0], yc[-1]], lw=2)
    
    ax.axis("equal")
    plt.show()
    plt.close()


    fig, ax = plt.subplots(1, 1)
    t = np.arange(1, xc.shape[0]) * 100
    MSD = np.zeros(xc.shape[0]-1)
    for j in range(MSD.size):
        dx2 = (xc_arr[j+1] - xc_arr[0]) ** 2 
        dy2 = (yc_arr[j+1] - yc_arr[0]) ** 2
        dis2 = dx2 + dy2
        # r_mean = np.mean(r_arr[:j+1], axis=0)
        # dis2 *= r_mean**2
        MSD[j] = np.mean(dis2)
    ax.plot(t, MSD)
    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0, 0, 0.8, slope=0.5)
    add_line(ax, 0, 0, 0.8, slope=1)
    plt.show()
    plt.close()


if __name__ == "__main__":
    folder = "D:/tmp/data"
    fin = f"{folder}/L2560_2560_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2079_t300_0.bin"

    # beg_frame = 0
    # end_frame = None
    # frames = read_one_file(fin, 2560, 2560, 4, 300, beg_frame, end_frame)
    # n_pad = 30
    # flag_nc_increase = False
    # nc_pre = None
    # for (t, rhoA, rhoB) in frames:
    #     flag_show = flag_nc_increase
    #     xc, yc, r = detect_cluster(rhoB, n_pad, sigma=4, flag_show=flag_show, rhoB_thresh=13)
    #     nc_cur = len(xc)
    #     if nc_pre is not None and nc_cur > nc_pre:
    #         flag_nc_increase = True
    #     else:
    #         flag_nc_increase = False
    #         if nc_pre is None or nc_cur < nc_pre:
    #             print(t//300, nc_cur)
    #     nc_pre = nc_cur
    
    beg_frame = 100
    # # beg_frame = 9849
    end_frame = 411

    n_frames = end_frame - beg_frame
    nc = 77
    xc_arr = np.zeros((n_frames, nc))
    yc_arr = np.zeros((n_frames, nc))
    r_arr = np.zeros((n_frames, nc))
    frames = read_one_file(fin, 2560, 2560, 4, 300, beg_frame, end_frame)
    nx = 2560 // 4
    ny = 2560 // 4
    n_pad = 30
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

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    xc_arr, yc_arr = untangle_multi_traj(xc_arr, yc_arr, nx, ny)
    for j in range(nc):
        xc = xc_arr[:, j]
        yc = yc_arr[:, j]
        # xc = untangle(xc, nx)
        # yc = untangle(yc, ny)
        ax.plot(xc, yc)
        ax.plot([xc[0], xc[-1]], [yc[0], yc[-1]], lw=2)
    
    ax.axis("equal")
    plt.show()
    plt.close()


    fig, ax = plt.subplots(1, 1)
    t = np.arange(1, xc.shape[0]) * 100
    MSD = np.zeros(xc.shape[0]-1)
    for j in range(MSD.size):
        dx2 = (xc_arr[j+1] - xc_arr[0]) ** 2 
        dy2 = (yc_arr[j+1] - yc_arr[0]) ** 2
        dis2 = dx2 + dy2
        # r_mean = np.mean(r_arr[:j+1], axis=0)
        # dis2 *= r_mean**2
        MSD[j] = np.mean(dis2)
    ax.plot(t, MSD)
    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0, 0, 0.8, slope=0.5)
    add_line(ax, 0, 0, 0.8, slope=1)
    plt.show()
    plt.close()
    

         


