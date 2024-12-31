import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import measure

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


if __name__ == "__main__":
    # folder = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333/coarse_grain_dx4"
    folder = "D:/tmp/data"
    fin = f"{folder}/L640_640_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"

    with np.load(fin, "r") as data:
        t, x, y, fields = data["t"], data["x"], data["y"], data["fields"]
        n_frames = fields.shape[0]
        ny, nx = fields.shape[2:]

        rhoA = fields[:, 0]
        rhoB = fields[:, 1]


        i_frame = 40

        for i_frame in range(40, n_frames):
            n_pad = 30
            fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
            rhoB_smooth = ndi.gaussian_filter(rhoB[i_frame], 4, mode="wrap")
            rhoB_pad = pad(rhoB_smooth, n_pad)

            ax.imshow(rhoB_pad, origin="lower")
            ax.axhline(n_pad, color="w")
            ax.axhline(n_pad + ny, color="w")
            ax.axvline(n_pad, color="w")
            ax.axvline(n_pad + nx, color="w")

            contours = measure.find_contours(rhoB_pad, 11)

            count = 0
            for contour in contours:
                x, y = contour[:, 1], contour[:, 0]
                xm, ym = find_center(x, y)
                xm_shift = xm - n_pad
                ym_shift = ym - n_pad

                if 0 <= xm_shift < nx and 0 <= ym_shift < ny:
                    ax.plot(x, y, linewidth=2)
                    area = shoelace_formula(x, y)
                    plot_circle(xm, ym, np.sqrt(area/np.pi), ax)
                    count += 1
            print("find", count, "clusters")
            plt.show()
            plt.close()


