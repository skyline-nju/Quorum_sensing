import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import feature
from skimage.segmentation import active_contour
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.draw import circle_perimeter

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def plot_circle(xc, yc, radii, ax):
    theta = np.linspace(0, np.pi * 2, 100)
    x = xc + radii * np.cos(theta)
    y = yc + radii * np.sin(theta)
    ax.plot(x, y, "-r")

if __name__ == "__main__":
    folder = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333/coarse_grain_dx4"
    fin = f"{folder}/L640_640_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"

    with np.load(fin, "r") as data:
        t, x, y, fields = data["t"], data["x"], data["y"], data["fields"]
        print(fields.shape)

        rhoA = fields[:, 0]
        rhoB = fields[:, 1]

        fig, axes = plt.subplots(2, 4, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)

        i_frame = 40
        im1 = axes[0, 0].imshow(rhoA[i_frame], origin="lower")
        im2 = axes[1, 0].imshow(rhoB[i_frame], origin="lower")

        rhoA_smooth = ndi.gaussian_filter(rhoA[i_frame], 4)
        rhoB_smooth = ndi.gaussian_filter(rhoB[i_frame], 4)

        edges1 = feature.canny(rhoA_smooth)
        edges2 = feature.canny(rhoB_smooth)


        axes[0, 1].imshow(rhoA_smooth, origin="lower")
        axes[1, 1].imshow(rhoB_smooth, origin="lower")
        axes[0, 2].imshow(edges1, origin="lower")
        axes[1, 2].imshow(edges2, origin="lower")


        hough_radii = np.arange(8, 10, 22)
        hough_res = hough_circle(edges2, hough_radii)
        accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=10)

        for center_y, center_x, radius in zip(cy, cx, radii):
            plot_circle(center_x, center_y, radius, axes[1, 2])
        plt.show()
        plt.close()


