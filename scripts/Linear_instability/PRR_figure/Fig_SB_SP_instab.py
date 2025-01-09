import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


def SB_instability():
    fig, axes = plt.subplots(4, 3, sharex="col", sharey=True, figsize=(9, 3.2), constrained_layout=True, width_ratios=[1, 2, 4])

    Lx_arr = [80, 160, 320]
    t_labels = [r"$t=0$", r"$t=4000$", r"$t=5000$", r"$t=10^4$"]
    labels = ["(a)", "(b)", "(c)"]
    for row in range(4):
        for col in range(3):
            fname = "fig/SB_instab/%d_%d.png" % (Lx_arr[col], row)
            im = mpimg.imread(fname)
            ax = axes[row, col]
            ax.imshow(im)
            ax.set_xticks([])
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(t_labels[row])
            if row == 0:
                ax.set_title(r"%s $L_x=%d$" % (labels[col], Lx_arr[col]))

    # plt.show()
    plt.savefig("fig/SB_instability.pdf")
    plt.close()


def SP_instability():
    fig, axes = plt.subplots(2, 6, sharex="row", sharey="row", figsize=(12, 4.8), gridspec_kw=dict(hspace=0.25, wspace=0, left=0., right=1, bottom=0.0, top=0.94))

    t_arr = [160, 200, 300, 400, 1000, 20000]
    titles = [r"(a) $t=160$", r"(b) $t=200$", r"(c) $t=300$", r"(d) $t=400$", r"(e) $t=10^3$", r"(f) $t=2\times 10^4$"]

    for i, ax in enumerate(axes[0]):
        fname = "fig/spiral_instab/%d.png" % t_arr[i]
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[i], fontsize="x-large")
    

    t_arr = [0, 500, 1000, 2020, 3000, 5000]
    titles = [r"(g) $t=0$", r"(h) $t=500$", r"(i) $t=1000$", r"(j) $t=2020$", r"(k) $t=3000$", r"(l) $t=5000$"]
    for i, ax in enumerate(axes[1]):
        fname = "fig/spiral_instab2/%d.png" % t_arr[i]
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[i], fontsize="x-large")
    plt.show()
    # plt.savefig("fig/SP_instability.pdf")
    plt.close()


if __name__ == "__main__":
    SP_instability()