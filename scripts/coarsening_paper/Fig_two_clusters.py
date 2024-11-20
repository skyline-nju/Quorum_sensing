import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg


if __name__ == "__main__":
    fig, axes = plt.subplots(4, 2, sharex="col", sharey="col",
                             figsize=(11, 7.5), constrained_layout=True, width_ratios=[2, 1])
    
    snaps_left = ["fig/snaps/ROR_t0.png", "fig/snaps/ROR_t3e5.png", "fig/snaps/ROR_t6e5.png", "fig/snaps/ROR_t1e6.png"]
    t_arr = [r"$t=0$", r"$t=3\times 10^5$", r"$t=6\times 10^5$", r"$t=1\times 10^6$"]

    for i, ax in enumerate(axes[:, 0]):
        im = mpimg.imread(snaps_left[i])
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(t_arr[i], fontsize="xx-large")
    
    axes[0, 0].set_title(r"(a) Reverse Ostwald Ripening", fontsize="xx-large")
    axes[0, 1].set_title(r"(b) Normal Ostwald Ripening", fontsize="xx-large")

    snaps_right = ["fig/snaps/FOR_t0.png", "fig/snaps/FOR_t6e5.png", "fig/snaps/FOR_t1.2e6.png", "fig/snaps/FOR_t1.7e6.png"]
    t_arr = [r"$t=0$", r"$t=6\times 10^5$", r"$t=1.2\times 10^6$", r"$t=1.7\times 10^6$"]

    for i, ax in enumerate(axes[:, 1]):
        im = mpimg.imread(snaps_right[i])
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(t_arr[i], fontsize="xx-large")

    # plt.show()
    plt.savefig("fig/Fig3.pdf")
    plt.close()