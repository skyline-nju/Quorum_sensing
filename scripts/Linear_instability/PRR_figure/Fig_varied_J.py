import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg


plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(8.2, 6.1), constrained_layout=True)

    imgs = ["fig/L40_40_Dr0.1_r20_e-2_J0.00_-1.00.jpeg",
            "fig/L40_40_Dr0.1_r20_e-2_J0.10_-0.90.jpeg",
            "fig/L40_40_Dr0.1_r20_e-2_J0.25_-0.75.jpeg",
            "fig/L40_40_Dr0.1_r20_e-2_J1.00_0.00.jpeg",
            "fig/L40_40_Dr0.1_r20_e-2_J0.90_-0.10.jpeg",
            "fig/L40_40_Dr0.1_r20_e-2_J0.75_-0.25.jpeg"]
    
    titles = [r"(a)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (0, -1),
              r"(b)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (0.1, -0.9),
              r"(c)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (0.25, -0.75),
              r"(d)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (1, 0),
              r"(e)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (0.9, -0.1),
              r"(f)$\,\eta^0_{AB}=%g, \eta^0_{BA}=%g$" % (0.75, -0.25)]

    label_fontsize = "xx-large" 
    for i, ax in enumerate(axes.flat):
        image = mpimg.imread(imgs[i])
        extent = [-0.25, 3.25, -0.25, 3.25]
        ax.imshow(image, extent=extent)
        for tick in ax.get_yticklabels():
            tick.set_rotation(90)
        ax.set_title(titles[i], fontsize="x-large")
        # if i % 3 == 0:
        #     ax.set_ylabel(r"$\bar{\rho}_B/\rho_0$", fontsize=label_fontsize)
        # if i // 3 == 1:
        #     ax.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize=label_fontsize)
        ax.text(0.75, 0.05, r"$\bar{\rho}_A/\rho_0$", fontsize=label_fontsize, transform=ax.transAxes)
        ax.text(0.02, 0.75, r"$\bar{\rho}_B/\rho_0$", fontsize=label_fontsize, transform=ax.transAxes, rotation=90)
    plt.show()
    # plt.savefig("fig/snaps_varied_J.pdf", dpi=150)
    plt.close()