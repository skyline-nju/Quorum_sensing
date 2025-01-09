import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

if __name__ ==  "__main__":
    fig = plt.figure(figsize=(8, 4), layout="constrained")
    subfigs = fig.subfigures(1, 2, wspace=0.0001, hspace=0.0001, width_ratios=[1, 1])

    (ax1, ax2) = subfigs[0].subplots(2, 1)
    ax3 = subfigs[1].subplots()
    c_list = ["tab:cyan", "tab:red", "tab:orange"]

    files = ["../data/lever_rule/Lx480_p7.71_10.5725.npz",
             "../data/lever_rule/Lx480_p8.625_13.36.npz",
             "../data/lever_rule/Lx480_p9.53875_16.1475.npz"]
    for i, fin in enumerate(files):
        with np.load(fin, "rb") as data:
            x = data["x"]
            rhoA = data["rhoA"]
            rhoB = data["rhoB"]
            # mA = data["mxA"] / data["rhoA"]
            # mB = data["mxB"] / data["rhoB"]
            ax1.plot(x, rhoA/10, c=c_list[i])
            ax1.plot(x, rhoB/10, ":", c=c_list[i])

    files = ["../data/lever_rule/Lx960_p7.6525_10.5675.npz",
             "../data/lever_rule/Lx960_p8.625_13.36.npz",
             "../data/lever_rule/Lx960_p9.597_16.152.npz"
        ]
    for i, fin in enumerate(files):
        with np.load(fin, "rb") as data:
            x = data["x"]
            rhoA = data["rhoA"]
            rhoB = data["rhoB"]
            mA = data["mxA"] / data["rhoA"]
            mB = data["mxB"] / data["rhoB"]
            ax2.plot(x, rhoA/10, c=c_list[i])
            ax2.plot(x, rhoB/10, ":", c=c_list[i])

    ms = 5
    ax3.plot(7.71/10, 10.5725/10, "o", fillstyle="none", c="tab:cyan", ms=ms)
    ax3.plot(8.625/10, 13.36/10, "o", fillstyle="none", c="tab:red", ms=ms)
    ax3.plot(9.53875/10, 16.1475/10, "o", fillstyle="none", c="tab:orange", ms=ms)

    ax3.plot(7.6525/10, 10.5675/10, "s", fillstyle="none", c="tab:cyan", ms=ms)
    ax3.plot(8.625/10, 13.36/10, "s", fillstyle="none", c="tab:red", ms=ms)
    ax3.plot(9.597/10, 16.152/10, "s", fillstyle="none", c="tab:orange", ms=ms)

    k = (2.352 - 0.112) / (1.171 - 0.392)

    x0 = 9.597/10
    y0 = 16.152/10

    x = 1
    y = y0 + k * (x-x0)

    ax3.plot(x, y, "s")
    print(x, y)

    x= 1.05
    y = y0 + k * (x-x0)
    ax3.plot(x, y, "s")
    print(x, y)

    ax3.plot([0.377, 1.111], [0.112, 2.346], ":<", c="tab:blue", fillstyle="full", label=r"$L_x=480$")
    ax3.plot([0.392, 1.171], [0.112, 2.352], "--<", c="tab:grey", fillstyle="full", label=r"$L_x=960$")

    line1 = ax1.axvline(-100, linestyle="-", c="k", label=r"$S=A$")
    line2 = ax1.axvline(-100, linestyle=":", c="k", label=r"$S=B$")
    ax1.legend(loc="lower center", fontsize="large", borderpad=0.2)
    line1 = ax2.axvline(-100, linestyle="-", c="k", label=r"$S=A$")
    line2 = ax2.axvline(-100, linestyle=":", c="k", label=r"$S=B$")
    ax2.legend(loc="lower center", fontsize="large", borderpad=0.2)

    ax1.set_xlim(0, 480)
    ax2.set_xlim(0, 960)
    ax3.set_xlim(0.355, 1.19)

    ax1.set_title(r"(a) $L_x=480$", fontsize="x-large")
    ax2.set_title(r"(b) $L_x=960$", fontsize="x-large")
    ax2.set_xlabel(r"$x$", fontsize="x-large")
    ax1.set_ylabel(r"$\langle \rho_S\rangle_{y, t}/\rho_0$", fontsize="x-large")
    ax2.set_ylabel(r"$\langle \rho_S\rangle_{y, t}/\rho_0$", fontsize="x-large")
    ax3.set_title(r"(c)", fontsize="x-large")
    ax3.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize="x-large")
    ax3.set_ylabel(r"$\bar{\rho}_B/\rho_0$", fontsize="x-large")
    ax3.plot(100, 1, "ko", fillstyle="none", label=r"$L_x=480$")
    ax3.plot(100, 1, "ks", fillstyle="none", label=r"$L_x=960$")

    ax3.legend(fontsize="large")
    plt.show()
    # plt.savefig("fig/lever_rule.pdf")
    plt.close()
