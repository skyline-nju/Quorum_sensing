import numpy as np
import matplotlib.pyplot as plt
from add_line import add_line


def read_data(fin, h=0.1):
    with open(fin, "r") as f:
        lines = f.readlines()
        n = len(lines)
        t = np.zeros(n)
        n_par = np.zeros(n, int)
        r = np.zeros(n)

        for i, line in enumerate(lines):
            s = line.rstrip("\n").split("\t")
            t[i] = int(s[0]) * h
            n_par[i] = int(s[1])
            r[i] = float(s[2])
    return t, n_par, r


def show_R(ax):
    fin = "/home/yduan/code/Quorum_sensing/Reduced_bubble/output2.dat"
    t, n_par, r = read_data(fin)
    ax.plot(t, r, label=r"$D(R)\propto R^{-1}$")

    fin = "/home/yduan/code/Quorum_sensing/Reduced_bubble/output3.dat"
    t, n_par, r = read_data(fin)
    ax.plot(t, r, label=r"$D(R)\propto R^{-2}$")

    fin = "/home/yduan/code/Quorum_sensing/Reduced_bubble/output4.dat"
    t, n_par, r = read_data(fin)
    ax.plot(t, r, label=r"$D(R)\propto R^{-5.7}$")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(fontsize="x-large")
    # ax.set_xlabel(r"$t$", fontsize="x-large")
    # ax.set_ylabel("mean domain size", fontsize="x-large")

    add_line(ax, 0.4, 0.4, 0.98, 1/3, label="1/3", xl=0.7, yl=0.6)
    add_line(ax, 0.4, 0.25, 0.98, 0.25, label="1/4", xl=0.7, yl=0.35)
    add_line(ax, 0.4, 0.15, 0.98, 0.13, label="0.13", xl=0.7, yl=0.2)


if __name__ == "__main__":
    pass