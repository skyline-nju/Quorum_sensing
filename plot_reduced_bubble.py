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


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)

    fin = "Reduced_bubble/output1.dat"
    t, n_par, r = read_data(fin)
    ax.plot(t, r, label="without equalising")

    fin = "Reduced_bubble/output2.dat"
    t, n_par, r = read_data(fin)
    ax.plot(t, r, label="with equalising")

    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.4, 0.4, 0.98, 1/3, label="1/3", xl=0.7, yl=0.6)
    ax.legend(fontsize="large")
    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel("mean domain size", fontsize="large")
    plt.show()
    plt.close()