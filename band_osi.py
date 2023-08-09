import os
import glob
import numpy as np
import matplotlib.pyplot as plt

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def locate_interfaces(x, y, y_threshold, Lx, acsending=True):
    xs = []
    for i in range(x.size):
        if (acsending and y[i-1] < y_threshold <= y[i]) or (not acsending and y[i-1] > y_threshold >= y[i]):
            x1, x2 = x[i-1], x[i]
            y1, y2 = y[i-1], y[i]
            if i == 0:
                x1 -= Lx
            x_threshold = x1 + (x2-x1) / (y2-y1) * (y_threshold - y1)
            if x_threshold < 0:
                x_threshold += Lx
            xs.append(x_threshold)
    return xs


if __name__ == "__main__":
    prefix = f"{root_sohrab}/QS_not_inc_self/L80_10_v1_1_J0.1/coarse_grain_dx1"
    basename = "L80_10_Dr6.000_r20_e-2.000_-2.000_J0.080000_-0.125000_h0.100_1000.npz"
    y_threshold = 30
    beg_frame = 199
    moving_direction = "+x"
    # basename = "L80_10_Dr6.000_r40_e-2.000_-2.000_J0.062500_-0.160000_h0.100_2000.npz"
    # y_threshold = 60
    # beg_frame = 0
    # moving_direction = "-x"

    fin = f"{prefix}/{basename}"
    with np.load(fin, "r") as data:
        t = data["t"]
        x = data["x"]

        rho_A = np.zeros((t.size, x.size))
        rho_B = np.zeros((t.size, x.size))

        fields = data["fields"]
        for i in range(t.size):
            rho_A[i] = np.mean(fields[i][0], axis=0)
            rho_B[i] = np.mean(fields[i][1], axis=0)

        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 6), sharex=True, sharey=True, constrained_layout=True)
        # ax1.imshow(rho_A, origin="lower", extent=[0, 80, 0, t[-1]], aspect="auto")
        # ax2.imshow(rho_B, origin="lower", extent=[0, 80, 0, t[-1]], aspect="auto")
        # # ax1.axis("a")

        # ax1.set_xlabel(r"$x$", fontsize="large")
        # ax2.set_xlabel(r"$x$", fontsize="large")
        # ax1.set_ylabel(r"$t$", fontsize="x-large")
        # ax1.set_title(r"(a) $\langle\rho_A(\mathbf{r}, t)\rangle_y$")
        # ax2.set_title(r"(b) $\langle\rho_B(\mathbf{r}, t)\rangle_y$")

        # plt.show()
        # plt.close()
        


        xAs = np.zeros((t.size-beg_frame, 2))
        xBs = np.zeros((t.size-beg_frame, 2))

        Lx = 80
        h = 0.1

        for i in range(beg_frame, t.size):
            # lineA, = plt.plot(x, rho_A[i])
            # lineB, = plt.plot(x, rho_B[i])

            if moving_direction == "+x":
                xA = locate_interfaces(x, rho_A[i], y_threshold=y_threshold, Lx=Lx, acsending=False)
                xB = locate_interfaces(x, rho_B[i], y_threshold=y_threshold, Lx=Lx, acsending=True)
            else:
                xA = locate_interfaces(x, rho_A[i], y_threshold=y_threshold, Lx=Lx, acsending=True)
                xB = locate_interfaces(x, rho_B[i], y_threshold=y_threshold, Lx=Lx, acsending=False)
            j = i - beg_frame
            # xAs[j, 0] = xA[0]
            # xAs[j, 1] = xA[1]
            # xBs[j, 0] = xB[0]
            # xBs[j, 1] = xB[1]
            if j == 0:
                xAs[0, 0] = xA[0]
                xAs[0, 1] = xA[1]
                xBs[0, 0] = xB[0]
                xBs[0, 1] = xB[1]
            else:
                if moving_direction == "+x":
                    dx1 = xA[0] - xAs[j-1, 0]
                    dx2 = xA[1] - xAs[j-1, 0]
                else:
                    dx1 = -xA[0] + xAs[j-1, 0]
                    dx2 = -xA[1] + xAs[j-1, 0]
                if dx1 < 0:
                    dx1 += Lx
                if dx2 < 0:
                    dx2 += Lx
                if dx1 < dx2:
                    xAs[j, 0] = xA[0]
                    xAs[j, 1] = xA[1]
                else:
                    xAs[j, 0] = xA[1]
                    xAs[j, 1] = xA[0]
                if moving_direction == "+x":
                    dx1 = xB[0] - xBs[j-1, 0]
                    dx2 = xB[1] - xBs[j-1, 0]
                else:
                    dx1 = -xB[0] + xBs[j-1, 0]
                    dx2 = -xB[1] + xBs[j-1, 0]
                if dx1 < 0:
                    dx1 += Lx
                if dx2 < 0:
                    dx2 += Lx
                if dx1 < dx2:
                    xBs[j, 0] = xB[0]
                    xBs[j, 1] = xB[1]
                else:
                    xBs[j, 0] = xB[1]
                    xBs[j, 1] = xB[0]
        # plt.plot(t[beg_frame:], xAs[:, 0], "o")
        # plt.plot(t[beg_frame:], xAs[:, 1], "o")

        if moving_direction == "+x":
            w1 = xAs[:, 0] - xBs[:, 1]
            w2 = xAs[:, 1] - xBs[:, 0]
            w1[w1<0] += Lx
            w2[w2<0] += Lx
            d1 = xBs[:, 0] - xAs[:, 0]
            d1[d1<0] += Lx
            d2 = xBs[:, 1] - xAs[:, 1]
            d2[d2<0] += Lx
        else:
            w1 = xBs[:, 0] - xAs[:, 1]
            w2 = xBs[:, 1] - xAs[:, 0]
            w1[w1<0] += Lx
            w2[w2<0] += Lx

            d1 = xAs[:, 1] - xBs[:, 1]
            d1[d1<0] += Lx
            d2 = xAs[:, 0] - xBs[:, 0]
            d2[d2<0] += Lx

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, constrained_layout=True)
        ax1.plot(t[beg_frame:] * h, w1)
        ax1.plot(t[beg_frame:] * h, w2)

        dx1 = xAs[1:, 0] - xAs[:-1, 0]
        dx2 = xAs[1:, 1] - xAs[:-1, 1]
        dx1[dx1<0] += Lx
        dx2[dx2<0] += Lx
        dt = (t[1] -t[0]) * h
        v1 = dx1 / dt
        v2 = dx2 / dt
        t_new = (t[beg_frame:-1] + t[beg_frame+1:]) * h / 2
        ax2.plot(t[beg_frame:]*h, d1)
        ax2.plot(t[beg_frame:]*h, d2)
        ax3.plot(t_new, v1)
        ax3.plot(t_new, v2)
        ax1.set_ylabel("band width", fontsize="x-large")
        ax2.set_ylabel("band spacing", fontsize="x-large")
        ax3.set_ylabel("band speed", fontsize="x-large")

        ax3.set_xlabel(r"$t$", fontsize='x-large')
        # ax3.set_xlim(1e6, 2e6)


        plt.show()
        plt.close()

        plt.plot(w1[:-1], v1, "o")
        plt.plot(w2[:-1], v2, "o")
        plt.xlabel("band width", fontsize="x-large")
        plt.ylabel("band speed", fontsize="x-large")
        plt.tight_layout()
        plt.show()
        plt.close()

        plt.plot(w2, d2)
        plt.show()
        plt.close()