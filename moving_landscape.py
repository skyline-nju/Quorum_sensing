import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from coarse_grain import coarse_grain_all_1D


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


if __name__ == "__main__":
    prefix = f"{root_sohrab}/moving_light_pat/L200_50"
    coarse_grain_all_1D(prefix, 0.2, False)

    Lx = 200
    Ly = 50
    h = 0.01
    c = 0.
    dx = 0.2
    v0 = 0.5
    v1 = 2
    prefix = f"{root_sohrab}/moving_light_pat/L200_50/profile_1D_dx0.2/"


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), constrained_layout=True)

    # clist = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    clist = [0, 0.1]
    # clist = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    for c in clist:
        fnpz = f"{prefix}/L{Lx}_{Ly}_r20_Dr0.1_Dt0_w20_v{v0:g}_{v1:g}_c{c:.4f}_h{h:.2f}_1000.npz"
        rho_x = np.zeros(int(Lx/dx))
        mx_x = np.zeros(int(Lx/dx))
        with np.load(fnpz, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]
            count = 0
            for i in range(15, t.size):
                shift = -int(t[i] * h * c / dx) + int(Lx/2/dx)
                rho = fields[i, 0]
                rho_x += np.roll(rho, shift)
                mx = fields[i, 2]
                mx_x += np.roll(mx, shift)
                count += 1
            rho_x /= count
            mx_x /= count

            color = cm.viridis((c-clist[0])/(clist[-1]-clist[0]))
            ax1.plot(x, rho_x, label=r"$%g$" % c, c=color, lw=1)
            ax2.plot(x, mx_x/rho_x, c=color, lw=1)

            mask = np.logical_and(x>100, x<=120)
            flux1 = np.sum(mx_x[mask] * v1)
            mask = np.logical_or(x<=100, x>120)
            flux2 = np.sum(mx_x[mask] * v0)
            flux = flux1 + flux2
            print(c, flux) 
    # plt.yscale("log")
    ax1.set_xlim(80, 140)
    ax2.set_xlim(80, 140)
    ax1.set_xlabel(r"$x-ct$", fontsize="large")
    ax2.set_xlabel(r"$x-ct$", fontsize="large")
    ax1.set_ylabel(r"$\langle\rho\rangle_{y,t}$", fontsize="large")
    ax2.set_ylabel(r"$\langle \mathbf{p}_x\rangle_{y, t} / \langle\rho\rangle_{y,t}$", fontsize="large")


    ax1.legend(title=r"$c=$", fontsize="large")
    ax1.set_yscale("log")
    plt.show()
    plt.close()



