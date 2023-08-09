import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from gsd import fl

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_cross_point(x, y, yc):
    for i in range(x.size):
        if y[i+1] > yc >= y[i]:
            xc = (yc - y[i]) / (y[i+1]-y[i]) * (x[i+1]-x[i]) + x[i]
            break
    return xc


def get_r_center(r_edges):
    r = np.zeros(r_edges.size - 1)
    r[0] = 0
    r[1:] = (r_edges[1:r_edges.size-1] + r_edges[2:])/2
    return r


def cal_rho_r(i_frame: int,
              fgsd: fl.GSDFile,
              dx: float,
              Lx: float,
              Ly: float):
    n = int(0.5 * Lx/dx)
    r_edges = np.arange(n+1) * dx
    rho = np.zeros(n)
    

    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    # try:
    #     id_type = fgsd.read_chunk(frame=i_frame, name="particles/typeid")
    # except KeyError:
    #     id_type = np.zeros(pos.shape[0], np.int32)
    # try:
    #     speed = fgsd.read_chunk(frame=i_frame, name="particles/charge")
    # except KeyError:
    #     print("Warning, not found speed[%d]" % i_frame)
    #     speed = np.zeros(id_type.size, np.single)

    x = pos[:,0]
    y = pos[:,1]
    r2 = x**2 + y**2

    for j in range(n):
        rmin2 = r_edges[j] ** 2
        rmax2 = r_edges[j+1] ** 2
        mask = np.logical_and(r2>=rmin2, r2<rmax2)
        rho[j] = np.sum(mask) / (np.pi * (rmax2 - rmin2))
    return r_edges, rho


def coarse_grain(fin: str, fout: str, dx: float, dt: float):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            rho_r_t0 = data["rho_r_t"]
            r_edges = data["r_edges"]
            existed_frames = t0.size
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            print(nframes, "frames found for", fin)
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]

            t = np.zeros(nframes, int)
            n = int(0.5 * Lx/dx)
            rho_r_t = np.zeros((nframes, n))
            if existed_frames > 0:
                t[:existed_frames] = t0
                rho_r_t[:existed_frames] = rho_r_t0

            for i_frame in range(existed_frames, nframes):
                try:
                    t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0] * dt
                except KeyError:
                    t[i_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                r_edges, rho_r_t[i_frame] = cal_rho_r(i_frame, fgsd, dx, Lx, Ly)

            np.savez_compressed(fout, t=t, r_edges=r_edges, rho_r_t=rho_r_t)

if __name__ == "__main__":

    Dr = 2.
    prefix = f"{root_sohrab}/BEC/L100_h0.01_SP"
    basename = f"L100_v1_M1_r100.000_1.000000_k0.010000_Dr{Dr:.6f}_h0.010_1000.gsd"
    fin = f"{prefix}/{basename}"
    # with fl.open(name=fin, mode="rb") as fgsd:
    #     r_edges, rho = cal_rho_r(900, fgsd, 0.5, 100, 100)
    #     plt.plot(get_r_center(r_edges), rho)
    #     plt.show()
    #     plt.close()
    fout_basename = basename.replace("gsd", "npz")
    fout = f"data/BEC/{fout_basename}"
    coarse_grain(fin, fout, 0.5, 0.01)

    # tc_arr = []
    # Dr_arr = [3.32, 3.33, 3.34, 3.35, 3.4, 3.5]
    # for Dr in [3.32, 3.33, 3.34, 3.35, 3.4, 3.5]:
    #     fnpz = f"data/BEC/L100_v1_M1_r100.000_1.000000_k0.010000_Dr{Dr:.6f}_h0.010_1000.npz"
    #     with np.load(fnpz, "r") as data:
    #         t = data["t"]
    #         r_edges = data["r_edges"]
    #         rho_r_t = data["rho_r_t"]
    #         tc = get_cross_point(t, rho_r_t[:, 0], 2000)
    #         # plt.axvline(tc)
    #         plt.plot(t - tc, rho_r_t[:, 0], label=r"$D_r=%g$" % Dr)
    #         tc_arr.append(tc)
    # # plt.plot(t, rho_r_t[:, 0] * 0.5 **2 + rho_r_t[:, 1] * (1**2 -0.5**2))
    # plt.legend(fontsize="large")
    # # plt.xscale("log")
    # plt.axhline(2000, linestyle="dashed", c="k")
    # plt.yscale("log")
    # plt.xlabel(r"$t-t_c$", fontsize="large")
    # plt.ylabel(r"$\rho(r=0, t)$", fontsize="large")
    # plt.tight_layout()
    # plt.show()
    # plt.close()

    # fig, ax = plt.subplots(1, 1, constrained_layout=True)
    # ax.plot(np.array(Dr_arr)-3.318, tc_arr, "-o")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_ylabel(r"$t_c$", fontsize='x-large')
    # ax.set_xlabel(r"$D_r - 3.318$", fontsize='x-large')
    
    # plt.show()

    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharey=True)
    ncut = [20, 200, 300, 300, 300, 300, 300]
    for i, Dr in enumerate([3.31]):
        fnpz = f"data/BEC/L100_v1_M1_r100.000_1.000000_k0.010000_Dr{Dr:.6f}_h0.010_1000.npz"
        with np.load(fnpz, "r") as data:
            t = data["t"]
            r_edges = data["r_edges"]
            rho_r_t = data["rho_r_t"]
            rho_r_m = np.mean(rho_r_t[ncut[i]:], axis=0)
            r = get_r_center(r_edges)
            rho_origin_t = rho_r_t[:, 0] * 0.5 **2 + rho_r_t[:, 1] * (1**2 -0.5**2)
            ax1.plot(t, rho_origin_t, label=r"$D_r=%g$" % Dr)
            ax2.plot(t, rho_origin_t, label=r"$D_r=%g$" % Dr)
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")

    ax1.axvspan(9e4, 1e5, color="tab:pink", alpha=0.5)

    ax2.set_xlim(9e4, 1e5)
    # plt.xscale("log")
    # plt.yscale("log")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\rho(r=0,t)$", fontsize="large")
    ax2.set_ylabel(r"$\rho(r=0,t)$", fontsize="large")

    plt.show()
    plt.close()


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    ncut = [100, 100, 300, 300]
    bin_size = [40, 40, 50, 50]
    for i, Dr in enumerate([2, 2.8, 3.3, 3.31]):
        fnpz = f"data/BEC/L100_v1_M1_r100.000_1.000000_k0.010000_Dr{Dr:.6f}_h0.010_1000.npz"
        with np.load(fnpz, "r") as data:
            t = data["t"]
            r_edges = data["r_edges"]
            rho_r_t = data["rho_r_t"]
            rho_origin_t = rho_r_t[:, 0] * 0.5 **2 + rho_r_t[:, 1] * (1**2 -0.5**2)
            # rho_origin_t = rho_r_t[:, 0]
            hist, edges = np.histogram(rho_origin_t[ncut[i]:], bins=bin_size[i], density=True)
            rho_m = np.mean(rho_origin_t[ncut[i]:])
            x = 0.5 * (edges[:-1] + edges[1:])
            ax1.plot(x, hist, "o", label=r"$D_r=%g$" % Dr, fillstyle="none")

            mean = np.mean(rho_origin_t[ncut[i]:])
            sigma = np.std(rho_origin_t[ncut[i]:])
            rho_rescaled = (rho_origin_t[ncut[i]:] - mean) / sigma
            hist, edges = np.histogram(rho_rescaled, bins=bin_size[i], density=True)
            x = 0.5 * (edges[:-1] + edges[1:])
            ax2.plot(x, hist * sigma, "o", label=r"$D_r=%g$" % Dr, fillstyle="none")
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")
    ax2.axvline(0, linestyle="dashed", c="k")

    # plt.xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\varrho\equiv\rho(r=0)$", fontsize="large")
    ax1.set_ylabel(r"$P(\varrho)$", fontsize="large")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$(\varrho-\bar{\varrho})/\sigma_\varrho$", fontsize="large")
    ax2.set_ylabel(r"$P(\varrho)\times \sigma_\varrho$", fontsize="large")
    plt.show()
    plt.close()