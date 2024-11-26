import os
import numpy as np
import matplotlib.pyplot as plt
from gsd import fl
import sys


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"
root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"


def find_xc(i_frame: int, fgsd: fl.GSDFile, Lx: float, which="both"):
    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    if which == "both":
        x = pos[:,0] / Lx * np.pi * 2
    else:
        type_id = fgsd.read_chunk(frame=i_frame, name="particles/typeid")
        if which == "A":
            mask = type_id == 0
        elif which == "B":
            mask = type_id == 1
        else:
            print("Error, which must be one of 'both', 'A' or 'B'.")
            sys.exit()
        x = pos[mask][:, 0] / Lx * np.pi * 2

    x_c = np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))) / (2 * np.pi) * Lx
    return x_c


def get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None, inc_v_mx=False):
    ncols = int(Lx / dx)
    bin_area = dx * Ly

    par_num = np.zeros((2, ncols), int)
    mx, my, v, vx, vy = np.zeros((5, 2, ncols))
    if inc_v_mx:
        v_mx = np.zeros_like(mx)

    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    try:
        id_type = fgsd.read_chunk(frame=i_frame, name="particles/typeid")
    except KeyError:
        id_type = np.zeros(pos.shape[0], np.int32)
    try:
        speed = fgsd.read_chunk(frame=i_frame, name="particles/charge")
    except KeyError:
        print("Warning, not found speed[%d]" % i_frame)
        speed = np.zeros(id_type.size, np.single)
    
    if center is None:
        x = pos[:,0] + Lx / 2
    else:
        x = pos[:, 0] - center + Lx/2
        x[x<0] += Lx
        x[x>=Lx] -= Lx

    ux = np.cos(pos[:, 2])
    uy = np.sin(pos[:, 2])

    col = (x/dx).astype(int)
    col[col < 0] += ncols
    col[col >= ncols] -= ncols
    for j in range(x.size):
        my_id = id_type[j]
        my_col = col[j]
        par_num[my_id, my_col] += 1
        mx[my_id, my_col] += ux[j]
        my[my_id, my_col] += uy[j]
        v[my_id, my_col] += speed[j]
        vx[my_id, my_col] += ux[j] * speed[j]
        vy[my_id, my_col] += uy[j] * speed[j]
    
    rho = par_num / bin_area
    mx /= bin_area
    my /= bin_area
    v /= bin_area
    vx /= bin_area
    vy /= bin_area
    if inc_v_mx:
        v_ave = v * bin_area
        mask = par_num > 0
        v_ave[mask] /= par_num[mask]
        v_mx = v_ave * mx
        return rho, mx, my, v, vx, vy, v_mx, v_ave
    else:
        return rho, mx, my, v, vx, vy


def get_L(fname):
    s = os.path.basename(fname).lstrip("L").split("_")
    Lx = int(s[0])
    Ly = int(s[1])
    return Lx, Ly


def time_ave(fname, beg_frame=20, end_frame=None, center_type="both", dx=1, inc_v_mx=False, folder="data/time_ave_profiles"):
    basename = os.path.basename(fname)
    print("handling", basename)
    fout = f"{folder}/{basename.replace('.gsd', '.npz')}"
    Lx, Ly = get_L(basename)
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            rho, mx, my, v, vx, vy = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]
            nframes = data["nframes"]
            rho *= nframes
            mx *= nframes
            my *= nframes
            v *= nframes
            vx *= nframes
            vy *= nframes
            if inc_v_mx:
                v_mx = data["v_mx"]
                v_ave = data["v_ave"]
    else:
        ncols = int(round(Lx / dx))
        rho, mx, my, v, vx, vy = np.zeros((6, 2, ncols))
        if inc_v_mx:
            v_mx = np.zeros((2, ncols))
            v_ave = np.zeros((2, ncols))
        nframes = 0
    with fl.open(name=fname, mode="rb") as fgsd:
        count = 0
        if end_frame is None:
            end_frame = fgsd.nframes
        for i in range(beg_frame + nframes, end_frame):
            if center_type is None:
                center = None
            else:
                center = find_xc(i, fgsd, Lx, which=center_type)
            if not inc_v_mx:
                rho_cur, mx_cur, my_cur, v_cur, vx_cur, vy_cur = get_profiles(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
            else:
                rho_cur, mx_cur, my_cur, v_cur, vx_cur, vy_cur, v_mx_cur, v_ave_cur = get_profiles(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center, inc_v_mx=inc_v_mx)
            rho += rho_cur
            mx += mx_cur
            my += my_cur
            v += v_cur
            vx += vx_cur
            vy += vy_cur
            if inc_v_mx:
                v_mx += v_mx_cur
                v_ave += v_ave_cur
            count +=1
            print("frame:%d/%d" % (i, fgsd.nframes))
    nframes += count
    rho /= nframes
    mx /= nframes
    my /= nframes
    v /= nframes
    vx /= nframes
    vy /= nframes

    if not inc_v_mx:
        np.savez_compressed(fout, rho=rho, mx=mx, my=my, v=v, vx=vx, vy=vy, nframes=nframes)
    else:
        v_mx /= nframes
        v_ave /= nframes
        np.savez_compressed(fout, rho=rho, mx=mx, my=my, v=v, vx=vx, vy=vy, v_mx=v_mx, v_ave=v_ave, nframes=nframes)


def update_time_ave_profiles():
    root = f"{root_tahmineh}/QS_MPI_noself"
    fname = f"{root}/L320_80_r10_10_10_h0.333/L320_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 40
    time_ave(fname, beg_frame)

    fname = f"{root}/L960_80_r10_10_10_h0.333/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 10
    time_ave(fname, beg_frame)

    fname = f"{root}/L480_80_r10_10_10_h0.333/L480_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 264
    time_ave(fname, beg_frame)

    fname = f"{root}/L960_160_r10_10_10_varied_h/L960_160_Dr0.100_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 60, None, center_type="B")

    fname = f"{root}/L960_160_r10_10_10_varied_h/L960_160_Dr0.100_k0.70_p6.25_7.15_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 30, None, center_type="B")

    fname = f"{root}/L960_160_r10_10_10_varied_h/L960_160_Dr0.100_k0.70_p6.25_8.15_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 60, None, center_type="B")

    fname = f"{root}/L960_160_r10_10_10_h0.333/L960_160_Dr0.100_k0.70_p8.625_13.36_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.gsd"
    time_ave(fname, 50, None, center_type="both", inc_v_mx=True, folder="data/time_ave_profiles/CCB_gas")

    fname = f"{root}/L960_160_r10_10_10_h0.333/L960_160_Dr0.100_k0.70_p7.6525_10.5675_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.gsd"
    time_ave(fname, 50, None, center_type="both", inc_v_mx=True, folder="data/time_ave_profiles/CCB_gas")

    fname = f"{root}/L960_160_r10_10_10_h0.333/L960_160_Dr0.100_k0.70_p9.597_16.152_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.gsd"
    time_ave(fname, 50, None, center_type="both", inc_v_mx=True, folder="data/time_ave_profiles/CCB_gas")

    fname = f"{root}/L480_160_r10_10_10_h0.333/L480_160_Dr0.100_k0.70_p19_8.8_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    time_ave(fname, 20, None, center_type="B")

    fname = f"{root}/L960_160_r10_10_10_h0.333/L960_160_Dr0.100_k0.70_p19_8.8_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    time_ave(fname, 20, None, center_type="B")

    fname = f"{root}/LB_G_L240_40_r20/L240_40_Dr0.400_k0.70_p7_28_r20_20_20_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 3, None, center_type="B")

    fname = f"{root}/LB_G_L120_40_r20/L120_40_Dr0.500_k0.70_p7.1_28_r20_20_20_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 10, None, center_type="B")

    fname = f"{root}/LB_G_L120_40_r20/L120_40_Dr1.000_k0.70_p8.3_28_r20_20_20_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 30, None, center_type="B")

    fname = f"{root}/LB_G_L120_40_r20/L120_40_Dr1.000_k0.70_p8_28_r20_20_20_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 20, None, center_type="B")

    fname = f"{root}/LB_G_L120_40_r20/L120_40_Dr1.000_k0.70_p8_28_r20_20_20_e-2.000_J2.000_-2.000_h0.100_2001.gsd"
    time_ave(fname, 20, None, center_type="B")

    fname = f"{root}/LB_G_L120_40_r10/L120_40_Dr0.100_k0.70_p3.6_12.5_r10_10_10_e-2.000_J0.500_-0.500_h0.100_2000.gsd"
    time_ave(fname, 20, None, center_type="B")

    fname = f"{root}/L1920_160_r10_10_10_varied_h/L1920_160_Dr0.100_k0.70_p6.25_8.15_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 10, None, center_type="B")

    fname = f"{root}/L1920_160_r10_10_10_varied_h/L1920_160_Dr0.100_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 10, None, center_type="B")

    fname = f"{root}/L3840_160_r10_10_10_varied_h/L3840_160_Dr0.100_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    time_ave(fname, 10, None, center_type="B")


def show_time_ave_profiles():
    # fnames = ["data/time_ave_profiles/L320_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
    #           "data/time_ave_profiles/L480_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
    #           "data/time_ave_profiles/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    # ]

    # fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    # for fname in fnames:
    #     with np.load(fname, "r") as data:
    #         rho = data["rho"]
    #         Lx, Ly = get_L(fname)
    #         dx = 1
    #         x = np.arange(Lx) * dx + dx / 2
    #     ax1.plot(x / Lx, rho[0], label=r"$L_x=%d$" % Lx)
    #     ax2.plot(x / Lx, rho[1])

    # ax2.set_xlabel(r"$x/L_x$")
    # ax1.set_ylabel(r"$\rho_A$")
    # ax2.set_ylabel(r"$\rho_B$")

    # ax1.set_xlim(0, 1)
    # ax1.legend()
    # plt.show()
    # plt.close()

    # fname = "data/time_ave_profiles/L960_160_Dr0.100_k0.70_p8.625_13.36_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.npz"
    # fname = "data/time_ave_profiles/L960_160_Dr0.100_k0.70_p7.6525_10.5675_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.npz"
    fname = "data/time_ave_profiles/CCB_gas/L960_160_Dr0.100_k0.70_p9.597_16.152_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.npz"

    with np.load(fname, "r") as data:
        print(data.files)
        rho = data["rho"]
        v_rho = data["v"]
        mx = data["mx"]
        vx = data["vx"]
        v_mx = data["v_mx"]
        v_ave = data["v_ave"]
        Lx, Ly = get_L(fname)
        dx = 1
        x = np.arange(Lx) * dx + dx / 2

    d_v_rho = np.zeros((2, x.size))
    print(d_v_rho.shape)
    for i in range(x.size):
        d_v_rho[:, i] = -(v_rho[:, (i+1)%x.size] - v_rho[:, i-1]) / (4 * 0.1)
    

    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    ax1.plot(x / Lx, rho[0], label=r"$\langle \rho_S\rangle_{y,t}$")
    ax2.plot(x / Lx, rho[1], label=r"$\langle \rho_S\rangle_{y,t}$")

    ax1.plot(x / Lx, v_rho[0], label=r"$\langle \rho_S v_S\rangle_{y,t}$")
    ax2.plot(x / Lx, v_rho[1], label=r"$\langle \rho_S v_S\rangle_{y,t}$")

    ax1.plot(x / Lx, v_ave[0] * rho[0], label=r"$\langle \rho_S\rangle_{y,t} \langle v_S\rangle_{y,t}$")
    ax2.plot(x / Lx, v_ave[1] * rho[1], label=r"$\langle \rho_S\rangle_{y,t} \langle v_S\rangle_{y,t}$")

    # ax1.plot(x / Lx, v_rho[0], label=r"$S=A$")
    # ax2.plot(x / Lx, v_rho[1], label=r"$S=B$")
    # ax1.plot(x / Lx, mx[0], label=r"$-\partial_x \langle v_S \phi_S\rangle_{y, t}/(2D_r)$")
    # ax2.plot(x / Lx, mx[1], label=r"$-\langle\frac{\partial_x (v_S \phi_S)}{2D_r}\rangle_{y, t}$")

    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$S=A$")
    ax2.set_ylabel(r"$S=B$")

    ax1.set_xlim(0, 1)
    ax1.legend()
    plt.show()
    plt.close()

    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    ax1.plot(x / Lx, d_v_rho[0], label=r"$\langle\mathbf{p}_s \cdot \hat{\mathbf{e}}_x\rangle_{y,t}$")
    ax2.plot(x / Lx, d_v_rho[1], label=r"$\langle\mathbf{p}_s \cdot \hat{\mathbf{e}}_x\rangle_{y,t}$")
    ax1.plot(x / Lx, mx[0], label=r"$-\partial_x \langle v_S \phi_S\rangle_{y, t}/(2D_r)$")
    ax2.plot(x / Lx, mx[1], label=r"$-\langle\frac{\partial_x (v_S \phi_S)}{2D_r}\rangle_{y, t}$")

    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$S=A$")
    ax2.set_ylabel(r"$S=B$")

    ax1.set_xlim(0, 1)
    ax1.legend(fontsize="large")
    plt.show()
    plt.close()


    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    ax1.plot(x / Lx, v_mx[0], label=r"$\langle v_S \mathbf{p}_S\rangle_{y,t}\cdot \hat{\mathbf{e}}_x$")
    ax2.plot(x / Lx, v_mx[1], label=r"$\langle v_S \mathbf{p}_S\rangle_{y,t}\cdot \hat{\mathbf{e}}_x$")
    ax1.plot(x / Lx, vx[0], label=r"$\langle \rho_{S}\mathbf{v}_S\rangle_{y, t}\cdot \hat{\mathbf{e}}_x$")
    ax2.plot(x / Lx, vx[1], label=r"$\langle \rho_{S}\mathbf{v}_S\rangle_{y, t}\cdot \hat{\mathbf{e}}_x$")

    ax1.plot(x / Lx, v_ave[0] * mx[0], label=r"$\langle v_S\rangle_{y,t} \langle \mathbf{p}_S\rangle_{y,t}\cdot \hat{\mathbf{e}}_x$")
    ax2.plot(x / Lx, v_ave[1] * mx[1], label=r"$\langle v_S\rangle_{y, t} \langle \mathbf{p}_S\rangle_{y,t}\cdot \hat{\mathbf{e}}_x$")

    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$S=A$")
    ax2.set_ylabel(r"$S=B$")

    ax1.set_xlim(0, 1)
    ax1.legend(fontsize="large")
    plt.show()
    plt.close()



def estimate_binodals():
    fname = "data/time_ave_profiles/CCB_gas/L960_160_Dr0.100_k0.70_p9.597_16.152_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2100.npz"
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    with np.load(fname, "r") as data:
        rho = data["rho"]
        Lx, Ly = get_L(fname)
        dx = 1
        x = np.arange(Lx) * dx + dx / 2
    ax1.plot(x / Lx, rho[0], label=r"$L_x=%d$" % Lx)
    ax2.plot(x / Lx, rho[1])

    xg = [0.15, 0.85]
    # xg = [0.1, 0.9]
    xl = [0.3, 0.7]
    # xl = [0.3, 0.7]
    mask_g = np.logical_or(x/Lx < xg[0], x/Lx > xg[1])
    mask_l = np.logical_and(x/Lx > xl[0], x/Lx < xl[1])

    rA_g = np.mean(rho[0][mask_g])
    rA_l = np.mean(rho[0][mask_l])
    rB_g = np.mean(rho[1][mask_g])
    rB_l = np.mean(rho[1][mask_l])

    ax1.axhline(rA_g, 0, xg[0], linestyle="dashed", c="tab:green")
    ax1.axhline(rA_g, xg[1], 1, linestyle="dashed", c="tab:green")
    ax1.axhline(rA_l, xl[0], xl[1], linestyle="dashed", c="tab:orange")

    ax2.axhline(rB_g, 0, xg[0], linestyle="dashed", c="tab:green")
    ax2.axhline(rB_g, xg[1], 1, linestyle="dashed", c="tab:green")
    ax2.axhline(rB_l, xl[0], xl[1], linestyle="dashed", c="tab:orange")

    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")

    rA_0 = 8.625
    rB_0 = 13.36
    print(rA_g, rB_g, rA_l, rB_l)
    print((rA_0-rA_g)/(rA_l-rA_g), (rB_0-rB_g)/(rB_l-rB_g))
    ax1.set_xlim(0, 1)
    ax1.legend()
    plt.show()
    plt.close()

    # phiA = 7
    # phiB = 4
    # plt.plot(phiA, phiB, "o")
    # plt.plot(rA_g, rB_g, "<")
    # plt.plot(rA_l, rB_l, ">")
    # plt.plot([rA_g, rA_l], [rB_g, rB_l], "--", c="tab:grey")


    # rA_new = 5.6

    # rB_new = (rA_new - phiA) * (rB_l - rB_g) / (rA_l - rA_g) + phiB

    # plt.plot(rA_new, rB_new, "s", fillstyle="none")
    # print(rA_new, rB_new)

    # plt.show()
    # plt.close()


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


def show_instant_profiles():
    fname = f"{root_tahmineh}/QS_MPI_noself/L1920_640_r10_10_10_h0.333/L1920_640_Dr1.000_k0.70_p6.48_7.2_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L960_160_r5_5_5_varied_Dt_h/L960_160_Dr0.100_Dt0_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_h0.100_2000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L1920_160_r5_5_5_varied_Dt_h/L1920_160_Dr0.100_Dt0_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_h0.100_2000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L960_160_r5_5_5_varied_Dt_h/L960_160_Dr0.100_Dt0_k0.70_p2.25_7.35_r5_5_5_e-2.000_J0.500_-0.500_h0.100_2000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L3840_160_r10_10_10_varied_h/L3840_160_Dr0.100_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    # fname = f"{root_sohrab}/QS_not_inc_self/L10_4_v1_1/L10_4_Dr0.100_r200_e0.000_0.000_J1.800000_-1.800000_h0.100_1001.gsd"

    # fname = f"{root_tahmineh}/QS_MPI_noself/L960_80_r10_10_10_h0.333/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2120.gsd"

    Lx, Ly = get_L(fname)
    rho0 = 10
    dx = 2
    if isinstance(dx, int):
        nx = Lx // dx
    else:
        nx = int(round(Lx/dx))
    x = np.arange(nx) * dx + dx / 2
    print(Lx, Ly, x.size)

    beg_frame = 7
    end_frame = None
    with fl.open(name=fname, mode="rb") as fgsd:
        if end_frame is None:
            end_frame = fgsd.nframes
        for i_frame in range(beg_frame, end_frame):
            print("frame: %d/%d" % (i_frame, end_frame))
            fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
            rho, mx, my, v, vx, vy = get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None)
            ax1.plot(x, rho[0]/rho0, label="A")
            ax1.plot(x, rho[1]/rho0, label="B")

            # rho_thresh = 12.5
            # ax2.axhline(rho_thresh, linestyle="dashed")
            # xs1 = locate_interfaces(x, rho[1], rho_thresh, Lx, True)
            # xs2 = locate_interfaces(x, rho[1], rho_thresh, Lx, False)

            # for xi in xs1:
            #     ax2.axvline(xi, linestyle="dotted", c="tab:green")
            # for xi in xs2:
            #     ax2.axvline(xi, linestyle="dotted", c="tab:red")

            ax1.set_ylim(0)
            ax2.plot(x, mx[0]/rho0)
            ax2.plot(x, mx[1]/rho0)
            ax2.set_xlabel(r"$x/L_x$")
            ax1.set_ylabel(r"$\rho/\rho_0$")
            ax2.set_ylabel(r"$m_x$")
            ax1.legend()
            ax1.set_xlim(0, Lx)
            plt.show()
            plt.close()
        

def show_profiles_varied_time():
    import matplotlib as mpl

    fname = f"{root_sohrab}/QS_not_inc_self/L10_4_v1_1/L10_4_Dr0.100_r200_e0.000_0.000_J1.800000_-1.800000_h0.100_1001.gsd"
    Lx, Ly = get_L(fname)
    rho0 = 200
    dx = 0.2
    if isinstance(dx, int):
        nx = Lx // dx
    else:
        nx = int(round(Lx/dx))
    x = np.arange(nx) * dx + dx / 2
    print(Lx, Ly, x.size)
    h = 0.1

    beg_frame =82
    end_frame =90
    sep = 1
    n_lines = (end_frame - beg_frame)//sep
    cmap = mpl.colormaps['plasma']
    colors = cmap(np.linspace(0, 1, n_lines))
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)

    t_arr = np.zeros(n_lines)
    rho_gas = np.zeros((n_lines, 2))
    with fl.open(name=fname, mode="rb") as fgsd:
        if end_frame is None:
            end_frame = fgsd.nframes
        for j, i_frame in enumerate(range(beg_frame, end_frame, sep)):
            print("frame: %d/%d" % (i_frame, end_frame))
            rho, mx, my, v, vx, vy = get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None)

            t = i_frame * h
            ax1.plot(x, rho[0]/rho0, label=r"$%g$" % t, c=colors[j], linestyle=":")
            ax2.plot(x, rho[1]/rho0, label=r"$%g$" % t, c=colors[j], linestyle=':')

            mask = np.logical_and(x >= 6, x <=7)
            rho_gas[j][0] = np.mean(rho[0][mask])
            rho_gas[j][1] = np.mean(rho[1][mask])
            t_arr[j] = t

    ax2.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$\rho_A/\rho_0$")
    ax2.set_ylabel(r"$\rho_B/\rho_0$")
    ax1.legend(title=r"$t=$")
    ax1.set_xlim(0, Lx)
    plt.show()
    plt.close()

    plt.plot(t_arr, rho_gas[:, 0], "-o", label=r"$S=A$")
    plt.plot(t_arr, rho_gas[:, 1], "-s", label=r"$S=B$")
    plt.xlabel(r"$t$", fontsize="x-large")
    plt.ylabel(r"$\langle \rho_S\rangle_{x\in [6, 7]}$", fontsize='x-large')
    plt.legend()
    plt.show()
    plt.close()



def update_space_time_profiles():
    from coarse_grain import coarse_grain_all_1D

    Lx = 960
    Ly = 40
    dx = 2
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10"
    coarse_grain_all_1D(prefix, dx, False)

    Lx = 240
    Ly = 20
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 240
    Ly = 40
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 480
    Ly = 40
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 240
    Ly = 80
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 240
    Ly = 160
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 120
    Ly = 40
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h"
    coarse_grain_all_1D(prefix, 1, False)

    Lx = 120
    Ly = 40
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r20_20_20_varied_h"
    coarse_grain_all_1D(prefix, 1, False, pat="*_p20_*.gsd")

    Lx = 1920
    Ly = 160
    dx = 2
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r5_5_5_varied_Dt_h"
    coarse_grain_all_1D(prefix, dx, False)

    Lx = 960
    Ly = 160
    dx = 2
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r5_5_5_varied_Dt_h"
    coarse_grain_all_1D(prefix, dx, False)

def space_time():
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_80_r10_10_10_h0.333"
    Lx = 1280
    seed = 2000
    fnpz =  f"{prefix}/profile_1D_dx1/L{Lx:g}_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"
    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
    
    fig = plt.figure(layout='constrained', figsize=(8, 6))
    subfig = fig.subfigures(1, 1)
    axes = subfig.subplots(2, 1, sharex=True)
    ax1, ax2 = axes
    with np.load(fnpz, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]

        rhoA = fields[:, 0]
        rhoB = fields[:, 1]
        im1 = ax1.imshow(rhoB[60:], origin="lower", aspect="auto", extent=[0, Lx, t[60]/3, t[-1]/3], vmax=30)
        im2 = ax2.imshow(rhoB[:61], origin="lower", aspect="auto", extent=[0, Lx, np.log10(t[0]/3), np.log10(t[60]/3)], vmax=30)

        ax1.axvline(875, linestyle="dashed", c="w")
        ax1.axvline(925, linestyle="dashed", c="w")

        ax1.axvline(450, linestyle="dotted", c="w")
        ax1.axvline(500, linestyle="dotted", c="w")

        ax1.axvline(1150, linestyle="dotted", c="w")
        ax1.axvline(1200, linestyle="dotted", c="w")
        ax1.axvline(200, linestyle="dotted", c="w")
        ax1.axvline(250, linestyle="dotted", c="w")

        ax1.set_ylabel(r"$t$", fontsize="x-large")
        ax2.set_ylabel(r"$\log t$", fontsize="x-large")
        ax2.set_xlabel(r"$x$", fontsize="x-large")

        subfig.colorbar(im1, shrink=0.6, ax=axes, extend="max")
        subfig.text(0.91, 0.82, r"$\langle \rho_B\rangle_{y}$", fontsize="x-large")
        # subfig.suptitle('Right plots', fontsize='x-large')
        plt.show()
        # plt.savefig("coarsening_1D_L1280.pdf")
        plt.close()


        fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
        ax1.plot(t[60:]/3, np.mean(rhoA[60:, 875:925]/10, axis=1), label="S=A")
        ax1.plot(t[60:]/3, np.mean(rhoB[60:, 875:925]/10, axis=1), label="S=B")

        ax1.plot(t[60:]/3, np.mean(rhoA[60:, 200:250]/10, axis=1), label="S=A")
        ax1.plot(t[60:]/3, np.mean(rhoB[60:, 200:250]/10, axis=1), label="S=B")

        ax2.plot(t[60:]/3, np.mean(rhoA[60:, 450:500]/10, axis=1), label="S=A")
        ax2.plot(t[60:]/3, np.mean(rhoB[60:, 450:500]/10, axis=1), label="S=B")

        ax2.plot(t[60:]/3, np.mean(rhoA[60:, 1150:1200]/10, axis=1), label="S=A")
        ax2.plot(t[60:]/3, np.mean(rhoB[60:, 1150:1200]/10, axis=1), label="S=B")
        ax1.set_xlim(1e6)
        ax1.legend()
        ax2.legend()
        ax2.set_xlabel(r"$t$", fontsize="x-large")
        ax1.set_ylabel(r"$\rho_S^{\rm gas}/\rho_0$", fontsize="x-large")
        ax2.set_ylabel(r"$\rho_S^{\rm liquid}/\rho_0$", fontsize="x-large")

        plt.show()
        plt.close()


def plot_space_time_varied_Ly(rA=9):
    Lx = 240
    if rA == 9:
        ncut_dict = {20: 0, 40: 600, 80: 1, 160: 1}
        seed_dict = {20: 1000, 40: 1000, 80: 2000, 160: 2000}
    elif rA == 8:
        ncut_dict = {40: 1, 80: 1, 160: 1}
        seed_dict = {40: 1000, 80: 2000, 160: 2000}
    titles = [r"(a) $L_y=%d$", r"(b) $L_y=%d$", r"(c) $L_y=%d$", r"(d) $L_y=%d$"]

    fig, axes = plt.subplots(1, len(ncut_dict), constrained_layout=True, sharey=False)
    for i, Ly in enumerate(ncut_dict.keys()):
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h/profile_1D_dx1"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p{rA}_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_{seed_dict[Ly]}.npz"


        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]

            ncut = ncut_dict[Ly]
            if ncut > 0:
                T = (t0[-1] - t0[ncut-1]) * 0.1
            else:
                T = (t0[-1]) * 0.1

            rhoA = fields0[ncut:, 0]
            rhoB = fields0[ncut:, 1]
            axes[i].imshow(rhoB, origin="lower", extent=[0, Lx, 0, T], aspect="auto")
            axes[i].set_title(titles[i] % Ly)
            axes[i].set_xlabel(r"$x$")

    axes[0].set_ylabel(r"$t$")
    plt.show()
    plt.close()


def plot_space_time_varied_Lx():
    Ly = 40
    ncut_dict = {120: 400, 240: 100, 480: 1}
    seed_dict = {120: 1000, 240: 1000, 480: 2000}
    titles = [r"(a) $L_x=%d$", r"(b) $L_x=%d$", r"(c) $L_x=%d$"]

    fig, axes = plt.subplots(1, len(ncut_dict), constrained_layout=True, sharey=False)
    for i, Lx in enumerate(ncut_dict.keys()):
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h/profile_1D_dx1"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p9_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_{seed_dict[Lx]}.npz"


        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]

            ncut = ncut_dict[Lx]
            if ncut > 0:
                T = (t0[-1] - t0[ncut-1]) * 0.1
            else:
                T = t0[-1] * 0.1

            rhoA = fields0[ncut:, 0]
            rhoB = fields0[ncut:, 1]
            axes[i].imshow(rhoB, origin="lower", extent=[0, Lx, 0, T], aspect="auto")
            axes[i].set_title(titles[i] % Lx)
            axes[i].set_xlabel(r"$x$")
            axes[i].set_ylabel(r"$t$")

    plt.show()
    plt.close()


def plot_space_time_varied_rhoA():
    Lx = 240
    Ly = 40
    ncut_dict = {8: 0, 9: 100, 10: 0, 12.5: 100}
    seed_dict = {8: 1000, 9: 1000, 10: 1000, 12.5: 1000}
    titles = [r"(a) $\rho_A=%g$", r"(b) $\rho_A=%g$", r"(c) $\rho_A=%g$", r"(d) $\rho_A=%g$"]

    fig, axes = plt.subplots(1, len(ncut_dict), constrained_layout=True, sharey=False)
    for i, rA in enumerate(ncut_dict.keys()):
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h/profile_1D_dx1"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p{rA:g}_34_r10_10_10_e-2.000_J0.500_-0.500_h0.100_{seed_dict[rA]}.npz"


        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]

            ncut = ncut_dict[rA]
            if ncut > 0:
                T = (t0[-1] - t0[ncut-1]) * 0.1
            else:
                T = t0[-1] * 0.1

            rhoA = fields0[ncut:, 0]
            rhoB = fields0[ncut:, 1]
            axes[i].imshow(rhoB, origin="lower", extent=[0, Lx, 0, T], aspect="auto")
            axes[i].set_title(titles[i] % rA)
            axes[i].set_xlabel(r"$x$")
            axes[i].set_ylabel(r"$t$")

    plt.show()
    plt.close()


def plot_space_time_varied_rhoB():
    Lx = 120
    Ly = 40
    rA = 9
    ncut_dict = {30: 100, 32: 100, 34: 400, 36: 0}
    seed_dict = {30: 1000, 32: 1000, 34: 1000, 36: 1000}
    titles = [r"(a) $\rho_B=%g$", r"(b) $\rho_B=%g$", r"(c) $\rho_B=%g$", r"(d) $\rho_B=%g$"]

    fig, axes = plt.subplots(1, len(ncut_dict), constrained_layout=True, sharey=False)
    for i, rB in enumerate(ncut_dict.keys()):
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_varied_h/profile_1D_dx1"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p{rA:g}_{rB:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.100_{seed_dict[rB]}.npz"


        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]

            ncut = ncut_dict[rB]
            dt = (t0[-1] - t0[-2]) * 0.1
            T = (t0[-1] - t0[ncut]) * 0.1 + dt
 
            rhoA = fields0[ncut:, 0]
            rhoB = fields0[ncut:, 1]
            axes[i].imshow(rhoB, origin="lower", extent=[0, Lx, 0, T], aspect="auto")
            axes[i].set_title(titles[i] % rB)
            axes[i].set_xlabel(r"$x$")
            axes[i].set_ylabel(r"$t$")

    plt.show()
    plt.close()

def plot_space_time_varied_rhoB_10():
    Lx = 120
    Ly = 40
    rA = 20
    ncut_dict = {50: 1, 68: 9, 76:1}
    seed_dict = {50: 1000, 68: 2000, 76:2000}
    titles = [r"(a) $\rho_B=%g$", r"(b) $\rho_B=%g$", r"(c) $\rho_B=%g$", r"(d) $\rho_B=%g$"]

    fig, axes = plt.subplots(1, len(ncut_dict), constrained_layout=True, sharey=False)
    for i, rB in enumerate(ncut_dict.keys()):
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r20_20_20_varied_h/profile_1D_dx1"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p{rA:g}_{rB:g}_r20_20_20_e-2.000_J0.500_-0.500_h0.100_{seed_dict[rB]}.npz"


        with np.load(fname, "r") as data:
            t0, x0, fields0 = data["t"], data["x"], data["fields"]

            ncut = ncut_dict[rB]
            dt = (t0[-1] - t0[-2]) * 0.1
            T = (t0[-1] - t0[ncut]) * 0.1 + dt
 
            rhoA = fields0[ncut:, 0]
            rhoB = fields0[ncut:, 1]
            axes[i].imshow(rhoB, origin="lower", extent=[0, Lx, 0, T], aspect="auto")
            axes[i].set_title(titles[i] % rB)
            axes[i].set_xlabel(r"$x$")
            axes[i].set_ylabel(r"$t$")

    plt.show()
    plt.close()

def plot_space_time_fields():
    Lx = 1920
    Ly = 160
    rA = 2.25
    rB = 7.35
    r0 = 5
    seed_dict = {160: 2000}
    ncut_dict = {160: 0}
    dx = 2

    if (Lx == 1920 or Lx == 960) and Ly == 160:
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r{r0:g}_{r0:g}_{r0:g}_varied_Dt_h/profile_1D_dx{dx:g}"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_Dt0_k0.70_p{rA:g}_{rB:g}_r{r0:g}_{r0:g}_{r0:g}_e-2.000_J0.500_-0.500_h0.100_{seed_dict[Ly]}.npz"
    else:
        if Lx == 960 and Ly == 40:
            folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r{r0:g}_{r0:g}_{r0:g}/profile_1D_dx{dx:g}"
        else:
            folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r{r0:g}_{r0:g}_{r0:g}_varied_h/profile_1D_dx{dx:g}"
        fname = f"{folder}/L{Lx}_{Ly}_Dr0.100_k0.70_p{rA:g}_{rB:g}_r{r0:g}_{r0:g}_{r0:g}_e-2.000_J0.500_-0.500_h0.100_{seed_dict[Ly]}.npz"

    fig, ax = plt.subplots(1, 1, constrained_layout=True, sharey=False)

    with np.load(fname, "r") as data:
        t0, x0, fields0 = data["t"], data["x"], data["fields"]

        ncut = ncut_dict[Ly]
        dt = (t0[ncut+1] - t0[ncut]) * 0.1
        T = (t0[-1] - t0[ncut]) * 0.1 + dt

        rhoA = fields0[ncut:, 0]
        rhoB = fields0[ncut:, 1]
        mask = t0 % 40000 == 0
        ax.imshow(rhoB[mask], origin="lower", extent=[0, Lx, 0, T], aspect="auto")
        # ax.set_title(titles[i] % Ly)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$t$")
    
    plt.show()
    plt.close()


def get_interface_locs(Lx, Ly=80, h=1/3.):
    if Lx == 1280:
        seed_arr = [1000, 2129042, 2030, 2240]
    elif Lx == 1120:
        seed_arr = [2000, 2003]
    elif Lx== 960:
        seed_arr = [2030, 2060, 2090, 2120]
    
    if Lx == 1120:
        dx = 8
    else:
        dx = 6
    for seed in seed_arr:
        fname = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_{Ly}_r10_10_10_h0.333/L{Lx}_{Ly}_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed}.gsd"
        x = np.arange(Lx//dx) * dx + dx / 2
        fout = f"data/interface_locs/{os.path.basename(fname).replace('.gsd', '.npz')}"

        if Lx == 960:
            beg_frame = 0
        elif Lx == 1120:
            beg_frame = 0
        elif Lx == 1280 and seed == 1000:
            beg_frame = 170
        else:
            beg_frame = 0
        end_frame = None

        with fl.open(name=fname, mode="rb") as fgsd:
            if end_frame is None:
                end_frame = fgsd.nframes
            print("find", end_frame, "in", os.path.basename(fname))
            xs = np.zeros((end_frame-beg_frame, 4))
            t_arr = np.zeros(end_frame-beg_frame)
            if os.path.exists(fout):
                with np.load(fout, "r") as data:
                    existed_frames = data["t"].size
                    t_arr[:existed_frames] = data["t"]
                    xs[:existed_frames] = data["xs"]
            else:
                existed_frames = 0
            for i_frame in range(beg_frame + existed_frames, end_frame):
                print("frame: %d/%d" % (i_frame, end_frame))
                rho_thresh = 12.5
                rho, mx, my, v, vx, vy = get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None)

                if Lx == 960 and seed == 2120:
                    rho = np.roll(rho, -8, axis=1)
                elif Lx == 1280 and seed == 1000:
                    rho = np.roll(rho, -30, axis=1)
                xs1 = locate_interfaces(x, rho[1], rho_thresh, Lx, True)
                xs2 = locate_interfaces(x, rho[1], rho_thresh, Lx, False)
                if i_frame == 0:
                    t_arr[0] = 0
                else:
                    t_arr[i_frame-beg_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0] * h
                if len(xs1) != 2 or len(xs2) != 2:
                    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
                    ax1.plot(x, rho[0])
                    ax2.plot(x, rho[1])
                    ax2.axhline(12.5, linestyle="dashed")

                    for xi in xs1:
                        ax2.axvline(xi, linestyle="dotted", c="tab:green")
                    for xi in xs2:
                        ax2.axvline(xi, linestyle="dotted", c="tab:red")

                    ax2.set_xlabel(r"$x/L_x$")
                    ax1.set_ylabel(r"$\rho_A$")
                    ax2.set_ylabel(r"$\rho_B$")
                    plt.show()
                    plt.close()
                    xs[i_frame-beg_frame] = xs[i_frame-beg_frame-1]
                else:
                    xs[i_frame-beg_frame][0] = xs1[0]
                    xs[i_frame-beg_frame][1] = xs2[0]
                    xs[i_frame-beg_frame][2] = xs1[1]
                    xs[i_frame-beg_frame][3] = xs2[1]
        np.savez_compressed(fout, t=t_arr, xs=xs)


def plot_spacing_vs_t(Lx=960):
    get_interface_locs(Lx)
    if Lx == 960:
        seed_arr = [2030, 2060, 2090, 2120]
    elif Lx == 1280:
        seed_arr = [1000, 2129042, 2030, 2240]
    elif Lx == 1120:
        seed_arr = [2000, 2003]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    for seed in seed_arr:
        f_interface = f"data/interface_locs/L{Lx}_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed}.npz"
        with np.load(f_interface, "r") as data:
            t, xs = data["t"], data["xs"]
            d1, d2 = np.zeros((2, t.size))
            w1, w2 = np.zeros((2, t.size))

            for i in range(t.size):
                d1[i] = xs[i][0] + Lx - xs[i][3]
                d2[i] = xs[i][2] - xs[i][1]

                w1[i] = xs[i][1] - xs[i][0]
                w2[i] = xs[i][3] - xs[i][2]
        
    
        line, = ax1.plot(t-t[0], d1, "o")
        ax1.plot(t-t[0], d2, "s", c=line.get_c())
        # line, = ax2.plot(t-t[0], np.abs(d2-d1), "o")
        line, = ax2.plot(t-t[0], d2-d1, "o")

        # ax2.plot(t, w2, "s", c=line.get_c())
    plt.show()
    plt.close()


if __name__ == "__main__":
    # update_time_ave_profiles()
    # show_time_ave_profiles()
    # show_instant_profiles()
    # show_profiles_varied_time()

    plot_spacing_vs_t(Lx=1280)

    # update_space_time_profiles()
    # plot_space_time_varied_Ly(rA=8)
    # plot_space_time_varied_rhoB_10()
    # plot_space_time_varied_Lx()
    # plot_space_time_fields()
    
    # estimate_binodals()
