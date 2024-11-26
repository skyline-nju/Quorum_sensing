import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_one_frame
from gsd import fl
from center_traj import find_center

root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_cg_fields(i_frame, fgsd, dx, Lx, Ly, center=None):
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)
    bin_area = dx ** 2

    par_num = np.zeros((2, nrows, ncols), int)
    mx, my, v = np.zeros((3, 2, nrows, ncols))
    vx, vy = np.zeros((2, 2, nrows, ncols))

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
        y = pos[:,1] + Ly / 2
    else:
        x = pos[:, 0] - center[0] + Lx/2
        y = pos[:, 1] - center[1] + Ly/2
        x[x<0] += Lx
        x[x>=Lx] -= Lx
        y[y<0] += Ly
        y[y>Ly] -= Ly

    ux = np.cos(pos[: ,2])
    uy = np.sin(pos[:, 2])

    col = (x/dx).astype(int)
    row = (y/dx).astype(int)
    col[col < 0] += ncols
    col[col >= ncols] -= ncols
    row[row < 0] += nrows
    row[row >= nrows] -= nrows

    for j in range(x.size):
        my_id = id_type[j]
        my_row = row[j]
        my_col = col[j]
        par_num[my_id, my_row, my_col] += 1
        mx[my_id, my_row, my_col] += ux[j]
        my[my_id, my_row, my_col] += uy[j]
        v[my_id, my_row, my_col] += speed[j]
        vx[my_id, my_row, my_col] += ux[j] * speed[j]
        vy[my_id, my_row, my_col] += uy[j] * speed[j]
    rho = par_num / bin_area
    mx /= bin_area
    my /= bin_area
    v /= bin_area
    vx /= bin_area
    vy /= bin_area
    return rho, mx, my, v, vx, vy


def get_L(fname):
    s = os.path.basename(fname).lstrip("L").split("_")
    Lx = int(s[0])
    Ly = int(s[1])
    return Lx, Ly


def time_ave(fname, beg_frame=20, center_type="both", dx=1):
    basename = os.path.basename(fname)
    print("handling", basename)
    fout = f"data/time_ave_fields/{basename.replace('.gsd', '.npz')}"
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
    else:
        nrows = int(round(Ly / dx))
        ncols = int(round(Lx / dx))
        rho, mx, my, v, vx, vy = np.zeros((6, 2, nrows, ncols))
        nframes = 0
    with fl.open(name=fname, mode="rb") as fgsd:
        count = 0
        for i in range(beg_frame + nframes, fgsd.nframes):
            if center_type is None:
                center = None
            else:
                x_c, y_c = find_center(i, fgsd, Lx, Ly, which=center_type)
                center = [x_c, y_c]
            rho_cur, mx_cur, my_cur, v_cur, vx_cur, vy_cur = get_cg_fields(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
            rho += rho_cur
            mx += mx_cur
            my += my_cur
            v += v_cur
            vx += vx_cur
            vy += vy_cur
            count +=1
            print("frame:%d/%d" % (i, fgsd.nframes))
    nframes += count
    rho /= nframes
    mx /= nframes
    my /= nframes
    v /= nframes
    vx /= nframes
    vy /= nframes

    np.savez_compressed(fout, rho=rho, mx=mx, my=my, v=v, vx=vx, vy=vy, nframes=nframes)


def CCB_profiles_varied_L():
    f1 = "data/time_ave_fields/L160_160_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    f2 = "data/time_ave_fields/L320_320_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    f3 = "data/time_ave_fields/L640_640_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    f4 = "data/time_ave_fields/L1280_1280_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"


    with np.load(f1, "r") as data:
        rho1, mx1, my1, v1, vx1, vy1 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f2, "r") as data:
        rho2, mx2, my2, v2, vx2, vy2 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f3, "r") as data:
        rho3, mx3, my3, v3, vx3, vy3 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]
    
    with np.load(f4, "r") as data:
        rho4, mx4, my4, v4, vx4, vy4 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.3), constrained_layout=True, sharex=True, sharey=True)
    ax1.imshow(rho3[0], origin="lower")
    ax2.imshow(rho3[1], origin="lower")
    ax1.set_title(r"$\rho_A$")
    ax2.set_title(r"$\rho_B$")

    ax1.set_xlabel(r"$x$")
    ax2.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")

    plt.show()
    plt.close()

    dx = 1
    x1 = np.arange(160) * dx + dx / 2
    x2 = np.arange(320) * dx + dx / 2
    x3 = np.arange(640) * dx + dx / 2
    x4 = np.arange(1280) * dx + dx / 2

    x1 /= 160
    x2 /= 320
    x3 /= 640
    # x4 /= 1280

    # mx1 /= rho1
    # mx2 /= rho2
    # mx3 /= rho3
    # rho1_x = np.mean(mx1[:, 75:85], axis=1)
    # rho2_x = np.mean(mx2[:, 150:170], axis=1)
    # rho3_x = np.mean(mx3[:, 300:340], axis=1)

    rho1_x = np.mean(rho1[:, 75:85], axis=1)
    rho2_x = np.mean(rho2[:, 150:170], axis=1)
    rho3_x = np.mean(rho3[:, 300:340], axis=1)
    rho4_x = np.mean(rho4[:, 600:680], axis=1)

    # rho1_x = np.mean(rho1[:, :], axis=1)
    # rho2_x = np.mean(rho2[:, :], axis=1)
    # rho3_x = np.mean(rho3[:, :], axis=1)

    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    line1, = ax1.plot(x1, rho1_x[0], label=r"$L=160$")
    line2, = ax1.plot(x2, rho2_x[0], label=r"$L=320$")
    line3, = ax1.plot(x3, rho3_x[0], label=r"$L=640$")
    line4, = ax1.plot(x4, rho4_x[0], label=r"$L=1280$")

    ax2.plot(x1, rho1_x[1])
    ax2.plot(x2, rho2_x[1])
    ax2.plot(x3, rho3_x[1])
    ax2.plot(x4, rho4_x[1])


    # mask1 = np.logical_or(x1 < 0.2, x1 >=0.8)
    # rm_1 = np.mean(rho1_x[:, mask1], axis=1)
    # ax1.axhline(rm_1[0], linestyle="dotted", c=line1.get_c())
    # ax2.axhline(rm_1[1], linestyle="dotted", c=line1.get_c())

    # mask2 = np.logical_or(x2 < 0.2, x2 >=0.8)
    # rm_2 = np.mean(rho2_x[:, mask2], axis=1)
    # ax1.axhline(rm_2[0], linestyle="dotted", c=line2.get_c())
    # ax2.axhline(rm_2[1], linestyle="dotted", c=line2.get_c())

    # mask3 = np.logical_or(x3 < 0.2, x3 >=0.8)
    # rm_3 = np.mean(rho3_x[:, mask3], axis=1)
    # ax1.axhline(rm_3[0], linestyle="dotted", c=line3.get_c())
    # ax2.axhline(rm_3[1], linestyle="dotted", c=line3.get_c())


    ax1.axvline(0.335, linestyle="dashed", c="tab:grey")
    ax1.axvline(0.665, linestyle="dashed", c="tab:grey")
    ax2.axvline(0.335, linestyle="dashed", c="tab:grey")
    ax2.axvline(0.665, linestyle="dashed", c="tab:grey")
    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$x/L$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")
    ax1.legend()
    plt.show()
    plt.close()



def SP_profiles_varied_L():
    f1 = "data/time_ave_fields/L160_160_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    f2 = "data/time_ave_fields/L240_240_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    f3 = "data/time_ave_fields/L320_320_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"


    with np.load(f1, "r") as data:
        rho1, mx1, my1, v1, vx1, vy1 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f2, "r") as data:
        rho2, mx2, my2, v2, vx2, vy2 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f3, "r") as data:
        rho3, mx3, my3, v3, vx3, vy3 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    plt.imshow(rho1[0], origin="lower")
    plt.show()
    plt.close()

    dx = 1
    x1 = np.arange(160) * dx + dx / 2
    x2 = np.arange(240) * dx + dx / 2
    x3 = np.arange(320) * dx + dx / 2

    x1 /= 160
    x2 /= 240
    x3 /= 320

    rho1_x = np.mean(rho1[:, 75:85], axis=1)
    rho2_x = np.mean(rho2[:, 112:128], axis=1)
    rho3_x = np.mean(rho3[:, 150:170], axis=1)

    # rho1_x = np.mean(rho1[:, :], axis=1)
    # rho2_x = np.mean(rho2[:, :], axis=1)
    # rho3_x = np.mean(rho3[:, :], axis=1)
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)

    line1, = ax1.plot(x1, rho1_x[0], label=r"$L=160$")
    line2, = ax1.plot(x2, rho2_x[0], label=r"$L=240$")
    line3, = ax1.plot(x3, rho3_x[0], label=r"$L=320$")
    ax2.plot(x1, rho1_x[1])
    ax2.plot(x2, rho2_x[1])
    ax2.plot(x3, rho3_x[1])

    # mask1 = np.logical_or(x1 < 0.2, x1 >=0.8)
    # rm_1 = np.mean(rho1_x[:, mask1], axis=1)
    # ax1.axhline(rm_1[0], linestyle="dotted", c=line1.get_c())
    # ax2.axhline(rm_1[1], linestyle="dotted", c=line1.get_c())

    # mask2 = np.logical_or(x2 < 0.2, x2 >=0.8)
    # rm_2 = np.mean(rho2_x[:, mask2], axis=1)
    # ax1.axhline(rm_2[0], linestyle="dotted", c=line2.get_c())
    # ax2.axhline(rm_2[1], linestyle="dotted", c=line2.get_c())

    # mask3 = np.logical_or(x3 < 0.2, x3 >=0.8)
    # rm_3 = np.mean(rho3_x[:, mask3], axis=1)
    # ax1.axhline(rm_3[0], linestyle="dotted", c=line3.get_c())
    # ax2.axhline(rm_3[1], linestyle="dotted", c=line3.get_c())

    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$x/L$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")
    ax1.legend()
    plt.show()
    plt.close()


def CCB_profiles_varied_L_r5():
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    Lx_arr = [320, 640, 960]
    Ly_arr = [320, 640, 320]
    dx = 1

    for i, Lx in enumerate(Lx_arr):
        Ly = Ly_arr[i]
        fname = f"data/time_ave_fields/L{Lx}_{Ly}_Dr0.100_k0.70_p2.5_1.25_r5_5_5_e-2.000_J0.500_-0.500_h0.333_2000.npz"
        with np.load(fname, "r") as data:
            rho, mx, my, v, vx, vy = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

        x = np.arange(Lx) * dx + dx / 2
        x /= Lx
        row1 = int(75 * Ly / 160)
        row2 = int(85 * Ly / 160)
        rho_x = np.mean(rho[:, row1: row2], axis=1)

        line, = ax1.plot(x, rho_x[0], label=r"$L_x=%d, L_y=%d$" % (Lx, Ly))
        ax2.plot(x, rho_x[1])

        mask = np.logical_or(x < 0.2, x >=0.8)
        rm = np.mean(rho_x[:, mask], axis=1)
        ax1.axhline(rm[0], linestyle="dotted", c=line.get_c())
        ax2.axhline(rm[1], linestyle="dotted", c=line.get_c())

    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")
    ax1.legend()
    plt.show()
    plt.close()


def CCB_profiles_large_NR():
    Lx_arr = [320]
    Ly_arr = [320]
    dx = 1
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)

    for i, Lx in enumerate(Lx_arr):
        Ly = Ly_arr[i]
        fname = f"data/time_ave_fields/L{Lx}_{Ly}_Dr0.100_Dt0_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.npz"
        with np.load(fname, "r") as data:
            rho, mx, my, v, vx, vy = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]
        
        # plt.imshow(rho[0], origin="lower")
        # plt.show()
        # plt.close()

        x = np.arange(Lx) * dx + dx / 2
        x /= Lx
        row1 = int(75 * Ly / 160)
        row2 = int(85 * Ly / 160)
        rho_x = np.mean(rho[:, row1: row2], axis=1)

        line, = ax1.plot(x, rho_x[0], label=r"$L_x=%d, L_y=%d$" % (Lx, Ly))
        ax2.plot(x, rho_x[1])

        mask = np.logical_or(x < 0.2, x >=0.8)
        rm = np.mean(rho_x[:, mask], axis=1)
        ax1.axhline(rm[0], linestyle="dotted", c=line.get_c())
        ax2.axhline(rm[1], linestyle="dotted", c=line.get_c())

    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")
    ax1.legend()
    plt.show()
    plt.close()


def CCB_LA_profiles():
    Lx_arr = [160]
    Ly_arr = [160]
    dx = 1
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)

    for i, Lx in enumerate(Lx_arr):
        Ly = Ly_arr[i]
        fname = f"data/time_ave_fields/L{Lx}_{Ly}_Dr0.100_k0.70_p22_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
        with np.load(fname, "r") as data:
            rho, mx, my, v, vx, vy = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]
        
        # plt.imshow(rho[0], origin="lower")
        # plt.show()
        # plt.close()

        x = np.arange(Lx) * dx + dx / 2
        x /= Lx
        row1 = int(75 * Ly / 160)
        row2 = int(85 * Ly / 160)
        rho_x = np.mean(rho[:, row1: row2], axis=1)

        line, = ax1.plot(x, rho_x[0], label=r"$L_x=%d, L_y=%d$" % (Lx, Ly))
        ax2.plot(x, rho_x[1])

        mask = np.logical_or(x < 0.2, x >=0.8)
        rm = np.mean(rho_x[:, mask], axis=1)
        ax1.axhline(rm[0], linestyle="dotted", c=line.get_c())
        ax2.axhline(rm[1], linestyle="dotted", c=line.get_c())

    ax2.set_xlim(0, 1)
    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")
    ax1.legend()
    plt.show()
    plt.close()


def show_instant_fields():
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.333"
    fname = f"{folder}/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"

    Lx, Ly = get_L(os.path.basename(fname))
    dx = 1
    x = np.arange(Lx) * dx + dx / 2
    i_frame = 50

    with fl.open(name=fname, mode="rb") as fgsd:
        rho5, mx5, my5, v5, vx5, vy5 = get_cg_fields(i_frame, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=None)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6), constrained_layout=True)

    im = ax.imshow(rho5[1], origin="lower", vmax=24, vmin=2)
    plt.colorbar(im, ax=ax)
    plt.show()
    plt.close()


def show_instant_profile():
    folder = f"{root_tahmineh}/QS_MPI_noself/L960_320_r10_10_10_h0.333"
    fname = f"{folder}/L960_320_Dr0.100_k0.70_p5_2.41_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"

    Lx, Ly = get_L(os.path.basename(fname))
    dx = 1
    x = np.arange(Lx) * dx + dx / 2
    i_frame = 37

    with fl.open(name=fname, mode="rb") as fgsd:
        x_c, y_c = find_center(i_frame, fgsd, Lx, Ly)
        center = [x_c, y_c]
        print(center)
        rho5, mx5, my5, v5, vx5, vy5 = get_cg_fields(i_frame, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
    rho5_x = np.mean(rho5, axis=1)

    if Lx == 320:
        delta_x = 80
    elif Lx == 960:
        delta_x = 300
    rho5_gas_1 = np.mean(rho5_x[:, 0:delta_x], axis=1)
    rho5_gas_2 = np.mean(rho5_x[:, Lx-delta_x:], axis=1)
    plt.plot(x, rho5_x[0])
    plt.plot(x, rho5_x[1])

    # plt.axhline(rho5_gas_1[0], linestyle="dotted")
    # plt.axhline(rho5_gas_2[0], linestyle="dotted", c="tab:green")
    # plt.axhline(rho5_gas_3[0], linestyle="dotted", c="tab:orange")

    plt.axhline(rho5_gas_1[1], linestyle="dotted", c="tab:pink")
    plt.axhline(rho5_gas_2[1], linestyle="dotted", c="tab:green")

    plt.axhline(rho5_gas_1[0], linestyle="dashed", c="tab:pink")
    plt.axhline(rho5_gas_2[0], linestyle="dashed", c="tab:green")

    plt.axvline(Lx/2, c="tab:green", linestyle="dashed")
    plt.axvline(delta_x, c="tab:green", linestyle="dashed")
    plt.axvline(Lx-delta_x, c="tab:green", linestyle="dashed")


    plt.xlim(0, Lx)
    plt.show()
    plt.close()

    beg_frame = 0
    with fl.open(name=fname, mode="rb") as fgsd:
        x_c, y_c = find_center(i_frame, fgsd, Lx, Ly)
        center = [x_c, y_c]

        r1, r2 = np.zeros([2, fgsd.nframes-beg_frame, 2])
        for i in range(beg_frame, fgsd.nframes):
            rho5, mx5, my5, v5, vx5, vy5 = get_cg_fields(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
            rho5_x = np.mean(rho5, axis=1)
            print("frame", i)

            r1[i-beg_frame] = np.mean(rho5_x[:, 0:delta_x], axis=1)
            r2[i-beg_frame] = np.mean(rho5_x[:, Lx-delta_x:], axis=1)
        t = np.arange(beg_frame, fgsd.nframes)
    dr = r1 - r2
    plt.plot(t, dr[:, 0], "o")
    plt.plot(t, dr[:, 1], "s")
    plt.show()
    plt.close()


def update_time_ave_fields():

    fname = f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_varied_Dt_h/L320_320_Dr0.100_Dt0_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_2000.gsd"
    beg_frame = 30
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/L160_r10_10_10_h0.333/L160_160_Dr0.100_k0.70_p22_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 46
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/MC_L160_r10_10_10_h0.333/L160_160_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 20
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/MC_L240_r10_10_10_h0.333/L240_240_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 10
    time_ave(fname, beg_frame, center_type="B")


    fname = f"{root_tahmineh}/QS_MPI_noself/MC_L320_r10_10_10_h0.333/L320_320_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 8
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/MC_L1280_r10_10_10_h0.333/L1280_1280_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 1
    time_ave(fname, beg_frame, center_type="B")


    '''
    fname = f"{root_tahmineh}/QS_MPI_noself/L160_r10_10_10_h0.333/L160_160_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 200
    time_ave(fname, beg_frame)

    fname = f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_h0.333/L320_320_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 100
    time_ave(fname, beg_frame)

    fname = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333/L640_640_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 50
    time_ave(fname, beg_frame)

    fname = f"{root_tahmineh}/QS_MPI_noself/L640_r5_5_5_h0.333/L640_640_Dr0.100_k0.70_p2.5_1.25_r5_5_5_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 20
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/L320_r5_5_5_h0.333/L320_320_Dr0.100_k0.70_p2.5_1.25_r5_5_5_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 30
    time_ave(fname, beg_frame, center_type="B")

    fname = f"{root_tahmineh}/QS_MPI_noself/L160_r5_5_5_h0.333/L160_160_Dr0.100_k0.70_p2.5_1.25_r5_5_5_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 0
    time_ave(fname, beg_frame, center_type="B")
 
    fname = f"{root_tahmineh}/QS_MPI_noself/L960_320_r5_5_5_h0.333/L960_320_Dr0.100_k0.70_p2.5_1.25_r5_5_5_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 30
    time_ave(fname, beg_frame, center_type="B")
    '''


if __name__ == "__main__":
    ### get time series of gas density around the CCB domains
    """
    with fl.open(name=f3, mode="rb") as fgsd:
        x_c, y_c = find_center(0, fgsd, Lx, Ly)
        nframes = fgsd.nframes

        rho_m = np.zeros((4, nframes, 2))
        t = np.zeros(nframes)
        for i in range(nframes):
            rho5, mx5, my5, v5, vx5, vy5 = get_cg_fields(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=[-168, 0])
            rho5_x = np.mean(rho5, axis=1)
            rho_m[0, i] = np.mean(rho5_x[:, 90:190], axis=1)
            rho_m[1, i] = np.mean(rho5_x[:, 450:550], axis=1)
            rho_m[2, i] = np.mean(rho5_x[:, 720:820], axis=1)
            rho_m[3, i] = np.mean(rho5_x[:, 1100:1200], axis=1)
            t[i] = i * 1e4
            print(i)
        fout = "data/tmp.npz"
        np.savez_compressed(fout, rho_m=rho_m, t=t)
    """

    ### show the time series of gas density around the CCB domains
    
    # with np.load("data/tmp.npz", "r") as data:
    #     t, rho_m = data["t"], data["rho_m"]

    # drho = (rho_m[2, :, 1] + rho_m[3, :, 1] - rho_m[1, :, 1] - rho_m[0, :, 1]) / 2
    # plt.plot(t, drho, "-o")
    # plt.show()
    # plt.close()


    # update_time_ave_fields()
    # CCB_profiles_varied_L_r5()

    # with np.load("data/time_ave_fields/L320_320_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz", "r") as data:
    #     rho = data["rho"]
    #     plt.imshow(rho[1], origin="lower")
    #     plt.show()
    #     plt.close()

    # CCB_profiles_varied_L()
    # SP_profiles_varied_L()

    # show_instant_profile()
    # show_instant_fields()
    # CCB_profiles_large_NR()
    CCB_LA_profiles()

    # fig, ax = plt.subplots(1, 1)
    # ax.plot([160, 320, 640], [4, 3.8, 3.5], "--o")
    # # ax.set_xscale("log")
    # # ax.set_yscale("log")

    # # from add_line import add_line
    # # add_line(ax, 0, 1, 1, slope=-1, scale="log")
    # plt.show()
    # plt.close()