import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_one_frame
from gsd import fl
from MSD import find_center

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
    s = fname.lstrip("L").split("_")
    Lx = int(s[0])
    Ly = int(s[1])
    return Lx, Ly


def time_ave():
    basename = "L160_160_Dr0.100_k0.70_p4_4_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    fname = f"{root_tahmineh}/QS_MPI_noself/MC_L160_r10_10_10_h0.333/{basename}"
    Lx, Ly = get_L(basename)
    dx = 1
    rho, mx, my, v, vx, vy = np.zeros((6, 2, Ly, Lx))
    beg_frame = 20
    with fl.open(name=fname, mode="rb") as fgsd:
        count = 0
        for i in range(beg_frame, fgsd.nframes):
            # x_c, y_c = find_center(i, fgsd, Lx, Ly)
            center = None
            rho_cur, mx_cur, my_cur, v_cur, vx_cur, vy_cur = get_cg_fields(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
            rho += rho_cur
            mx += mx_cur
            my += my_cur
            v += v_cur
            vx += vx_cur
            vy += vy_cur
            count +=1
            print("frame:%d/%d" % (i, fgsd.nframes))
    rho /= count
    mx /= count
    my /= count
    v /= count
    vx /= count
    vy /= count

    fout = f"data/time_ave_fields/{basename.replace('.gsd', '.npz')}"
    np.savez_compressed(fout, rho=rho, mx=mx, my=my, v=v, vx=vx, vy=vy)


def varied_system_size():
    f1 = "data/time_ave_fields/L160_160_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    f2 = "data/time_ave_fields/L320_320_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    f3 = "data/time_ave_fields/L640_640_Dr0.100_k0.70_p5_2.96_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"


    with np.load(f1, "r") as data:
        rho1, mx1, my1, v1, vx1, vy1 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f2, "r") as data:
        rho2, mx2, my2, v2, vx2, vy2 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    with np.load(f3, "r") as data:
        rho3, mx3, my3, v3, vx3, vy3 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    plt.imshow(rho1[0], origin="lower")
    plt.show()
    plt.close()
    # Lx = 
    # Ly = 320
    dx = 1
    # x = np.arange(Lx) * dx + dx / 2

    rho1_x = np.mean(rho1[:, 75:85], axis=1)
    x1 = np.arange(160) * dx + dx / 2

    rho2_x = np.mean(rho2[:, 150:170], axis=1)
    x2 = np.arange(320) * dx + dx / 2

    rho3_x = np.mean(rho3[:, 300:340], axis=1)
    x3 = np.arange(640) * dx + dx / 2

    plt.plot(x1, rho1_x[0])
    plt.plot(x2/2, rho2_x[0])
    plt.plot(x3/4, rho3_x[0])

    # plt.xlim(0, 160)
    plt.show()
    plt.close()



def show_rect_snap():
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_320_r10_10_10_h0.333"
    f3 = f"{folder}/L1280_320_Dr0.100_k0.70_p5_2.41_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"

    Lx, Ly = 1280, 320
    dx = 1
    x = np.arange(Lx) * dx + dx / 2

    with fl.open(name=f3, mode="rb") as fgsd:
        x_c, y_c = find_center(0, fgsd, Lx, Ly)
        rho5, mx5, my5, v5, vx5, vy5 = get_cg_fields(26, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=[-168, 0])
    rho5_x = np.mean(rho5, axis=1)


    rho5_gas_1 = np.mean(rho5_x[:, 90:190], axis=1)
    rho5_gas_2 = np.mean(rho5_x[:, 450:550], axis=1)
    rho5_gas_3 = np.mean(rho5_x[:, 720:820], axis=1)
    rho5_gas_4 = np.mean(rho5_x[:, 1100:1200], axis=1)
    plt.plot(x, rho5_x[0])
    plt.plot(x, rho5_x[1])

    # plt.axhline(rho5_gas_1[0], linestyle="dotted")
    # plt.axhline(rho5_gas_2[0], linestyle="dotted", c="tab:green")
    # plt.axhline(rho5_gas_3[0], linestyle="dotted", c="tab:orange")

    plt.axhline(rho5_gas_1[1], linestyle="dotted")
    plt.axhline(rho5_gas_2[1], linestyle="dotted", c="tab:green")
    plt.axhline(rho5_gas_3[1], linestyle="dotted", c="tab:orange")
    plt.axhline(rho5_gas_4[1], linestyle="dotted", c="tab:red")

    plt.axvline(90, c="tab:green", linestyle="dashed")
    plt.axvline(190, c="tab:green", linestyle="dashed")
    plt.axvline(450, c="tab:green", linestyle="dashed")
    plt.axvline(550, c="tab:green", linestyle="dashed")

    plt.axvline(720, c="tab:green", linestyle="dashed")
    plt.axvline(820, c="tab:green", linestyle="dashed")
    plt.axvline(1100, c="tab:green", linestyle="dashed")
    plt.axvline(1200, c="tab:green", linestyle="dashed")

    drho = (rho5_gas_3[1] + rho5_gas_4[1] - rho5_gas_1[1] - rho5_gas_2[1]) / 2
    print(drho * 280 * 320)


    rho5_gas_1 = np.mean(rho5_x[:, 500:750], axis=1)
    rho5_gas_2 = (np.sum(rho5_x[:, 1180:], axis=1) + np.sum(rho5_x[:, :150], axis=1))/250
    # plt.axhline(rho5_gas_1[1], linestyle="dotted", c="tab:blue")
    # plt.axhline(rho5_gas_2[1], linestyle="dotted", c="tab:red")
    print((rho5_gas_2[1] - rho5_gas_1[1]) * 280 * 320) 


    plt.axvline(320)
    plt.axvline(960)
    plt.axvline(640)

    plt.xlim(0, Lx)
    plt.show()
    plt.close()
    

if __name__ == "__main__":
    # f1 = "data/time_ave_fields/L640_320_Dr0.100_k0.70_p4.775_1.9825_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    # f2 = "data/time_ave_fields/L640_320_Dr0.100_k0.70_p5.225_2.8375_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"

    # with np.load(f1, "r") as data:
    #     rho1, mx1, my1, v1, vx1, vy1 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    # with np.load(f2, "r") as data:
    #     rho2, mx2, my2, v2, vx2, vy2 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    # Lx = 640
    # Ly = 320
    # dx = 1
    # x = np.arange(Lx) * dx + dx / 2

    # rho1_x = np.mean(rho1[:, :, :], axis=1)
    # rho2_x = np.mean(rho2[:, :, :], axis=1)



    # folder = f"{root_tahmineh}/QS_MPI_noself/L640_320_r10_10_10_h0.333"
    # f3 = f"{folder}/L640_320_Dr0.100_k0.70_p4.775_1.9825_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    # f4 = f"{folder}/L640_320_Dr0.100_k0.70_p5.225_2.8375_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    # folder = f"{root_tahmineh}/QS_MPI_noself/L1280_320_r10_10_10_h0.333"
    # f3 = f"{folder}/L640_320_Dr0.100_k0.70_p4.77375_1.98125_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    # f4 = f"{folder}/L640_320_Dr0.100_k0.70_p5.22625_2.83875_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"

    # with fl.open(name=f3, mode="rb") as fgsd:
    #     x_c, y_c = find_center(0, fgsd, Lx, Ly)
    #     rho3, mx3, my3, v3, vx3, vy3 = get_cg_fields(0, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=[x_c, y_c])
    # with fl.open(name=f4, mode="rb") as fgsd:
    #     x_c, y_c = find_center(0, fgsd, Lx, Ly)
    #     rho4, mx4, my4, v4, vx4, vy4 = get_cg_fields(0, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=[x_c, y_c])
    # rho3_x = np.mean(rho3, axis=1)
    # rho4_x = np.mean(rho4, axis=1)

    # plt.xlim(0, Lx*2)
    # plt.plot(x, rho1_x[0])
    # plt.plot(x, rho1_x[1])
    # plt.plot(x + Lx, rho2_x[0])
    # plt.plot(x + Lx, rho2_x[1])
    # plt.plot(x, rho3_x[0])
    # plt.plot(x, rho3_x[1])
    # plt.plot(x, rho4_x[0])
    # plt.plot(x, rho4_x[1])

    # rho3_gas = (np.mean(rho3_x[:, 0:140], axis=1) + np.mean(rho3_x[:, 500:], axis=1))/2
    # plt.axhline(rho3_gas[0], 0, 0.5, linestyle="dashed")
    # plt.axhline(rho3_gas[1], 0, 0.5, linestyle="dashed")

    # rho4_gas = (np.mean(rho4_x[:, 0:140], axis=1) + np.mean(rho4_x[:, 500:], axis=1))/2
    # plt.axhline(rho4_gas[0], 0.5, 1, linestyle="dashed")
    # plt.axhline(rho4_gas[1], 0.5, 1, linestyle="dashed")
    # print((rho4_gas[1] - rho3_gas[1]) * 280 * 320)

    # rho1_gas = (np.mean(rho1_x[:, 0:200], axis=1) + np.mean(rho1_x[:, 440:], axis=1))/2
    # plt.axhline(rho1_gas[0], 0, 0.5, linestyle="dotted")
    # plt.axhline(rho1_gas[1], 0, 0.5, linestyle="dotted")

    # rho2_gas = (np.mean(rho2_x[:, 0:200], axis=1) + np.mean(rho2_x[:, 440:], axis=1))/2
    # plt.axhline(rho2_gas[0], 0.5, 1, linestyle="dotted")
    # plt.axhline(rho2_gas[1], 0.5, 1, linestyle="dotted")
    # plt.show()
    # plt.close()

    
    
    
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

    # time_ave()


    f1 = "data/time_ave_fields/L160_160_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    with np.load(f1, "r") as data:
        rho1, mx1, my1, v1, vx1, vy1 = data["rho"], data["mx"], data["my"], data["v"], data["vx"], data["vy"]

    Lx = 160
    Ly = 160
    dx = 1
    x = np.arange(Lx) * dx + dx / 2

    plt.imshow(rho1[1], origin="lower")
    # plt.axhline(180, c="w")
    # plt.axhline(220, c="w")
    plt.show()
    plt.close()

    # rho1_x = np.mean(rho1[:, 180:220], axis=1)
    # plt.plot(x, rho1_x[0])
    # plt.plot(x, rho1_x[1])
    # plt.plot(x+Lx, rho1_x[0])
    # plt.plot(x+Lx, rho1_x[1])

    # rho_m1 = np.mean(rho1_x[:, 400:600], axis=1)
    # rho_m2 = np.mean(rho1_x[:, 600:800], axis=1)

    # rho_m3 = np.mean(rho1_x[:, 20:200], axis=1)
    # rho_m4 = np.mean(rho1_x[:, 1100:1280], axis=1)

    # plt.axhline(rho_m1[0], linestyle="dashed", c="tab:blue")
    # plt.axhline(rho_m2[0], linestyle="dotted", c="tab:blue")

    # plt.axhline(rho_m3[0], linestyle="dashed", c="tab:red")
    # plt.axhline(rho_m4[0], linestyle="dotted", c="tab:red")

    # plt.axhline(rho_m1[1], linestyle="dashed", c="tab:blue")
    # plt.axhline(rho_m2[1], linestyle="dotted", c="tab:blue")

    # plt.axhline(rho_m3[1], linestyle="dashed", c="tab:red")
    # plt.axhline(rho_m4[1], linestyle="dotted", c="tab:red")

    # plt.show()
    # plt.close()