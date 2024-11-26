import os
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_frames
from gsd import hoomd


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_polarity(snap: hoomd.Snapshot, dt=0.1):
    t = snap.configuration.step * dt
    theta = snap.particles.position[:, 2]

    px = np.mean(np.cos(theta))
    py = np.mean(np.sin(theta))

    p = np.sqrt(px **2 + py ** 2)
    return t, p, px, py


def get_moment_vel(snap: hoomd.Snapshot, dt=0.1):
    t = snap.configuration.step * dt
    theta = snap.particles.position[:, 2]
    v = snap.particles.charge

    px = np.mean(np.cos(theta))
    py = np.mean(np.sin(theta))

    mx = np.mean(np.cos(theta) * v)
    my = np.mean(np.sin(theta) * v)
    return t, px, py, mx, my


if __name__ == "__main__":
    # basename = "L480_160_Dr0.100_k0.70_p10_23_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_h0.333/{basename}"


    # frames = read_frames(fname)

    # fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    # t_arr, p_arr, px_arr, py_arr = [], [], [], []
    # for frame in frames:
    #     t, p, px, py = get_polarity(frame)
    #     t_arr.append(t)
    #     p_arr.append(p)
    #     px_arr.append(px)
    #     py_arr.append(py)
    
    # ax.plot(t_arr, px_arr)
    # ax.plot(t_arr, py_arr)
    # ax.set_xlabel("time", fontsize="large")
    # ax.set_ylabel("global polarity", fontsize="large")

    # plt.show()
    # plt.close()

    # basename = "L120_40_Dr0.100_k0.70_p20_44_r20_20_20_e-2.000_J0.500_-0.500_h0.100_1000.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20_varied_h/{basename}"
    basename = "L480_160_Dr0.100_k0.70_p10_23_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    fname = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_h0.333/{basename}"

    frames = read_frames(fname)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    t_arr, px_arr, py_arr, mx_arr, my_arr = [], [], [], [], []
    for frame in frames:
        t, px, py, mx, my = get_moment_vel(frame)
        t_arr.append(t)
        px_arr.append(px)
        py_arr.append(py)
        mx_arr.append(mx)
        my_arr.append(my)    
    ax.plot(t_arr, px_arr, label=r"$\langle \cos(\theta_i)\rangle_i$")
    ax.plot(t_arr, mx_arr, label=r"$\langle v_i\cos(\theta_i)\rangle_i$")
    ax.set_xlabel("time", fontsize="x-large")
    ax.set_ylabel("global order", fontsize="x-large")
    ax.legend(fontsize="x-large")
    plt.show()
    plt.close()


