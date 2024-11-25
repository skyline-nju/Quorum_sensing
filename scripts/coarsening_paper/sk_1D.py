import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
sys.path.append("../../")

from add_line import add_line

root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_Lx(fname):
    s = os.path.basename(fname).rstrip(".npz").split("_")
    Lx = int(s[0].lstrip("L"))
    return Lx


def get_q1(q_radial, Sq):
    """ Get the first moment of s(q, t). """
    I0 = np.sum(Sq)
    I1 = np.sum(q_radial * Sq)
    q1 = I1 / I0
    return q1


def orientation_average(S, q_radial):
    n = S.size // 2
    S_new = np.zeros_like(q_radial)
    for i in range(q_radial.size):
        S_new[i] = 0.5 * (S[n+1+i] + S[n-1-i])
    return S_new


def cal_Sq_t_rho(f_fields, dx=1, dt=1/3., refresh_all=False):
    folder, basename = os.path.split(f_fields)
    fout = f"data/Sqt_rho_1D/{basename}"

    if os.path.exists(fout) and os.path.getmtime(fout) > os.path.getmtime(f_fields):
        with np.load(fout) as data:
            t, q_radial, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
        print("loading data from", fout)
    else:
        L = get_Lx(f_fields)
        with np.load(f_fields, "r") as data:
            t, x = data["t"] * dt, data["x"]
            fields = data["fields"]
            
        n = int(L / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[1:qx.size//2]
        qx = np.fft.fftshift(qx)
        # qx_ij, qy_ij = np.meshgrid(qx, qy)
        # q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        rho_Sqt = np.zeros((3, t.size, q_radial.size))
        rho_var = np.zeros((3, t.size))

        beg = 0
        if os.path.exists(fout) and not refresh_all:
            with np.load(fout, "r") as data:
                t_pre, q_pre, rho_Sqt_pre, rho_var_pre = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
                beg = t_pre.size
                rho_Sqt[:, :beg, :] = rho_Sqt_pre
                rho_var[:, :beg] = rho_var_pre
                print("loading %d frames from" % beg, fout)
        for j, t_i in enumerate(t[beg:]):
            i = j + beg
            print("processing %d/%d frame, t=%g" % (i, t.size, t_i))
            rhoA_m, rhoB_m = np.mean(fields[i][0]), np.mean(fields[i][1])
            rhoA_q = np.fft.fft(fields[i][0]-rhoA_m, norm="ortho")
            rhoB_q = np.fft.fft(fields[i][1]-rhoB_m, norm="ortho")

            rhoA_q = np.fft.fftshift(rhoA_q)
            rhoB_q = np.fft.fftshift(rhoB_q)

            S_A = np.abs(rhoA_q) ** 2
            S_B = np.abs(rhoB_q) ** 2
            S_tot = np.abs(rhoA_q + rhoB_q) ** 2
            
            rho_Sqt[0][i] = orientation_average(S_A, q_radial)
            rho_Sqt[1][i] = orientation_average(S_B, q_radial)
            rho_Sqt[2][i] = orientation_average(S_tot, q_radial)

            rho_var[0][i] = np.var(fields[i][0])
            rho_var[1][i] = np.var(fields[i][1])
            rho_var[2][i] = np.var(fields[i][0] + fields[i][1])

        np.savez_compressed(fout, t=t, q=q_radial, rho_Sqt=rho_Sqt, rho_var=rho_var)
        print("updating cache:", fout)
    return t, q_radial, rho_Sqt, rho_var


def get_r(fname):
    t, q, rho_Sqt, rho_var = cal_Sq_t_rho(fname)
    rA = np.zeros_like(t)
    rB = np.zeros_like(t)

    for i in range(t.size):
        q1_A = get_q1(q, rho_Sqt[0, i])
        q1_B = get_q1(q, rho_Sqt[1, i])
        rA[i] = np.pi * 2 / q1_A
        rB[i] = np.pi * 2 / q1_B
    return t, rA, rB


def plot_R_t_1D(ax, ms=4):
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1"
    fname = f"{folder}/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"

    folder2 = f"{root_tahmineh}/QS_MPI_noself/L5120_80_r10_10_10_h0.333/profile_1D_dx1"
    fname2 = f"{folder2}/L5120_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    t, rA, rB = get_r(fname)
    # ax.plot(t, rA, "o")
    ax.plot(t, rB, "o", fillstyle="none", label=r"$L_x=1280$", ms=ms)

    t, rA, rB = get_r(fname2)
    print(t)
    print(t.size)
    # ax.plot(t, rA, "o")
    ax.plot(t, rB, "s", fillstyle="none", label=r"$L_x=5120$", ms=ms)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize="x-large")

    add_line(ax, 0.1, 0.45, 0.95, 0.13, label="0.13", xl=0.6, yl=0.6, fontsize="large")
    # add_line(ax, 0.1, 0.5, 0.95, 0.33, label="0.33", xl=0.8, yl=0.5, fontsize="large")
    add_line(ax, 0.1, 0.3, 0.95, 1, label="1", xl=0.15, yl=0.7, fontsize="large")

if __name__ == "__main__":
    # rho = np.random.rand(12)
    # rho_q = np.fft.fft(rho)

    # qx = np.fft.fftfreq(12, d=1/(2 * np.pi))

    # q_radial = qx[1:qx.size//2]

    # print(rho_q)
    # print(qx)
    # print(q_radial)
    folder = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/profile_1D_dx1"
    fname = f"{folder}/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"


    folder2 = f"{root_tahmineh}/QS_MPI_noself/L5120_80_r10_10_10_h0.333/profile_1D_dx1"
    fname2 = f"{folder2}/L5120_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"



    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    t, rA, rB = get_r(fname)
    ax.plot(t, rA, "o")
    ax.plot(t, rB, "s")

    t, rA, rB = get_r(fname2)
    print(t)
    print(t.size)
    ax.plot(t, rA, "o")
    ax.plot(t, rB, "s")

    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.1, 0.5, 0.95, 0.13, label="0.13", xl=0.8, yl=0.5, fontsize="large")
    # add_line(ax, 0.1, 0.5, 0.95, 0.33, label="0.33", xl=0.8, yl=0.5, fontsize="large")
    add_line(ax, 0.1, 0.3, 0.95, 1, label="1", xl=0.2, yl=0.6, fontsize="large")

    plt.show()
    plt.close()