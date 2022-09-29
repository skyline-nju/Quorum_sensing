import numpy as np
import matplotlib.pyplot as plt
import os
from corr_length import get_Lx

root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def orientation_average(S, q_radial, q_module):
    """ Cal orientation-averaged structure factor"""
    Sq = np.zeros_like(q_radial)
    half_dq = (q_radial[1]-q_radial[0]) / 2

    for i, q in enumerate(q_radial):
        mask = np.logical_and(q_module >= q-half_dq, q_module < q+half_dq)
        Sq[i] = np.mean(S[mask])
    return Sq


def get_q1(q_radial, Sq):
    """ Get the first moment of s(q, t). """
    I0 = np.sum(Sq)
    I1 = np.sum(q_radial * Sq)
    q1 = I1 / I0
    return q1

def get_Sq_t_rho(fin, dx=1, dt=0.01):
    """ cal domain length Rt for A, B particles respectively"""

    folder, basename = os.path.split(fin)
    fout = f"data/Sqt_rho/{basename}"

    if os.path.exists(fout) and os.path.getmtime(fout) > os.path.getmtime(fin):
        with np.load(fout) as data:
            t, q, R_t_A, R_t_B, R_t_tot = data["t"], data["q"], data["S_rho_A"], data["S_rho_B"], data["S_rho_tot"]
        print("loading data from", fout)
    else:
        L = get_Lx(fin)
        with np.load(fin, "r") as data:
            t, x, y = data["t"] * dt, data["x"], data["y"]
            fields = data["fields"]
            
        n = int(L / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[1:qx.size//2-1]
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)
        qx_ij, qy_ij = np.meshgrid(qx, qy)
        q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        S_rho_A, S_rho_B, S_rho_tot = np.zeros((3, t.size, q_radial.size))
        for i, t_i in enumerate(t):
            print("processing %d/%d frame, t=%g" % (i, t.size, t_i))
            rhoA_m, rhoB_m = np.mean(fields[i][0]), np.mean(fields[i][1])
            rhoA_q = np.fft.fft2(fields[i][0]-rhoA_m, norm="ortho")
            rhoB_q = np.fft.fft2(fields[i][1]-rhoB_m, norm="ortho")

            rhoA_q = np.fft.fftshift(rhoA_q)
            rhoB_q = np.fft.fftshift(rhoB_q)

            S_A = np.abs(rhoA_q) ** 2
            S_B = np.abs(rhoB_q) ** 2
            S_tot = np.abs(rhoA_q + rhoB_q) ** 2
            
            S_rho_A[i] = orientation_average(S_A, q_radial, q_module)
            S_rho_B[i] = orientation_average(S_B, q_radial, q_module)
            S_rho_tot[i] = orientation_average(S_tot, q_radial, q_module)

        np.savez_compressed(fout, t=t, q=q_radial, S_rho_A=S_rho_A, S_rho_B=S_rho_B, S_rho_tot=S_rho_tot)
        print("updating cache:", fout)
    return t, q_radial, S_rho_A, S_rho_B, S_rho_tot


# def collapse_

if __name__ == "__main__":
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10"
    # folder = f"{prefix}/coarse_grain_dx1"
    # fin = f"{folder}/L640_640_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
    # get_Sq_t_rho(fin)


    # prefix = f"{root_rudabeh}/QS_MPI_noself/L1280_r10_10_10"
    # folder = f"{prefix}/coarse_grain_dx1"
    # fin = f"{folder}/L1280_1280_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
    # get_Sq_t_rho(fin)

    fin = "data/Sqt_rho/L640_640_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
    with np.load(fin, "r") as data:
        t, q, S_rho_A, S_rho_B, S_rho_tot = data["t"], data["q"], data["S_rho_A"], data["S_rho_B"], data["S_rho_tot"]
    
    for i, t_i in enumerate(t):
        print(i, t_i)
    # for i, t_i in enumerate(t):
    #     line1, = plt.plot(q, S_rho_A[i])
    #     line2, = plt.plot(q, S_rho_B[i])
    #     plt.plot(q, S_rho_tot[i])
    #     q1_A = get_q1(q, S_rho_A[i])
    #     q1_B = get_q1(q, S_rho_B[i])

    #     plt.axvline(q1_A, linestyle="dashed", c=line1.get_c())
    #     plt.axvline(q1_B, linestyle="dashed", c=line2.get_c())

    #     plt.xscale("log")
    #     plt.yscale("log")
    #     plt.title(r"$t=%g$" % t_i)
    #     plt.show()
    #     plt.close()

    frames = [8, 21, 33, 54, 78]
    fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex="row", sharey="row", constrained_layout=True)
    for i in frames:
        q1_A = get_q1(q, S_rho_A[i])
        q1_B = get_q1(q, S_rho_B[i])

        axes[0, 0].plot(q, S_rho_A[i])
        axes[0, 1].plot(q, S_rho_B[i])
        axes[1, 0].plot(q / q1_A, S_rho_A[i] * q1_A**2)
        axes[1, 1].plot(q / q1_B, S_rho_B[i] * q1_B**2)

    
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_yscale("log")
    plt.show()
    plt.close()

    fin = "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p6_3_r10_10_10_e-2.000_J0.500_-0.500_1000.npz"
    with np.load(fin, "r") as data:
        t, q, S_rho_A, S_rho_B, S_rho_tot = data["t"], data["q"], data["S_rho_A"], data["S_rho_B"], data["S_rho_tot"]
    
    frames = [8, 21, 34, 57, 95]
    fig, axes = plt.subplots(2, 3, figsize=(8, 6), sharex="col", sharey="col", constrained_layout=True)
    ms = 2.5
    for i in frames:
        axes[0, 0].plot(q, S_rho_A[i], "o", fillstyle="none", ms=ms)
        axes[1, 0].plot(q, S_rho_B[i], "o", fillstyle="none", ms=ms)

        q1_A = get_q1(q, S_rho_A[i])
        q1_B = get_q1(q, S_rho_B[i])
        axes[0, 1].plot(q / q1_A, S_rho_A[i] * q1_A**2, "o", fillstyle="none", ms=ms)
        axes[1, 1].plot(q / q1_B, S_rho_B[i] * q1_B**2, "o", fillstyle="none", ms=ms)

        mask = q < 0.3
        q1_A = get_q1(q[mask], S_rho_A[i][mask])
        q1_B = get_q1(q[mask], S_rho_B[i][mask])
        axes[0, 2].plot(q / q1_A, S_rho_A[i] * q1_A**2, "o", fillstyle="none", ms=ms)
        axes[1, 2].plot(q / q1_B, S_rho_B[i] * q1_B**2, "o", fillstyle="none", ms=ms)
    
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")
    axes[0, 2].set_xscale("log")
    axes[0, 2].set_yscale("log")
    plt.show()
    plt.close()