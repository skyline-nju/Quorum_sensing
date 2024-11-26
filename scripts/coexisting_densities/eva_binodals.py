import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("../")
from plot_profile import get_center, shift_profile

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

def get_time_ave_profile(fin, beg_frame, end_frame=None, h=0.1, inversed_A=False, inversed_B=True, shifting_mode="separated"):
    Lx = int(os.path.basename(fin).lstrip("L").split("_")[0])
    with np.load(fin, "r") as data:
        t, x = data["t"] * h, data["x"]
        fields = data["fields"]
    if end_frame is None:
        end_frame = t.size
    n_frames = end_frame - beg_frame
    rho_x_A = np.zeros((n_frames, x.size))
    rho_x_B = np.zeros_like(rho_x_A)
    for i, i_frame in enumerate(np.arange(beg_frame, end_frame)):
        rho_A = fields[i_frame][0]
        rho_B = fields[i_frame][1]
        xm_A = get_center(x, rho_A, Lx, inversed=inversed_A)
        xm_B = get_center(x, rho_B, Lx, inversed=inversed_B)
        if shifting_mode == "A":
            xm_B = xm_A
        elif shifting_mode == "B":
            xm_A = xm_B
        rho_x_A[i] = shift_profile(rho_A, xm_A, Lx)
        rho_x_B[i] = shift_profile(rho_B, xm_B, Lx)
    rA = np.mean(rho_x_A, axis=0)
    rB = np.mean(rho_x_B, axis=0)
    return x, rA, rB


if __name__ == "__main__":
    prefix = f"{root_tahmineh}/QS_MPI_noself/surface_instability/L120_40_r10_MR/profile_1D_dx1"

    basename = "L120_40_Dr1.000_k0.70_p30_14_r10_e-2.000_J0.100_0.100_h0.100_2000.npz"
    fin = f"{prefix}/{basename}"
    beg_frame = 50
    x, rA, rB = get_time_ave_profile(fin, beg_frame, inversed_A=False, inversed_B=False, shifting_mode="A")

    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    plt.plot(x, rA)
    plt.plot(x, rB)

    plt.show()
    plt.close()
