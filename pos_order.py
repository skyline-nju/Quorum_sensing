import numpy as np
import matplotlib.pyplot as plt
from add_line import add_line

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def cal_pos_order(rho_A, rho_B):
    ''' J = rho_A \grad rho_B - rho_B \grad rho_A
    '''
    def diff(f):
        dx_f = 0.5 * (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1))
        dy_f = 0.5 * (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0))
        return dx_f, dy_f

    dx_rho_A, dy_rho_A = diff(rho_A)
    dx_rho_B, dy_rho_B = diff(rho_B)
    Jx = rho_A * dx_rho_B - rho_B * dx_rho_A
    Jy = rho_A * dy_rho_B - rho_B * dy_rho_A
    Jm = np.sqrt(np.mean(Jx) ** 2 + np.mean(Jy) ** 2)
    return Jm


def get_Jt(fnpz):
    with np.load(fnpz, "rb") as data:
        t, fields = data["t"], data["fields"]
        n_frames = t.size
        J = np.zeros(n_frames)
        for i in range(t.size):
            rho_A, rho_B = fields[i, 0:2]
            J[i] = cal_pos_order(rho_A, rho_B)
    return t, J

if __name__ == "__main__":
    L_arr = [40, 80, 160, 320, 640]
    J_mean1 = []
    J_mean2 = []
    J_mean3 = []


    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), sharex=False, sharey=False, constrained_layout=True)
    for L in L_arr:
        if L == 640:
            seed = 1000
        else:
            seed = 1000
        fnpz = f"{root_tahmineh}/QS_MPI_noself/L{L}_r10_10_10_h0.333/coarse_grain_dx1/L{L}_{L}_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"

        t0, J0 = get_Jt(fnpz)
        J0 /= (10 * 28)
        ax1.plot(t0/3, J0, "-o")

        if L <= 160:
            mask = t0 >= 2e5 * 3
            J_mean1.append(np.mean(J0[mask]))
        else:
            mask = t0 >= 7e4
            J_mean1.append(np.mean(J0[mask]))
    
    
    seed = 1000
    for L in L_arr[:4]:
        fnpz = f"{root_tahmineh}/QS_MPI_noself/L{L}_r10_10_10_h0.333/coarse_grain_dx1/L{L}_{L}_Dr0.100_k0.70_p17_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"
        t0, J0 = get_Jt(fnpz)
        J0 /= (17 * 22.5)
        ax2.plot(t0/3, J0, "-o")
        J_mean2.append(np.mean(J0))

    seed = 1000
    for L in L_arr[:5]:
        fnpz = f"{root_tahmineh}/QS_MPI_noself/L{L}_r10_10_10_h0.333/coarse_grain_dx1/L{L}_{L}_Dr0.100_k0.70_p4.5_14.7_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed:d}.npz"
        t0, J0 = get_Jt(fnpz)
        J0 /= (4.5 * 14.7)
        ax3.plot(t0/3, J0, "-o")
        if L < 640:
            J_mean3.append(np.mean(J0))
        else:
            mask = t0 >= 3e5
            J_mean3.append(np.mean(J0[mask]))
    plt.show()
    plt.close()

    fig, ax = plt.subplots(1, 1)
    ax.plot(L_arr[:], J_mean1, "o")
    ax.plot(L_arr[:4], J_mean2, "s")
    # ax.plot(L_arr, J_mean3, ">")

    ax.set_xscale("log")
    ax.set_yscale("log")
    add_line(ax, 0, 1, 1, -1)
    plt.show()

    print(J_mean1)
    print(J_mean2)
    print(J_mean3)
