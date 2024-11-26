import numpy as np
import matplotlib.pyplot as plt
from gsd import fl


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def untangle(x, Lx):
    x_new = np.copy(x)
    phase_shift = 0 
    for i in range(1, x.size):
        dx = x[i] - x[i-1]
        if dx < -Lx/2:
            phase_shift += Lx
        elif dx > Lx/2:
            phase_shift -= Lx
        x_new[i] += phase_shift
    return x_new

def read_traj(f_traj, Lx, Ly, n_cut=100, which="B"):
    with open(f_traj, "r") as f:
        lines = f.readlines()[n_cut:]
        t, x, y = np.zeros((3, len(lines)))
        for i, line in enumerate(lines):
            s = line.rstrip("\n").split("\t")
            t[i] = float(s[0]) / 3
            if which == "A":
                x[i] = float(s[1])
                y[i] = float(s[2])
            else:
                x[i] = float(s[3])
                y[i] = float(s[4])
        x += Lx/2
        y += Ly/2
        x[x>=Lx/2] -= Lx
        y[y>=Ly/2] -= Ly
        x = untangle(x, Lx)
        y = untangle(y, Ly)
    return t, x, y


def find_center(i_frame: int, fgsd: fl.GSDFile, Lx: float, Ly: float, which="both"):
    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    if which == "both":
        x = pos[:,0] / Lx * np.pi * 2
        y = pos[:,1] / Ly * np.pi * 2
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
        y = pos[mask][:, 1] / Ly * np.pi * 2

    x_c = np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))) / (2 * np.pi) * Lx
    y_c = np.arctan2(np.mean(np.sin(y)), np.mean(np.cos(y))) / (2 * np.pi) * Ly
    return x_c, y_c


if __name__ == "__main__":

    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(6, 6))
    L_arr = [160, 240, 320, 480]
    phiA = 3
    phiB = 6
    for L in L_arr:
        f_traj = f"{root_tahmineh}/QS_MPI_noself/MC_L{L:d}_r10_10_10_h0.333/mc_L{L:d}_{L:d}_Dr0.100_k0.70_p{phiA:g}_{phiB:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.dat"
        t, x, y = read_traj(f_traj, L, L)
        ax.plot(x, y, label=r"$L=%d$" % L)
    ax.legend()
    plt.show()
    plt.close()