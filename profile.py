import numpy as np
import matplotlib.pyplot as plt


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_center(x, y, Lx, inversed=False):
    theta = x / Lx * np.pi * 2
    px = np.mean(np.cos(theta) * y)
    py = np.mean(np.sin(theta) * y)
    theta_m = np.arctan2(py, px)
    if inversed:
        theta_m += np.pi
    if theta_m < 0:
        theta_m += np.pi * 2
    elif theta_m >= np.pi * 2:
        theta_m -= np.pi * 2
    xm = theta_m / (np.pi * 2) * Lx
    return xm


def shift_profile(y, xc, Lx):
    shift = int(Lx / 2 - np.around(xc))
    y_new = np.roll(y, shift)
    return y_new


def plot_instant_profile(t, x, frame, normed=False):
    fig, axes = plt.subplots(5, 1, figsize=(6, 6), constrained_layout=True, sharex=True)
    
    fields_x = np.array([np.mean(i, axis=0) for i in frame])
    
    if normed:
        fields_x[2] /= fields_x[0]
        fields_x[4] /= fields_x[0]
        fields_x[3] /= fields_x[1]
        fields_x[5] /= fields_x[1]
    
    for i, ax in enumerate(axes[:4]):
        ax.plot(x, fields_x[2*i])
        ax.plot(x, fields_x[2*i+1])

    rv_A = np.mean(frame[0] * frame[6], axis=0)
    rv_B = np.mean(frame[1] * frame[7], axis=0)
    axes[4].plot(x, rv_A)
    axes[4].plot(x, rv_B)

    fig.suptitle(r"$t=%g$" % t)
    plt.show()
    plt.close()


if __name__ == "__main__":
    prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20/coarse_grain_dx1"
    fname = f"{prefix}/L120_40_Dr0.100_k0.70_p5_33_r20_20_20_e-2.000_J2.000_-2.000_1000.npz"

    dt = 0.01
    Lx = 120
    with np.load(fname, "r") as data:
        t, x, y = data["t"] * dt, data["x"], data["y"]
        fields = data["fields"]
    
    i_frame = 103
    rho_A = np.mean(fields[i_frame][0], axis=0)
    rho_B = np.mean(fields[i_frame][1], axis=0)

    lineA, = plt.plot(x, rho_A)
    xm_A = get_center(x, rho_A, Lx)
    plt.axvline(xm_A, linestyle="dashed", c=lineA.get_c())

    lineB, = plt.plot(x, rho_B)
    xm_B = get_center(x, rho_B, Lx, inversed=True)
    plt.axvline(xm_B, linestyle="dashed", c=lineB.get_c())

    plt.title(r"$t=%g$" %t[i_frame])
    plt.show()
    plt.close()

    rho_A_new = shift_profile(rho_A, xm_A, Lx)
    rho_B_new = shift_profile(rho_B, xm_B, Lx)

    plt.plot(x, rho_A_new)
    plt.plot(x, rho_B_new)
    plt.show()
    plt.close()

    plot_instant_profile(t[i_frame], x, fields[i_frame], normed=True)

    i_frames = np.arange(102, 210)
    rho_x_A = np.zeros((i_frames.size, x.size))
    rho_x_B = np.zeros_like(rho_x_A)

    for i, i_frame in enumerate(i_frames):
        rho_A = np.mean(fields[i_frame][0], axis=0)
        rho_B = np.mean(fields[i_frame][1], axis=0)
        xm_A = get_center(x, rho_A, Lx)
        xm_B = get_center(x, rho_B, Lx, inversed=True)
        rho_x_A[i] = shift_profile(rho_A, xm_A, Lx)
        rho_x_B[i] = shift_profile(rho_B, xm_B, Lx)
    
    rA = np.mean(rho_x_A, axis=0)
    rB = np.mean(rho_x_B, axis=0)

    lineA, = plt.plot(x, rA)
    lineB, = plt.plot(x, rB)

    x_edge_1 = [35, 85]
    x_edge_2 = [50, 70]

    plt.axvline(x_edge_1[0], linestyle="dashed", c="tab:grey")
    plt.axvline(x_edge_1[1], linestyle="dashed", c="tab:grey")
    plt.axvline(x_edge_2[0], linestyle="dotted", c="tab:grey")
    plt.axvline(x_edge_2[1], linestyle="dotted", c="tab:grey")

    mask1 = np.logical_or(x <= x_edge_1[0], x > x_edge_1[1])
    r1 = [np.mean(rA[mask1]), np.mean(rB[mask1])]
    mask2 = np.logical_and(x > x_edge_2[0], x <= x_edge_2[1])
    r2 = [np.mean(rA[mask2]), np.mean(rB[mask2])]

    plt.axhline(r1[0], c=lineA.get_c(), linestyle="dashed", label="%.3f" % r1[0])
    plt.axhline(r1[1], c=lineB.get_c(), linestyle="dashed", label="%.3f" % r1[1])

    plt.axhline(r2[0], c=lineA.get_c(), linestyle="dotted", label="%.3f" % r2[0])
    plt.axhline(r2[1], c=lineB.get_c(), linestyle="dotted", label="%.3f" % r2[1])

    plt.legend()

    phiA = 5
    phiB = 33
    
    l1_A = Lx * (phiA - r2[0]) / (r1[0] - r2[0])
    l1_B = Lx * (phiB - r2[1]) / (r1[1] - r2[1])

    print(l1_A, l1_B)

    plt.axvline(l1_A/2, c="k", linestyle="dashdot")

    drhoA = r2[0] - r1[0]
    drhoB = r2[1] - r1[1]

    dl = 20
    phiA_new = phiA + drhoA * dl / Lx
    phiB_new = phiB + drhoB * dl / Lx
    print(phiA_new, phiB_new)

    plt.xlim(0, Lx)
    plt.show()
    plt.close()


