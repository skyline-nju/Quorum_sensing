import os
import numpy as np
import matplotlib.pyplot as plt
from MSD_1D import find_xc
from gsd import fl


root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None):
    ncols = int(Lx / dx)
    bin_area = dx * Ly

    par_num = np.zeros((2, ncols), int)
    mx, my, v, vx, vy = np.zeros((5, 2, ncols))

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

    ux = np.cos(pos[: ,2])
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
    return rho, mx, my, v, vx, vy


def get_L(fname):
    s = os.path.basename(fname).lstrip("L").split("_")
    Lx = int(s[0])
    Ly = int(s[1])
    return Lx, Ly


def time_ave(fname, beg_frame=20, center_type="both", dx=1):
    basename = os.path.basename(fname)
    print("handling", basename)
    fout = f"data/time_ave_profiles/{basename.replace('.gsd', '.npz')}"
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
        ncols = int(round(Lx / dx))
        rho, mx, my, v, vx, vy = np.zeros((6, 2, ncols))
        nframes = 0
    with fl.open(name=fname, mode="rb") as fgsd:
        count = 0
        for i in range(beg_frame + nframes, fgsd.nframes):
            if center_type is None:
                center = None
            else:
                center = find_xc(i, fgsd, Lx, which=center_type)
            rho_cur, mx_cur, my_cur, v_cur, vx_cur, vy_cur = get_profiles(i, fgsd, dx=dx, Lx=Lx, Ly=Ly, center=center)
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


def update_time_ave_profiles():
    root = f"{root_tahmineh}/QS_MPI_noself"
    fname = f"{root}/L320_80_r10_10_10_h0.333/L320_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 40
    time_ave(fname, beg_frame)

    root = f"{root_tahmineh}/QS_MPI_noself"
    fname = f"{root}/L960_80_r10_10_10_h0.333/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.gsd"
    beg_frame = 10
    time_ave(fname, beg_frame)

    root = f"{root_tahmineh}/QS_MPI_noself"
    fname = f"{root}/L480_80_r10_10_10_h0.333/L480_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    beg_frame = 264
    time_ave(fname, beg_frame)


def show_time_ave_profiles():
    fnames = ["data/time_ave_profiles/L320_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
              "data/time_ave_profiles/L480_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
              "data/time_ave_profiles/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2000.npz"
    ]
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    for fname in fnames:
        with np.load(fname, "r") as data:
            rho = data["rho"]
            Lx, Ly = get_L(fname)
            dx = 1
            x = np.arange(Lx) * dx + dx / 2
        ax1.plot(x / Lx, rho[0], label=r"$L_x=%d$" % Lx)
        ax2.plot(x / Lx, rho[1])

    ax2.set_xlabel(r"$x/L_x$")
    ax1.set_ylabel(r"$\rho_A$")
    ax2.set_ylabel(r"$\rho_B$")

    ax1.set_xlim(0, 1)
    ax1.legend()
    plt.show()
    plt.close()


def show_instant_profiles():
    # fname = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2129042.gsd"
    # fname = f"{root_tahmineh}/QS_MPI_noself/L1280_80_r10_10_10_h0.333/L1280_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.gsd"
    fname = f"{root_tahmineh}/QS_MPI_noself/L960_80_r10_10_10_h0.333/L960_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_2030.gsd"


    Lx, Ly = get_L(fname)
    dx = 1
    x = np.arange(Lx) * dx + dx / 2

    beg_frame = 0
    end_frame = None
    with fl.open(name=fname, mode="rb") as fgsd:
        if end_frame is None:
            end_frame = fgsd.nframes
        for i_frame in range(beg_frame, end_frame):
            print("frame: %d/%d" % (i_frame, end_frame))
            fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
            rho, mx, my, v, vx, vy = get_profiles(i_frame, fgsd, dx, Lx, Ly, center=None)
            ax1.plot(x, rho[0])
            ax2.plot(x, rho[1])
            ax2.set_xlabel(r"$x/L_x$")
            ax1.set_ylabel(r"$\rho_A$")
            ax2.set_ylabel(r"$\rho_B$")
            plt.show()
            plt.close()
    

def space_time():
    from coarse_grain import coarse_grain_all_1D
    Lx=1280
    seed=1000
    prefix = f"{root_tahmineh}/QS_MPI_noself/L{Lx}_80_r10_10_10_h0.333"
    coarse_grain_all_1D(prefix, 1, False)


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


if __name__ == "__main__":
    # update_time_ave_profiles()
    # show_time_ave_profiles()
    show_instant_profiles()
    # space_time()
