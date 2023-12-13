import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from gsd import fl

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def coarse_grain_one_frame(i_frame: int,
                           fgsd: fl.GSDFile,
                           fields: np.ndarray,
                           dx: float,
                           Lx: float,
                           Ly: float,
                           vxvy: bool = False,
                           only_rho: bool = False):
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)
    bin_area = dx ** 2

    par_num = np.zeros((2, nrows, ncols), int)
    mx, my, v = np.zeros((3, 2, nrows, ncols))

    if vxvy:
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

    x = pos[:,0] + Lx / 2
    y = pos[:,1] + Ly / 2
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
        if vxvy:
            vx[my_id, my_row, my_col] += ux[j] * speed[j]
            vy[my_id, my_row, my_col] += uy[j] * speed[j]
    
    fields[i_frame][0:2] = par_num / bin_area
    if not only_rho:
        fields[i_frame][2:4] = mx / bin_area
        fields[i_frame][4:6] = my / bin_area
        mask = par_num > 0
        fields[i_frame][6:8][mask] = v[mask] / par_num[mask]
        if vxvy:
            # fields[i_frame][8:10][mask] = vx[mask] / par_num[mask]
            # fields[i_frame][10:12][mask] = vy[mask] / par_num[mask]
            fields[i_frame][8:10][mask] = vx[mask]
            fields[i_frame][10:12][mask] = vy[mask]

def coarse_grain(fin: str, fout: str, dx: float, vxvy:bool=False, only_rho:bool=False):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            y0 = data["y"]
            fields0 = data["fields"]
            fields_name0 = data["fields_name"]
            existed_frames = t0.size
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            print(nframes, "frames found for", fin)
            # print(fgsd.find_matching_chunk_names(""))
            # for j in range(fgsd.nframes):
            #     print(j, fgsd.chunk_exists(frame=j, name='particles/charge'))
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
            ncols = int(Lx / dx)
            nrows = int(Ly / dx)
            x = np.arange(ncols) * dx + dx / 2
            y = np.arange(nrows) * dx + dx / 2
            
            t = np.zeros(nframes, int)

            if only_rho:
                fields = np.zeros((nframes, 2, nrows, ncols), np.single)
            elif not vxvy:
                # one frame = [rho_A, rho_B, mx_A, mx_B, my_A, my_B, v_A, v_B]
                fields = np.zeros((nframes, 8, nrows, ncols), np.single)
            else:
                # one frame = [rho_A, rho_B, mx_A, mx_B, my_A, my_B, v_A, v_B, vx_A, vx_B, vy_A, vy_B]
                fields = np.zeros((nframes, 12, nrows, ncols), np.single)

            if existed_frames > 0:
                t[:existed_frames] = t0
                fields[:existed_frames] = fields0

            for i_frame in range(existed_frames, nframes):
                try:
                    t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                except KeyError:
                    t[i_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                coarse_grain_one_frame(i_frame, fgsd, fields, dx, Lx, Ly, vxvy=vxvy, only_rho=only_rho)
            
            if only_rho:
                fields_name = np.array(["rho_A", "rho_B"])
            elif not vxvy:
                fields_name = np.array(["rho_A", "rho_B", "mx_A", "mx_B", "my_A", "my_B", "v_A", "v_B"])
            else:
                fields_name = np.array(["rho_A", "rho_B", "mx_A", "mx_B", "my_A", "my_B", "v_A", "v_B",
                    "vx_A", "vx_B", "vy_A", "vy_B"])

            np.savez_compressed(fout, t=t, x=x, y=y, fields=fields, fields_name=fields_name)


def coarse_grain_all(prefix, dx, forced_updating=False, pat="*.gsd", vxvy=False, only_rho=False):
    out_folder = f"{prefix}/coarse_grain_dx{dx:g}"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    fins = glob.glob(f"{prefix}/{pat}")
    for fin in fins:
        basename = os.path.basename(fin)
        out_basename = basename.replace(".gsd", ".npz")
        fout = f"{out_folder}/{out_basename}"
        if os.path.exists(fout):
            if not forced_updating and os.path.getmtime(fout) > os.path.getmtime(fin):
                print("skipping", fin)
                continue
        coarse_grain(fin, fout, dx, vxvy=vxvy, only_rho=only_rho)
            

def coarse_grain_one_frame_1D(i_frame: int,
                              fgsd: fl.GSDFile,
                              fields: np.ndarray,
                              dx: float,
                              Lx: float,
                              Ly: float):
    ncols = int(Lx / dx)
    bin_area = dx * Ly

    par_num = np.zeros((2, ncols), int)
    mx, my, v = np.zeros((3, 2, ncols))

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

    x = pos[:,0] + Lx / 2
    # y = pos[:,1] + Ly / 2
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
    
    fields[i_frame][0:2] = par_num / bin_area
    fields[i_frame][2:4] = mx / bin_area
    fields[i_frame][4:6] = my / bin_area
    mask = par_num > 0
    fields[i_frame][6:8][mask] = v[mask] / par_num[mask]


def coarse_grain_1D(fin: str, fout: str, dx: float):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            fields0 = data["fields"]
            fields_name0 = data["fields_name"]
            existed_frames = t0.size
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            print(nframes, "frames found for", fin)
            # print(fgsd.find_matching_chunk_names(""))
            # for j in range(fgsd.nframes):
            #     print(j, fgsd.chunk_exists(frame=j, name='particles/charge'))
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
            ncols = int(Lx / dx)
            x = np.arange(ncols) * dx + dx / 2
            t = np.zeros(nframes, int)

            # one frame = [rho_A, rho_B, mx_A, mx_B, my_A, my_B, v_A, v_B]
            fields = np.zeros((nframes, 8, ncols), np.single)

            if existed_frames > 0:
                t[:existed_frames] = t0
                fields[:existed_frames] = fields0

            for i_frame in range(existed_frames, nframes):
                try:
                    t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                except KeyError:
                    t[i_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                coarse_grain_one_frame_1D(i_frame, fgsd, fields, dx, Lx, Ly)
            
            fields_name = np.array(["rho_A", "rho_B", "mx_A", "mx_B", "my_A", "my_B", "v_A", "v_B"])
            np.savez_compressed(fout, t=t, x=x, fields=fields, fields_name=fields_name)


def coarse_grain_all_1D(prefix, dx, forced_updating=False, pat="*.gsd"):
    out_folder = f"{prefix}/profile_1D_dx{dx:g}"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    fins = glob.glob(f"{prefix}/{pat}")
    for fin in fins:
        basename = os.path.basename(fin)
        out_basename = basename.replace(".gsd", ".npz")
        fout = f"{out_folder}/{out_basename}"
        if os.path.exists(fout):
            if not forced_updating and os.path.getmtime(fout) > os.path.getmtime(fin):
                print("skipping", fin)
                continue
        coarse_grain_1D(fin, fout, dx)


if __name__ == "__main__":
    from datetime import datetime

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L80_Turing_v3_1"
    # coarse_grain_all(prefix, dx=1, pat="*.gsd", vxvy=True)
    
    # prefix = f"{root_tahmineh}/QS_MPI_noself/L80_Turing"
    # coarse_grain_all(prefix, dx=1, pat="*.gsd", vxvy=True)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333"
    # coarse_grain_all(prefix, dx=1, pat="*_p20_*.gsd", vxvy=True)


    # prefix = f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_varied_Dt_h"
    # coarse_grain_all(prefix, dx=1, vxvy=True, pat="L320_320_Dr0.100_Dt0_k0.70_p*_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.100_1000.gsd")

    # prefix = f"{root_rudabeh}/QS_MPI_noself/L120_40_r80_12_60_h0.333"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_40_r40_40_40_varied_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r40_40_40_varied_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r20_20_20_varied_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_40_r20_20_20_varied_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L120_40_r10_10_10_varied_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_varied_Dt_h"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L240_160_r10_10_10"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_h0.333"
    # coarse_grain_all(prefix, dx=1, vxvy=True)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L960_160_r10_10_10_h0.333"
    # coarse_grain_all(prefix, dx=1, vxvy=True)

    # prefix = f"{root_sohrab}/QS_not_inc_self/L10_4_v1_1"
    # coarse_grain_all(prefix, pat="*_1002.gsd", dx=0.25, vxvy=True)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L40_16_v1_1"
    # coarse_grain_all(prefix, dx=1, vxvy=True)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L80_32_v1_1"
    # coarse_grain_all(prefix, dx=0.5, vxvy=True)
    # prefix = f"{root_sohrab}/QS_not_inc_self/L30_10_v8_1_J0.25"
    # coarse_grain_all(prefix, dx=0.5, vxvy=True, pat="L30_10_Dr6_12_v8_1_r60*.gsd")

    # prefix = f"{root_sohrab}/QS_not_inc_self/L80_10_v1_1_J0.1"
    # coarse_grain_all(prefix, dx=1, vxvy=True)


    # # prefix = f"{root_sohrab}/QS_not_inc_self/L80_10_v1_1_J0.1/coarse_grain_dx1"
    # prefix = f"{root_sohrab}/QS_not_inc_self/L30_10_v8_1_J0.25/coarse_grain_dx0.5"
    # J_AB = -1
    # Lx = 80
    # # basename = "L80_10_Dr0.100_r20_e-2.000_-2.000_J%.6f_%.6f_h0.100_1000.npz" % (J_AB, -0.25**2/J_AB)
    # basename = "L30_10_Dr6_12_v8_1_r60_e-0.9725_-1.5091_J%.6f_%.6f_h0.025_2000.npz" % ((J_AB, -0.25**2/J_AB))
    # fin = f"{prefix}/{basename}"
    # with np.load(fin, "r") as data:
    #     t = data["t"]
    #     x = data["x"]

    #     rho_A = np.zeros((t.size, x.size))
    #     rho_B = np.zeros((t.size, x.size))

    #     fields = data["fields"]
    #     for i in range(t.size):
    #         rho_A[i] = np.mean(fields[i][0], axis=0)
    #         rho_B[i] = np.mean(fields[i][1], axis=0)

    #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 6), sharex=True, sharey=True, constrained_layout=True)
    #     ax1.imshow(rho_A, origin="lower", extent=[0, Lx, 0, t[-1]], aspect="auto", vmax=80, vmin=20)
    #     ax2.imshow(rho_B, origin="lower", extent=[0, Lx, 0, t[-1]], aspect="auto", vmax=80, vmin=20)
    #     # ax1.axis("a")

    #     ax1.set_xlabel(r"$x$", fontsize="large")
    #     ax2.set_xlabel(r"$x$", fontsize="large")
    #     ax1.set_ylabel(r"$t$", fontsize="x-large")
    #     ax1.set_title(r"(a) $\langle\rho_A(\mathbf{r}, t)\rangle_y$")
    #     ax2.set_title(r"(b) $\langle\rho_B(\mathbf{r}, t)\rangle_y$")

    #     plt.show()
    #     plt.close()
        
        # for i in range(500, 505):
        #     plt.plot(x, rho_A[i])
        #     plt.plot(x, rho_B[i])
        #     plt.show()
        #     plt.close()


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

        ax2.plot(t[60:]/3, np.mean(rhoA[60:, 450:500]/10, axis=1), label="S=A")
        ax2.plot(t[60:]/3, np.mean(rhoB[60:, 450:500]/10, axis=1), label="S=B")

        ax1.set_xlim(1e6)
        ax1.legend()
        ax2.legend()
        ax2.set_xlabel(r"$t$", fontsize="x-large")
        ax1.set_ylabel(r"$\rho_S^{\rm gas}/\rho_0$", fontsize="x-large")
        ax2.set_ylabel(r"$\rho_S^{\rm liquid}/\rho_0$", fontsize="x-large")

        plt.show()
        plt.close()
