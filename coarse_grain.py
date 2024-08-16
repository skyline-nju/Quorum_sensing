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

    prefix = f"{root_tahmineh}/QS_MPI_noself/L40_r10_10_10_h0.333"
    coarse_grain_all(prefix, dx=1, pat="*.gsd", only_rho=True)


    prefix = f"{root_tahmineh}/QS_MPI_noself/L80_r10_10_10_h0.333"
    coarse_grain_all(prefix, dx=1, pat="*.gsd", only_rho=True)

    prefix = f"{root_tahmineh}/QS_MPI_noself/L160_r10_10_10_h0.333"
    coarse_grain_all(prefix, dx=1, pat="*_p10_28_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p17_22.5_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p17_22.5_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p4.5_14.7_*.gsd", only_rho=True)


    prefix = f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_h0.333"
    coarse_grain_all(prefix, dx=1, pat="*_p10_28_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p17_22.5_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p4.5_14.7_*.gsd", only_rho=True)


    prefix = f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333"
    coarse_grain_all(prefix, dx=1, pat="*_p10_28_*.gsd", only_rho=True)
    coarse_grain_all(prefix, dx=1, pat="*_p4.5_14.7_*.gsd", only_rho=True)

    prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_varied_Dt_h"
    coarse_grain_all(prefix, dx=1, vxvy=False, pat="*_p18_14.5_*.gsd")

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r5_5_5"
    # coarse_grain_all_1D(prefix, 1, pat="*_p2.25_7.35_*.gsd")

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L240_160_r10_10_10"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10"
    # coarse_grain_all(prefix, dx=1, vxvy=False)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L480_160_r10_10_10_h0.333"
    # coarse_grain_all(prefix, dx=1, vxvy=True)

    # prefix = f"{root_tahmineh}/QS_MPI_noself/L960_160_r10_10_10_h0.333"
    # coarse_grain_all(prefix, dx=1, vxvy=True)

