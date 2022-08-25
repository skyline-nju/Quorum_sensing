import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_para, read_one_frame
from gsd import fl


def coarse_grain_one_frame(i_frame: int,
                           fgsd: fl.GSDFile,
                           fields: np.ndarray,
                           dx: float,
                           Lx: float,
                           Ly: float):
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)
    bin_area = dx ** 2

    par_num = np.zeros((2, nrows, ncols), int)
    mx, my, v = np.zeros((3, 2, nrows, ncols))

    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    speed = fgsd.read_chunk(frame=i_frame, name="particles/charge")
    id_type = fgsd.read_chunk(frame=i_frame, name="particles/typeid")

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
    
    fields[i_frame][0:2] = par_num / bin_area
    fields[i_frame][2:4] = mx / bin_area
    fields[i_frame][4:6] = my / bin_area
    mask = par_num > 0
    fields[i_frame][6:8][mask] = v[mask] / par_num[mask]


def coarse_grain(fin: str, fout: str, dx: float):
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
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
            ncols = int(Lx / dx)
            nrows = int(Ly / dx)
            x = np.arange(ncols) * dx + dx / 2
            y = np.arange(nrows) * dx + dx / 2
            
            t = np.zeros(nframes, int)

            # one frame = [rho_A, rho_B, mx_A, mx_B, my_A, my_B, v_A, v_B]
            fields = np.zeros((nframes, 8, nrows, ncols), np.single)

            if existed_frames > 0:
                t[:existed_frames] = t0
                fields[:existed_frames] = fields0

            for i_frame in range(existed_frames, nframes):
                t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                print("processing frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                coarse_grain_one_frame(i_frame, fgsd, fields, dx, Lx, Ly)
            
            fields_name = np.array(["rho_A", "rho_B", "mx_A", "mx_B", "my_A", "my_B", "v_A", "v_B"])
            np.savez_compressed(fout, t=t, x=x, y=y, fields=fields, fields_name=fields_name)


def coarse_grain_all(prefix, dx, forced_updating=False):
    out_folder = f"{prefix}/coarse_grain_dx{dx:g}"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    fins = glob.glob(f"{prefix}/*.gsd")
    for fin in fins:
        basename = os.path.basename(fin)
        out_basename = basename.replace(".gsd", ".npz")
        fout = f"{out_folder}/{out_basename}"
        if os.path.exists(fout):
            if not forced_updating and os.path.getmtime(fout) > os.path.getmtime(fin):
                print("skipping", fin)
                continue
        coarse_grain(fin, fout, dx)
            


if __name__ == "__main__":
    from datetime import datetime
    
    # prefix = "/scratch03.local/yduan/QS_MPI/L480_160"
    prefix = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan/QS_MPI/L320_e-2.0_J0.5"
    # prefix = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan//QS_MPI_noself/L640_r20_20_20"

    coarse_grain_all(prefix, dx=1)

    # basename = "L480_160_Dr0.100_k0.70_p40_40_r40_e-1.000_J1.000_-1.000_1000.gsd"
    # fin = f"{prefix}/{basename}"
    # out_basename = basename.replace(".gsd", ".npz")
    # fout = f"{prefix}/coarse_grain_dx4/{out_basename}"






