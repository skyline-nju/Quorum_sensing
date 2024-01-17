import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
from gsd import fl
from add_line import add_line
from MSD import untangle


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

def find_xc(i_frame: int, fgsd: fl.GSDFile, Lx: float, which="both"):
    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    if which == "both":
        x = pos[:,0] / Lx * np.pi * 2
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

    x_c = np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))) / (2 * np.pi) * Lx
    return x_c


def get_xc_series(fin: str, fout: str, beg_frame=0):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            existed_frames = t0.size
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            if beg_frame + existed_frames < nframes:
                Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
                
                t = np.zeros(nframes-beg_frame, int)
                x = np.zeros(nframes-beg_frame)

                if existed_frames > 0:
                    t[:existed_frames] = t0
                    x[:existed_frames] = x0

                for i_frame in range(beg_frame + existed_frames, nframes):
                    try:
                        t[i_frame-beg_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                    except KeyError:
                        t[i_frame-beg_frame] = 0
                        print("Warning, not found t[%d]" % i_frame)
                    print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame-beg_frame]))
                    x[i_frame-beg_frame] = find_xc(i_frame, fgsd, Lx)
                np.savez_compressed(fout, t=t, x=x)
            else:
                print("find no new frames")



if __name__ == "__main__":
    L_arr = [320, 480, 640, 960]

    seed_dict = {320: 1000, 480: 1000, 640: 1000, 960: 2000}
    beg_frame_dict = {320: 30, 480: 260, 640: 200, 960: 10}

    for Lx in L_arr:
        folder = f"{root_tahmineh}/QS_MPI_noself/L{Lx:d}_80_r10_10_10_h0.333"
        basename = f"L{Lx:d}_80_Dr0.100_k0.70_p7.5_10_r10_10_10_e-2.000_J0.500_-0.500_h0.333_{seed_dict[Lx]:d}.gsd"
        fin = f"{folder}/{basename}"
        fout = f"data/center_traj_1D/{basename.replace('.gsd', '.npz')}"

        get_xc_series(fin, fout, beg_frame=beg_frame_dict[Lx])

        with np.load(fout, "r") as data:
            t, xc = data["t"], data["x"]
            xc = untangle(xc, Lx)

    
        plt.plot(t-t[0], xc-xc[0], "-o", label=r"$L_x=%d$" % Lx)
    plt.legend()
    plt.show()
    plt.close()



