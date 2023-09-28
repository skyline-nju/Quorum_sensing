import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from gsd import fl
from structure_factor import orientation_average, get_q1
from add_line import add_line


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def coarse_grain_density_one_frame(i_frame: int,
                                   fgsd: fl.GSDFile,
                                   dx: float,
                                   Lx: float,
                                   Ly: float):
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)
    bin_area = dx ** 2

    par_num = np.zeros((nrows, ncols), int)
    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    x = pos[:,0] + Lx / 2
    y = pos[:,1] + Ly / 2

    col = (x/dx).astype(int)
    row = (y/dx).astype(int)
    col[col < 0] += ncols
    col[col >= ncols] -= ncols
    row[row < 0] += nrows
    row[row >= nrows] -= nrows
    
    for j in range(x.size):
        my_row = row[j]
        my_col = col[j]
        par_num[my_row, my_col] += 1
    return  par_num / bin_area


def coarse_grain_density(fin: str, fout: str, dx: float):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            y0 = data["y"]
            density_fields0 = data["density"]
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
            density_fields = np.zeros((nframes, nrows, ncols), np.single)
            if existed_frames > 0:
                t[:existed_frames] = t0
                density_fields[:existed_frames] = density_fields0

            for i_frame in range(existed_frames, nframes):
                try:
                    t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                except KeyError:
                    t[i_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                density_fields[i_frame] = coarse_grain_density_one_frame(i_frame, fgsd, dx, Lx, Ly)
            np.savez_compressed(fout, t=t, x=x, y=y, density=density_fields)


def coarse_grain_density_all(prefix, dx, forced_updating=False, pat="*.gsd"):
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
        coarse_grain_density(fin, fout, dx)


def cal_Sq_t_rho(density_file, dx, dt, L):
    basename = os.path.basename(density_file)
    fout = f"data/Sqt_rho_IQS/{basename}"

    if os.path.exists(fout) and os.path.getmtime(fout) > os.path.getmtime(density_file):
        with np.load(fout) as data:
            t, q_radial, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
        print("loading data from", fout)
    else:
        with np.load(density_file, "r") as data:
            t, x, y, rho = data["t"] * dt, data["x"], data["y"], data["density"]
        n = int(L / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[1:qx.size//2-1]
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)
        qx_ij, qy_ij = np.meshgrid(qx, qy)
        q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        rho_Sqt = np.zeros((t.size, q_radial.size))
        rho_var = np.zeros(t.size)

        beg = 0
        if os.path.exists(fout):
            with np.load(fout, "r") as data:
                t_pre, q_pre, rho_Sqt_pre, rho_var_pre = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
                beg = t_pre.size
                rho_Sqt[:beg, :] = rho_Sqt_pre
                rho_var[:beg] = rho_var_pre
                print("loading %d frames from" % beg, fout)

        for j, t_i in enumerate(t[beg:]):
            i = j + beg
            print("processing %d/%d frame, t=%g" % (i, t.size, t_i))
            rho_m = np.mean(rho[i])
            rho_q = np.fft.fft2(rho[i]-rho_m, norm="ortho")
            rho_q = np.fft.fftshift(rho_q)
            S_q = np.abs(rho_q) ** 2
            
            rho_Sqt[i] = orientation_average(S_q, q_radial, q_module)
            rho_var[i] = np.var(rho[i])
        np.savez_compressed(fout, t=t, q=q_radial, rho_Sqt=rho_Sqt, rho_var=rho_var)
        print("updating cache:", fout)
    return t, q_radial, rho_Sqt, rho_var    


if __name__ == "__main__":
    # prefix = "/scratch03.local/yduan/IQS_MPI/L320"
    # coarse_grain_density_all(prefix, dx=1)

    Dt = 0
    h = 0.1
    dx = 1
    L = 320
    prefix = f"{root_tahmineh}/IQS_MPI/L{L:d}/coarse_grain_dx1"

    Dr = 0.1
    rho0 = 10
    phi = 10
    eta = -4
    gamma = 2
    seed = 1000

    basename = f"L{L:d}_{L:d}_Dr{Dr:.3f}_Dt{Dt:g}_r{rho0:g}_p{phi:g}_e{eta:.3f}_g{gamma:.3f}_h{h:.3f}_{seed:d}.npz"
    fnpz = f"{prefix}/{basename}"
    # with np.load(fnpz) as data:
    #     t, x, y, rho = data["t"] * h, data["x"], data["y"], data["density"]

    t, q, rho_Sqt, rho_var = cal_Sq_t_rho(fnpz, dx, h, L)

    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    xi = np.zeros(t.size)
    for i in range(t.size):
        xi[i] = 2 * np.pi / get_q1(q, rho_Sqt[i])
    
    ax.plot(t, xi, "o", fillstyle="none")

    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.1, 0.1, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax, 0.2, 0.1, 0.95, 0.25, label="0.25", xl=0.60, yl=0.75, fontsize="large")


    plt.show()
    plt.close()
    


    

