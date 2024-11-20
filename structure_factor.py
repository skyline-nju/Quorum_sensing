import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from corr_length import get_Lx, get_para
from coarse_grain import coarse_grain_all
from add_line import add_line


root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def orientation_average(S, q_radial, q_module):
    """ Cal orientation-averaged structure factor"""
    Sq = np.zeros_like(q_radial)
    half_dq = (q_radial[1]-q_radial[0]) / 2

    for i, q in enumerate(q_radial):
        mask = np.logical_and(q_module >= q-half_dq, q_module < q+half_dq)
        Sq[i] = np.mean(S[mask])
    return Sq


def get_q1(q_radial, Sq):
    """ Get the first moment of s(q, t). """
    I0 = np.sum(Sq)
    I1 = np.sum(q_radial * Sq)
    q1 = I1 / I0
    return q1


def cal_Sq_t_rho(f_fields, dx=1, dt=0.01, refresh_all=False):
    folder, basename = os.path.split(f_fields)
    fout = f"data/Sqt_rho/{basename}"

    if os.path.exists(fout) and os.path.getmtime(fout) > os.path.getmtime(f_fields):
        with np.load(fout) as data:
            t, q_radial, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
        print("loading data from", fout)
    else:
        L = get_Lx(f_fields)
        with np.load(f_fields, "r") as data:
            t, x, y = data["t"] * dt, data["x"], data["y"]
            fields = data["fields"]
            
        n = int(L / dx)
        qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
        q_radial = qx[1:qx.size//2-1]
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)
        qx_ij, qy_ij = np.meshgrid(qx, qy)
        q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)
        rho_Sqt = np.zeros((3, t.size, q_radial.size))
        rho_var = np.zeros((3, t.size))

        beg = 0
        if os.path.exists(fout) and not refresh_all:
            with np.load(fout, "r") as data:
                t_pre, q_pre, rho_Sqt_pre, rho_var_pre = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
                beg = t_pre.size
                rho_Sqt[:, :beg, :] = rho_Sqt_pre
                rho_var[:, :beg] = rho_var_pre
                print("loading %d frames from" % beg, fout)
        for j, t_i in enumerate(t[beg:]):
            i = j + beg
            print("processing %d/%d frame, t=%g" % (i, t.size, t_i))
            rhoA_m, rhoB_m = np.mean(fields[i][0]), np.mean(fields[i][1])
            rhoA_q = np.fft.fft2(fields[i][0]-rhoA_m, norm="ortho")
            rhoB_q = np.fft.fft2(fields[i][1]-rhoB_m, norm="ortho")

            rhoA_q = np.fft.fftshift(rhoA_q)
            rhoB_q = np.fft.fftshift(rhoB_q)

            S_A = np.abs(rhoA_q) ** 2
            S_B = np.abs(rhoB_q) ** 2
            S_tot = np.abs(rhoA_q + rhoB_q) ** 2
            
            rho_Sqt[0][i] = orientation_average(S_A, q_radial, q_module)
            rho_Sqt[1][i] = orientation_average(S_B, q_radial, q_module)
            rho_Sqt[2][i] = orientation_average(S_tot, q_radial, q_module)

            rho_var[0][i] = np.var(fields[i][0])
            rho_var[1][i] = np.var(fields[i][1])
            rho_var[2][i] = np.var(fields[i][0] + fields[i][1])

        np.savez_compressed(fout, t=t, q=q_radial, rho_Sqt=rho_Sqt, rho_var=rho_var)
        print("updating cache:", fout)
    return t, q_radial, rho_Sqt, rho_var


def get_npz_file_name(para):
    prefix = f"data/Sqt_rho/L{para['Lx']:g}_{para['Ly']}_Dr{para['Dr']:.3f}"
    if para["Dt"] > 0:
        prefix = f"{prefix}_Dt{para['Dt']:g}"
    prefix = f"{prefix}_k{para['k']:.2f}_p{para['pA']:g}_{para['pB']:g}_" \
        f"r{para['rho0']:g}_{para['bar_rho_A']}_{para['bar_rho_B']}_" \
        f"e{para['eta']:.3f}_J{para['eta_AB']:.3f}_{-para['eta_AB']:.3f}"
    if para['h'] != 0.01:
        prefix = f"{prefix}_h{para['h']:.3f}"
    fname = f"{prefix}_{para['seed']}.npz"
    return fname


def get_Sq_t_rho(f_fields=None, f_Sqt=None, dx=1, dt=0.01, refresh_all=False):
    if f_fields is not None:
        t, q, rho_Sqt, rho_var = cal_Sq_t_rho(f_fields, dx, dt, refresh_all)
    elif f_Sqt is not None:
        with np.load(f_Sqt, "r") as data:
                t, q, rho_Sqt, rho_var = data["t"], data["q"], data["rho_Sqt"], data["rho_var"]
    return t, q, rho_Sqt, rho_var



def get_h(fname):
    basename = os.path.basename(fname)
    s_list = basename.rstrip(".npz").split("_")
    h = 0.01
    for s in s_list:
        if s[0] == "h":
            h = float(s.lstrip("h"))
        if h == 0.333:
            h = 1 / 3.
    return h


def update_density_Sqt(dx=1):
    prefix_list = [
        # f"{root_tahmineh}/QS_MPI_noself/L320_r10_10_10_h0.333",
        # f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333",
        f"{root_tahmineh}/QS_MPI_noself/L640_r5_5_5",
        f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333",
        # f"{root_tahmineh}/QS_MPI_noself/L1280_r5_5_5_h0.333",
        f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.333",
        # f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.100_Dt0.01",
        # f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.100_Dt0.001",
        # f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_h0.100_Dt0.005",
        # f"{root_tahmineh}/QS_MPI_noself/L2560_r5_5_5_h0.333",
        # f"{root_tahmineh}/QS_MPI_noself/L2560_r10_10_10_h0.333",
        f"{root_tahmineh}/QS_MPI_noself/L1280_r10_10_10_varied_Dt_h"
        # f"{root_tahmineh}/QS_MPI_noself/L160_r40_6_30",
        # f"{root_tahmineh}/QS_MPI_noself/L320_r40_6_30"
    ]
    for prefix in prefix_list:
        L = int(os.path.basename(prefix).lstrip("L").split("_")[0])
        print("L =", L)
        if L == 2560:
            only_rho = True
        else:
            only_rho = False
        coarse_grain_all(prefix, dx=dx, only_rho=only_rho)
        folder = f"{prefix}/coarse_grain_dx{dx:g}"
        fins = glob.glob(f"{folder}/*.npz")
        for fin in fins:
            h = get_h(fin)
            get_Sq_t_rho(f_fields=fin, dx=dx, dt=h)
    
    # prefix_list = [
    #     f"{root_tahmineh}/QS_MPI_noself/L80_Turing",
    #     f"{root_tahmineh}/QS_MPI_noself/L80_Turing_v3_1",
    #     f"{root_tahmineh}/QS_MPI_noself/L640_r10_10_10_h0.333"
    # ]
    # for prefix in prefix_list:
    #     folder = f"{prefix}/coarse_grain_dx{dx:g}"
    #     fins = glob.glob(f"{folder}/*.npz")
    #     for fin in fins:
    #         h = get_h(fin)
    #         get_Sq_t_rho(f_fields=fin, dx=dx, dt=h)


def collapse_rho_Sq(fsqt, frames=[8, 21, 34, 57, 95], q_thresh=None):
    para = get_para(fsqt)
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fsqt, dt=para["h"])
    
    S_rho_A, S_rho_B = rho_Sqt[:2]
    if q_thresh is None:
        fig, axes = plt.subplots(2, 2, figsize=(6.5, 6), sharex="col", sharey="col", constrained_layout=True)
    else:
        fig, axes = plt.subplots(2, 3, figsize=(8, 6), sharex="col", sharey="col", constrained_layout=True)

    ms = 2
    for i in frames:
        axes[0, 0].plot(q, S_rho_A[i], "o", fillstyle="none", ms=ms, label=r"$%g$" % t[i])
        axes[1, 0].plot(q, S_rho_B[i], "o", fillstyle="none", ms=ms, label=r"$%g$" % t[i])

        q1_A = get_q1(q, S_rho_A[i])
        q1_B = get_q1(q, S_rho_B[i])
        # axes[0, 1].plot(q / q1_A, S_rho_A[i] * q1_A**2 / rho_var[0][i], "o", fillstyle="none", ms=ms)
        # axes[1, 1].plot(q / q1_B, S_rho_B[i] * q1_B**2 / rho_var[0][i], "o", fillstyle="none", ms=ms)
        axes[0, 1].plot(q / q1_A, S_rho_A[i] * q1_A**2, "o", fillstyle="none", ms=ms)
        axes[1, 1].plot(q / q1_B, S_rho_B[i] * q1_B**2, "o", fillstyle="none", ms=ms)
        if q_thresh is not None:
            if not isinstance(q_thresh, list):
                q_thresh = [q_thresh, q_thresh]
            mask_A = q < q_thresh[0]
            mask_B = q < q_thresh[1]
            q1_A = get_q1(q[mask_A], S_rho_A[i][mask_A])
            q1_B = get_q1(q[mask_B], S_rho_B[i][mask_B])
            axes[0, 2].plot(q / q1_A, S_rho_A[i] * q1_A**2, "o", fillstyle="none", ms=ms)
            axes[1, 2].plot(q / q1_B, S_rho_B[i] * q1_B**2, "o", fillstyle="none", ms=ms)
    
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")

    axes[0, 0].set_ylabel(r"$S_A(q, t)$")
    axes[1, 0].set_ylabel(r"$S_B(q, t)$")

    axes[0, 1].set_ylabel(r"$S_A(q,t)\times (q^1_\mathrm{A,mix})^2$")
    axes[1, 1].set_ylabel(r"$S_B(q,t)\times (q^1_\mathrm{B,mix})^2$")



    axes[0, 0].set_xlabel(r"$q$")
    axes[0, 1].set_xlabel(r"$q / q^1_\mathrm{A,mix}(t)$")
    axes[1, 0].set_xlabel(r"$q$")
    axes[1, 1].set_xlabel(r"$q / q^1_\mathrm{B,mix}(t)$")

    axes[0, 0].legend(title=r"$t=$",fontsize="small")
    if q_thresh is not None:
        axes[0, 2].set_xscale("log")
        axes[0, 2].set_yscale("log")
        axes[0, 2].set_ylabel(r"$S_A(q,t)\times (q^1_\mathrm{A,L})^2$")
        axes[1, 2].set_ylabel(r"$S_B(q,t)\times (q^1_\mathrm{B,L})^2$")
        axes[0, 2].set_xlabel(r"$q / q^1_\mathrm{A,L}(t)$")
        axes[1, 2].set_xlabel(r"$q / q^1_\mathrm{B,L}(t)$")

    title = r"$L=%g, \eta=%g,\eta_{AB}=-\eta_{BA}=%g, D_r=%g, \rho_0=\bar{\rho}_{A,B}=%g,\phi_A=%g,\phi_B=%g$" % (
        para["Lx"], para["eta"], para["eta_AB"], para["Dr"], para["rho0"], para["pA"], para["pB"]
    )
    fig.suptitle(title)
    plt.show()
    plt.close()


def h_one_third_varied_L():
    para = {}
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["Dt"] = 0
    para["eta"] = -2
    para["eta_AB"] = 0.5

    para["h"] = 1 / 3.

    para["Lx"] = 1280
    para["Ly"] = para["Lx"]
    para["pA"] = 5.5
    para["pB"] = 3.85

    composition_pairs = [(3, 6), (5, 6), (5.5, 3.85), (6, 3)]
    L_list = [[1280], [1280], [640, 1280, 2560], [1280]]
    ms = 3.5
    titles = [r"(a) $\phi_A=%g,\phi_B=%g$",
              r"(b) $\phi_A=%g,\phi_B=%g$",
              r"(c) $\phi_A=%g,\phi_B=%g$",
              r"(d) $\phi_A=%g,\phi_B=%g$"]
    fig, axes = plt.subplots(2, 2, figsize=(8, 8.2), sharex=True, sharey=True, constrained_layout=True)
    for j, ax in enumerate(axes.flat):
        para["pA"], para["pB"] = composition_pairs[j]
        for L in L_list[j]:
            para["Lx"] = para["Ly"] = L
            f_Sqt = get_npz_file_name(para)
            t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
            q1_A, q1_B = np.zeros((2, t.size))
            for i in range(t.size):
                q1_A[i] = get_q1(q, rho_Sqt[0][i])
                q1_B[i] = get_q1(q, rho_Sqt[1][i])
            line, = ax.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms)
            ax.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c())
        ax.set_title(titles[j] % (para['pA'], para['pB']))
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")

    add_line(axes[0, 0], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(axes[0, 1], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(axes[0, 1], 0.1, 0.4, 0.95, 0.25, label="0.25", xl=0.40, yl=0.75, fontsize="large")

    add_line(axes[1, 0], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(axes[1, 1], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    plt.show()
    plt.close()


def L1280_varied_h():
    para = {}
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["Dt"] = 0
    para["eta"] = -2
    para["eta_AB"] = 0.5
    para["Lx"] = 1280
    para["Ly"] = para["Lx"]

    para["h"] = 1 / 3.
    para["pA"] = 5.5
    para["pB"] = 3.85

    composition_pairs = [(5, 6), (5.5, 3.85), (6, 3)]
    h_list = [[0.01, 0.1, 1/3.], [0.01, 0.1, 1/3.], [0.01, 0.1, 1/3.]]
    ms = 3.5
    titles = [r"(a) $\phi_A=%g,\phi_B=%g$",
              r"(b) $\phi_A=%g,\phi_B=%g$",
              r"(c) $\phi_A=%g,\phi_B=%g$",
              r"(d) $\phi_A=%g,\phi_B=%g$"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharex=True, sharey=True, constrained_layout=True)
    for j, ax in enumerate(axes.flat):
        para["pA"], para["pB"] = composition_pairs[j]
        for h in h_list[j]:
            para['h'] = h
            f_Sqt = get_npz_file_name(para)
            t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
            q1_A, q1_B = np.zeros((2, t.size))
            for i in range(t.size):
                q1_A[i] = get_q1(q, rho_Sqt[0][i])
                q1_B[i] = get_q1(q, rho_Sqt[1][i])
            line, = ax.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms)
            ax.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c())
        ax.set_title(titles[j] % (para['pA'], para['pB']))
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")

    add_line(axes[0], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(axes[1], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(axes[2], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    plt.show()
    plt.close()


def L1280_varied_Dt():
    para = {}
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["eta"] = -2
    para["eta_AB"] = 0.5
    para["Lx"] = 1280
    para["Ly"] = para["Lx"]

    para["h"] = 0.1
    para["pA"] = 5.5
    para["pB"] = 3.85
    para["Dt"] = 0
    Dt_arr = np.array([0, 0.001, 0.005, 0.01])

    ms = 3.5

    fig, ax = plt.subplots(1, 1, figsize=(4, 4.2), sharex=True, sharey=True, constrained_layout=True)

    for para["Dt"] in Dt_arr:
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])

        line, = ax.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:D_t=%g$" % para["Dt"])
        ax.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:D_t=%g$" % para["Dt"])

    # ax.set_title(titles[j] % (para['pA'], para['pB']))
    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax, 0.5, 0.55, 0.98, 0.14, label="0.14", xl=0.5, yl=0.58, fontsize="large")
    plt.show()
    plt.close()


def L1280_varied_composition():
    para = {}
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["eta"] = -2
    para["eta_AB"] = 0.5
    para["Lx"] = 1280
    para["Ly"] = para["Lx"]

    para["h"] = 1 / 3.
    para["pA"] = 5.5
    para["pB"] = 3.85
    para["Dt"] = 0
    ms = 3.5

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4.2), sharex=True, sharey=True, constrained_layout=True)
    composition_pairs = [(5.5, 3.85), (7.5, 10)]
    para["rho0"] = para["bar_rho_A"] = para["bar_rho_B"] = 10
    for (para["pA"], para["pB"]) in composition_pairs:
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        line, = ax1.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:\phi_A=%g,\phi_B=%g$" % (para["pA"], para["pB"]))
        ax1.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:\phi_A=%g,\phi_B=%g$" % (para["pA"], para["pB"]))
    
    # composition_pairs = [(3.75, 5), (7.5, 10)]
    # r0_arr = [5, 10]
    # for i, (para["pA"], para["pB"]) in enumerate(composition_pairs):
    #     para["rho0"] = para["bar_rho_A"] = para["bar_rho_B"] = r0_arr[i]
    #     f_Sqt = get_npz_file_name(para)
    #     t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
    #     q1_A, q1_B = np.zeros((2, t.size))
    #     for i in range(t.size):
    #         q1_A[i] = get_q1(q, rho_Sqt[0][i])
    #         q1_B[i] = get_q1(q, rho_Sqt[1][i])
    #     line, = ax2.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:\rho_0=%g$" % (para["rho0"]))
    #     ax2.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:\rho_0=%g$" % (para["rho0"]))
    composition_pairs = [(3.75, 5), (4.8, 8.075), (5, 8.7)]
    r0_arr = [5, 5, 5]
    for i, (para["pA"], para["pB"]) in enumerate(composition_pairs):
        para["rho0"] = para["bar_rho_A"] = para["bar_rho_B"] = r0_arr[i]
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        line, = ax2.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:\rho_0=%g$" % (para["rho0"]))
        ax2.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:\rho_0=%g$" % (para["rho0"]))
    para["pA"] = 3.75
    para["pB"] = 5
    para["rho0"] = para["bar_rho_A"] = para["bar_rho_B"] = 5
    L_arr = [1280, 2560]
    for i, para["Lx"] in enumerate(L_arr):
        para["Ly"] = para["Lx"]
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        line, = ax3.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:\rho_0=%g$" % (para["rho0"]))
        ax3.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:\rho_0=%g$" % (para["rho0"]))
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    add_line(ax1, 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax2, 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax2, 0.1, 0.4, 0.95, 0.25, label="0.25", xl=0.60, yl=0.75, fontsize="large")
    add_line(ax2, 0.1, 0.4, 0.95, 0.14, label="0.14", xl=0.80, yl=0.75, fontsize="large")
    add_line(ax2, 0.1, 0.05, 0.95, 0.33, label="0.33", xl=0.80, yl=0.75, fontsize="large")
    add_line(ax2, 0.4, 0.05, 0.95, 0.33, label="0.33", xl=0.80, yl=0.75, fontsize="large")

    add_line(ax2, 0.1, 0.01, 0.95, 0.5, label="0.5", xl=0.80, yl=0.75, fontsize="large")

    add_line(ax3, 0.1, 0.01, 0.95, 0.33, label="0.33", xl=0.40, yl=0.4, fontsize="large")
    add_line(ax3, 0.4, 0.1, 0.95, 0.13, label="0.13", xl=0.40, yl=0.05, fontsize="large")
    add_line(ax3, 0.5, 0.4, 0.95, 0.2, label="0.2", xl=0.40, yl=0.3, fontsize="large")
    # add_line(ax3, 0.6, 0.4, 0.95, 0.2, label="0.2", xl=0.40, yl=0.3, fontsize="large")


    # add_line(ax3, 0.1, 0.1, 0.95, 0.14, label="0.14", xl=0.80, yl=0.75, fontsize="large")

    # add_line(axes[1], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    # add_line(axes[2], 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    plt.show()
    plt.close()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.2), constrained_layout=True)
    composition_pairs = [(5.5, 3.85), (7.5, 10)]
    para["rho0"] = para["bar_rho_A"] = para["bar_rho_B"] = 10
    para["Lx"] = para["Ly"] = 1280

    for i, (para["pA"], para["pB"]) in enumerate(composition_pairs):
        line, = ax1.plot(para["pA"], para["pB"], "o", label=r"$\phi_A=%g,\phi_B=%g$" % (para["pA"], para["pB"]))
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
        ax2.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_{\rm mix}$", c=line.get_c())
        mask = t > 600
        t_new = t[mask]
        q1_L_A, q1_S_A = np.zeros((2, t_new.size))
        # q1_L_B, q1_S_B = np.zeros((2, t_new.size))
        q_thresh = 0.25
        mask_L = q < q_thresh
        mask_S = q >= q_thresh
        for i in range(t_new.size):
            q1_L_A[i] = get_q1(q[mask_L], rho_Sqt[0][i+t.size-t_new.size][mask_L])
            q1_S_A[i] = get_q1(q[mask_S], rho_Sqt[0][i+t.size-t_new.size][mask_S])
            # q1_L_B[i] = get_q1(q[mask_L], rho_Sqt[1][i+t.size-t_new.size][mask_L])
            # q1_S_B[i] = get_q1(q[mask_S], rho_Sqt[1][i+t.size-t_new.size][mask_S])
        ax2.plot(t_new, 2 * np.pi / q1_L_A, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_L$")
        ax2.plot(t_new, 2 * np.pi / q1_S_A, "v", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_S$")
        # ax2.plot(t_new, 2 * np.pi / q1_L_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_L$")
        # ax2.plot(t_new, 2 * np.pi / q1_S_B, "v", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_S$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax1.set_title(r"(a)")
    ax2.set_title(r"(b) $L=1280$")
    ax1.set_xlabel(r"$\phi_A$")
    ax1.set_ylabel(r"$\phi_B$")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$R_A$")

    ax1.plot([3.8935, 10.5367], [1.2001, 21.6776], linestyle="dotted", c="tab:grey")
    ax1.plot(3.8935, 1.2001, "v", fillstyle="none", c="tab:grey", label="liquid binodal")
    ax1.plot(10.5367, 21.6776, "^", fillstyle="none", c="tab:brown", label="gas binodal")
    ax1.legend()
    ax2.legend()

    fig.suptitle(r"$D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=0.5,\rho_0=\bar{\rho}_{A,B}=10$")

    add_line(ax2, 0.4, 0.4, 0.95, 0.13, label="0.13", xl=0.70, yl=0.45, fontsize="large")
    add_line(ax2, 0.4, 0.73, 0.95, 0.10, label="0.10", xl=0.70, yl=0.75, fontsize="large")
    plt.show()
    plt.close()



def SP_S1D1_S2D1():
    para = {}
    para["rho0"] = 10
    para["bar_rho_A"] = 10
    para["bar_rho_B"] = 10
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["Dt"] = 0
    para["eta"] = -2
    para["eta_AB"] = 0.5

    para["h"] = 1 / 3.

    para["Lx"] = 1280
    para["Ly"] = para["Lx"]
    para["pA"] = 5.5
    para["pB"] = 3.85

    composition_pairs = [(3, 6), (5.5, 3.85), (6, 3)]
    L_list = [[1280], [1280, 2560], [1280]]
    ms = 3.5
    titles = [r"(a) $\phi_A=%g,\phi_B=%g$",
              r"(b) $\phi_A=%g,\phi_B=%g$",
              r"(c) $\phi_A=%g,\phi_B=%g$"]
    fig, axes = plt.subplots(2, 3, figsize=(10, 7), sharex="col", sharey=True, constrained_layout=True)
    for col in range(3):
        para["pA"], para["pB"] = composition_pairs[col]
        for L in L_list[col]:
            para["Lx"] = para["Ly"] = L
            f_Sqt = get_npz_file_name(para)
            t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
            if col == 1:
                print(t[-1])
                mask = t <= 1e7
                t = t[mask]
                rho_Sqt = rho_Sqt[:, mask]
            q1_A, q1_B = np.zeros((2, t.size))
            for i in range(t.size):
                q1_A[i] = get_q1(q, rho_Sqt[0][i])
                q1_B[i] = get_q1(q, rho_Sqt[1][i])
            if col >= 1:
                line, = axes[0, col].plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_{\rm mix}:L=%d$" % L)
                axes[1, col].plot(t, 2 * np.pi/q1_B, "o", fillstyle="none", ms=ms, label=r"$R_{\rm mix}:L=%d$" % L)
                mask = t > 600
                t_new = t[mask]
                q1_L_A, q1_S_A = np.zeros((2, t_new.size))
                q1_L_B, q1_S_B = np.zeros((2, t_new.size))
                
                q_thresh = 0.25
                mask_L = q < q_thresh
                mask_S = q >= q_thresh
                for i in range(t_new.size):
                    q1_L_A[i] = get_q1(q[mask_L], rho_Sqt[0][i+t.size-t_new.size][mask_L])
                    q1_S_A[i] = get_q1(q[mask_S], rho_Sqt[0][i+t.size-t_new.size][mask_S])
                    q1_L_B[i] = get_q1(q[mask_L], rho_Sqt[1][i+t.size-t_new.size][mask_L])
                    q1_S_B[i] = get_q1(q[mask_S], rho_Sqt[1][i+t.size-t_new.size][mask_S])
                
                # if col == 1:
                #     mask = t_new <= 3e6
                #     t_new = t_new[mask]
                #     q1_L_A = q1_L_A[mask]
                #     q1_S_A = q1_S_A[mask]
                #     q1_L_B = q1_L_B[mask]
                #     q1_S_B = q1_S_B[mask]
                axes[0, col].plot(t_new, 2 * np.pi / q1_L_A, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_L:L=%d$" % L)
                axes[0, col].plot(t_new, 2 * np.pi / q1_S_A, "v", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_S:L=%d$" % L)
                axes[1, col].plot(t_new, 2 * np.pi / q1_L_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_L:L=%d$" % L)
                axes[1, col].plot(t_new, 2 * np.pi / q1_S_B, "v", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_S:L=%d$" % L)
            else:
                axes[0, col].plot(t, 2 * np.pi/q1_A, "s", fillstyle="none", ms=ms, label=r"$L=%d$" % L)
                axes[1, col].plot(t, 2 * np.pi/q1_B, "s", fillstyle="none", ms=ms, label=r"$L=%d$" % L)

        axes[0, col].set_title(titles[col] % (para['pA'], para['pB']))
        axes[1, col].set_xlabel(r"$t$")
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_xscale("log")
    axes[0, 2].set_xscale("log")

    axes[0, 0].set_ylabel(r"$R_A$")
    axes[1, 0].set_ylabel(r"$R_B$")
    axes[0, 0].legend()
    axes[1, 1].legend(loc=(0.55, 0.2), fontsize="small")
    axes[1, 2].legend(loc=(0.55, 0.2), fontsize="small")

    # axes[0, 1].axvline(600, linestyle="dashed", c="tab:red")
    fig.suptitle(r"$D_r=0.1,\eta=-2,\eta_{AB}=-\eta_{BA}=0.5,\rho_0=\bar{\rho}_{A,B}=10$")
    add_line(axes[0, 0], 0.3, 0.4, 0.95, 0.33, label="0.33", xl=0.50, yl=0.7, fontsize="large")
    add_line(axes[1, 0], 0.3, 0.4, 0.95, 0.33, label="0.33", xl=0.50, yl=0.7, fontsize="large")

    add_line(axes[0, 1], 0.4, 0.4, 0.95, 0.13, label="0.13", xl=0.70, yl=0.45, fontsize="large")
    add_line(axes[0, 1], 0.4, 0.70, 0.95, 0.10, label="0.10", xl=0.70, yl=0.72, fontsize="large")
    add_line(axes[1, 1], 0.4, 0.58, 0.95, 0.13, label="0.13", xl=0.70, yl=0.65, fontsize="large")
    add_line(axes[1, 1], 0.4, 0.82, 0.95, 0.10, label="0.10", xl=0.6, yl=0.9, fontsize="large")

    add_line(axes[0, 2], 0.4, 0.4, 0.95, 0.13, label="0.13", xl=0.70, yl=0.45, fontsize="large")
    add_line(axes[0, 2], 0.4, 0.82, 0.95, 0.10, label="0.10", xl=0.6, yl=0.9, fontsize="large")
    add_line(axes[0, 2], 0.4, 0.6, 0.95, 0.14, label="0.14", xl=0.6, yl=0.7, fontsize="large")

    add_line(axes[1, 2], 0.4, 0.6, 0.95, 0.13, label="0.13", xl=0.70, yl=0.65, fontsize="large")
    add_line(axes[1, 2], 0.4, 0.82, 0.95, 0.10, label="0.10", xl=0.6, yl=0.9, fontsize="large")

    plt.show()
    plt.close()


def cluster_phase_separation():
    para = {}
    para["rho0"] = 40
    para["bar_rho_A"] = 6
    para["bar_rho_B"] = 30
    para["seed"] = 1000
    para["Dr"] = 0.1
    para["k"] = 0.7
    para["eta"] = -2
    para["eta_AB"] = 2
    para["Lx"] = 160
    para["Ly"] = para["Lx"]

    para["h"] = 1/3.
    para["pA"] = 4
    para["pB"] = 46
    para["Dt"] = 0

    ms = 3.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.2), sharex=True, sharey=True, constrained_layout=True)

    # ax1
    for para["Lx"] in [160, 320]:
        para["Ly"] = para["Lx"]
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])

        line, = ax1.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$L=%g$" % para["Lx"])
        ax1.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$L=%g$" % para["Lx"])
        # if para["Lx"] == 160:
        #     line, = ax2.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:D_t=%g$" % para["Dt"])
        #     ax2.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:D_t=%g$" % para["Dt"])

    
    # ax2
    para["pA"] = 5
    para["pB"] = 34
    para["h"] = 1 / 3.
    for para["Lx"] in [160, 320]:
        para["Ly"] = para["Lx"]
        f_Sqt = get_npz_file_name(para)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=para['h'])
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])

        line, = ax2.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$R_A:D_t=%g$" % para["Dt"])
        ax2.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, c=line.get_c(), label=r"$R_B:D_t=%g$" % para["Dt"])

    # ax.set_title(titles[j] % (para['pA'], para['pB']))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend()

    add_line(ax1, 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax1, 0.5, 0.55, 0.98, 0.25, label="0.25", xl=0.5, yl=0.58, fontsize="large")
    add_line(ax1, 0.5, 0.55, 0.98, 0.14, label="0.14", xl=0.5, yl=0.58, fontsize="large")

    add_line(ax2, 0.1, 0.4, 0.95, 0.33, label="0.33", xl=0.40, yl=0.75, fontsize="large")
    add_line(ax2, 0.5, 0.55, 0.98, 0.25, label="0.25", xl=0.5, yl=0.58, fontsize="large")
    add_line(ax2, 0.5, 0.55, 0.98, 0.14, label="0.14", xl=0.5, yl=0.58, fontsize="large")
    plt.show()
    plt.close()


def fig2():
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(6, 4), layout="constrained")
    
    gs1 = gridspec.GridSpec(2, 6, figure=fig)
    ax1 = fig.add_subplot(gs1[0, :2])
    ax2 = fig.add_subplot(gs1[1, :2])
    ax3 = fig.add_subplot(gs1[:, 2:])

    f0= "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f0, dt=1/3.)
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])

    # frame 40 t=1e5
    ax1.plot(q, rho_Sqt[0, 40], "o", ms=2, fillstyle="none", c="tab:grey")
    ax2.plot(q, rho_Sqt[1, 40], "o", ms=2, fillstyle="none", c="tab:grey")
    ax3.plot(t, 2 * np.pi/q1_B, "o", fillstyle="none", ms=3, c="tab:grey", label="SPS")

    f1 = "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f1, dt=1/3.)
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])

    # frame 40 t=1e5
    ax1.plot(q, rho_Sqt[0, 40], "o", ms=2, fillstyle="none")
    ax2.plot(q, rho_Sqt[1, 40], "o", ms=2, fillstyle="none")
    line1, = ax3.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=3, label="DCSP")

    f2 = "data/Sqt_rho/L640_640_Dr0.100_k0.70_p20_20_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f2, dt=1/3.)
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        print(t[i])
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
    # frame 40 t=1e5
    ax1.plot(q, rho_Sqt[0, 30], "o", ms=2, fillstyle="none")
    ax2.plot(q, rho_Sqt[1, 30], "o", ms=2, fillstyle="none")
    ax3.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=3, label="DC")

    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")

    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax3.set_xscale("log")

    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax3.set_yscale("log")

    # ax1.set_xlabel(r"$q$")
    ax2.set_xlabel(r"$q$")
    ax3.set_xlabel(r"$t$")

    ax1.set_ylabel(r"$S_q(A)$")
    ax2.set_ylabel(r"$S_q(B)$")
    ax3.set_ylabel(r"$R_1$")

    ax3.legend(fontsize="large")
    add_line(ax3, 0.1, 0.2, 0.95, 1/3, label="0.33", yl=0.6)
    add_line(ax3, 0.3, 0.4, 0.95, 0.13, label="0.13", xl=0.7, yl=0.5, c=line1.get_c())

    plt.show()
    plt.close()


def CCB_G_large_eta():
    Lx = 1280
    h = 0.1
    ms = 3.5

    f_Sqt = "data/Sqt_rho/L1280_1280_Dr0.100_Dt0_k0.70_p7_4_r10_10_10_e-2.000_J2.000_-2.000_h0.100_1000.npz"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=h)
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    ax.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$L=%g$" % Lx)
    ax.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, label=r"$L=%g$" % Lx)
    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax, 0.1, 0.4, 0.95, 0.13, label="0.13", xl=0.15, yl=0.5, fontsize="large")
    add_line(ax, 0.1, 0., 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")

    plt.show()
    plt.close()


def CCB_G_small_phi_A():
    Lx = 1280
    h = 0.1
    ms = 3.5

    f_Sqt = "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p4.65_5.6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
    t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=f_Sqt, dt=h)
    q1_A, q1_B = np.zeros((2, t.size))
    for i in range(t.size):
        q1_A[i] = get_q1(q, rho_Sqt[0][i])
        q1_B[i] = get_q1(q, rho_Sqt[1][i])

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    ax.plot(t, 2 * np.pi/q1_A, "o", fillstyle="none", ms=ms, label=r"$L=%g$" % Lx)
    ax.plot(t, 2 * np.pi/q1_B, "^", fillstyle="none", ms=ms, label=r"$L=%g$" % Lx)
    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax, 0.1, 0.4, 0.95, 0.13, label="0.13", xl=0.15, yl=0.5, fontsize="large")
    add_line(ax, 0.1, 0., 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")

    plt.show()
    plt.close()

def CCB_LA_varied_L():
    phi_A, phi_B = 22, 6
    # phi_A, phi_B = 7.5, 3
    fins = [
        f"data/Sqt_rho/L320_320_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
        f"data/Sqt_rho/L640_640_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
        f"data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p{phi_A:g}_{phi_B:g}_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
        ]
    
    Lx = [320, 640, 1280]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)
    for j, fin in enumerate(fins):
        h = get_h(fin)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fin, dt=h)
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
        ax1.plot(t, 2 * np.pi /q1_A, "o", fillstyle="none", label=r"$L=%g$" % Lx[j])
        ax2.plot(t, 2 * np.pi /q1_B, "o", fillstyle="none", label=r"$L=%g$" % Lx[j])

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\zeta_A$", fontsize="large")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax2.set_ylabel(r"$\zeta_B$", fontsize="large")
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")
    plt.suptitle(r"$\bar{\rho}_A/\rho_0=%g, \bar{\rho}_B/\rho_0=%g$" % (phi_A/10, phi_B/10), fontsize="x-large")

    add_line(ax1, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.4, 0.95, 0.13, label="0.13", xl=0.15, yl=0.5, fontsize="large")
    add_line(ax2, 0.1, 0.1, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax2, 0.1, 0.4, 0.95, 0.13, label="0.13", xl=0.15, yl=0.5, fontsize="large")

    plt.show()
    plt.close()


def CCB_LB_G():
    fins = ["data/Sqt_rho/L640_640_Dr0.100_k0.70_p4.5_14.7_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p20_20_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)

    phi_A = [4.5, 10, 20]
    phi_B = [10, 28, 20]
    for j, fin in enumerate(fins):
        h = get_h(fin)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fin, dt=h)
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
        ax1.plot(t, 2 * np.pi /q1_A, "o", fillstyle="none", label=r"$(%g, %g)$" % (phi_A[j], phi_B[j]))
        ax2.plot(t, 2 * np.pi /q1_B, "o", fillstyle="none", label=r"$(%g, %g)$" % (phi_A[j], phi_B[j]))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\zeta_A$", fontsize="large")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax2.set_ylabel(r"$\zeta_B$", fontsize="large")
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")

    add_line(ax1, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.2, 0.95, 0.5, label="0.5", xl=0.15, yl=0.2, fontsize="large")


    plt.show()
    plt.close()


def NJP_coarsening():
    fins = ["data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p4.5_14.7_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p4.2_7.8_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p17_22.5_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            # "data/Sqt_rho/L640_640_Dr0.100_k0.70_p20_20_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"]
    # fins = ["data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
    #         "data/Sqt_rho/L640_640_Dr0.100_k0.70_p10_28_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
    #         "data/Sqt_rho/L640_640_Dr0.100_k0.70_p20_20_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
    #         "data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p3_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)

    labels = ["G+CCB", "G+LB+CCB", "G+LB+CCB", "LB+CCB", "CCB", "G+LB"]
    # phi_A = [4.5, 10, 20]
    # phi_B = [10, 28, 20]
    for j, fin in enumerate(fins):
        h = get_h(fin)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fin, dt=h)
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
        ax1.plot(t, 2 * np.pi /q1_A, "o", fillstyle="none", label="%s" % labels[j])
        ax2.plot(t, 2 * np.pi /q1_B, "o", fillstyle="none", label="%s" % labels[j])
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\zeta_A$", fontsize="large")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax2.set_ylabel(r"$\zeta_B$", fontsize="large")
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")

    add_line(ax1, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.2, 0.95, 0.5, label="0.5", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.2, 0.95, 1, label="1", xl=0.15, yl=0.2, fontsize="large")



    plt.show()
    plt.close()


def three_CCB():
    fins = ["data/Sqt_rho/L1280_1280_Dr0.100_k0.70_p5.5_3.85_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p22_6_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz",
            "data/Sqt_rho/L640_640_Dr0.100_k0.70_p7.5_3_r10_10_10_e-2.000_J0.500_-0.500_h0.333_1000.npz"
            ]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True, constrained_layout=True)

    labels = ["G+CCB", "LA+CCB", "G+LA+CCB"]
    # phi_A = [4.5, 10, 20]
    # phi_B = [10, 28, 20]
    for j, fin in enumerate(fins):
        h = get_h(fin)
        t, q, rho_Sqt, rho_var = get_Sq_t_rho(f_Sqt=fin, dt=h)
        q1_A, q1_B = np.zeros((2, t.size))
        for i in range(t.size):
            q1_A[i] = get_q1(q, rho_Sqt[0][i])
            q1_B[i] = get_q1(q, rho_Sqt[1][i])
    
        ax1.plot(t, 2 * np.pi /q1_A, "o", fillstyle="none", label="%s" % labels[j])
        ax2.plot(t, 2 * np.pi /q1_B, "o", fillstyle="none", label="%s" % labels[j])
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax1.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\zeta_A$", fontsize="large")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax2.set_ylabel(r"$\zeta_B$", fontsize="large")
    ax1.legend(fontsize="large")
    ax2.legend(fontsize="large")

    add_line(ax1, 0.1, 0.2, 0.95, 0.33, label="0.33", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.2, 0.95, 0.5, label="0.5", xl=0.15, yl=0.2, fontsize="large")
    add_line(ax1, 0.1, 0.2, 0.95, 1, label="1", xl=0.15, yl=0.2, fontsize="large")



    plt.show()
    plt.close()
if __name__ == "__main__":
    # update_density_Sqt()
    # cluster_phase_separation()
    # h_one_third_varied_L()
    # L1280_varied_h()
    # L1280_varied_Dt()
    L1280_varied_composition()

    # SP_S1D1_S2D1()

    # CCB_G_large_eta()

    # CCB_LA_varied_L()

    # CCB_G_small_phi_A()

    # three_CCB()

    # NJP_coarsening()