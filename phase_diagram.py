import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from read_gsd import read_one_frame, read_frames


def read_txt(txtfile:str, para_names:list) -> dict:
    data = {}
    with open(txtfile, "r") as f:
        lines = f.readlines()
        n = len(lines)
        for para_name in para_names:
            data[para_name] = np.zeros(n)
        for i, line in enumerate(lines):
            str_list = line.rstrip("\n").split("\t")
            for j, para_name in enumerate(para_names):
                data[para_name][i] = float(str_list[j])
    return data


def get_paras(para_names:list, fname: str)->dict:
    basename = os.path.basename(fname)
    s_list = basename.rstrip('.gsd').split("_")
    paras = {}
    for para_name in para_names:
        for s in s_list:
            if para_name in s:
                s_new = s.lstrip(para_name)
                try:
                    paras[para_name] = float(s_new)
                except ValueError:
                    print("could not convert string to float:", s)
        if para_name not in paras:
            print("Could not get", para_name)
    return paras


def Dr_alpha_eta0():
    alpha = np.array(
        [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4])
    Dr = np.array(
        [0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.2, 0.4, 0.6, 0.8, 1.0, 2, 3])
    alpha_c = np.array(
        [0, 0.3, 0.3, 0.4, 0.4, 0.4, 0.6, 0.9, 1.0, 1.2, 1.4, 3, 4, 4])

    alpha_dict, Dr_dict = {}, {}
    ms = {0: "o", 1: "s"}
    color = {0: "tab:blue", 1: "tab:red"}
    label = {0: "H", 1: "DP"}

    for j in range(Dr.size):
        for i in range(alpha.size):
            if alpha[i] > alpha_c[j]:
                key = 1
            else:
                key = 0
            if key in alpha_dict:
                alpha_dict[key].append(alpha[i])
                Dr_dict[key].append(Dr[j])
            else:
                alpha_dict[key] = [alpha[i]]
                Dr_dict[key] = [Dr[j]]
    for key in sorted(alpha_dict.keys()):
        plt.plot(alpha_dict[key],
                 Dr_dict[key],
                 ms[key],
                 c=color[key],
                 label=label[key])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\tilde{\alpha}$", fontsize="x-large")
    plt.ylabel(r"$D_r$", fontsize="x-large")
    title = r"$L_x=20,L_y=5,\phi_A=\phi_B=\rho_0=100,v_0=1," \
        + r"\tilde{\eta}=0,\kappa=1$"
    plt.title(title)
    plt.legend(loc="lower right", fontsize="x-large")
    plt.tight_layout()
    plt.show()
    plt.close()


def eta_alpha_Dr():
    """
        Phase diagram on (eta, alpha) plane with Dr=0.02
    """
    alpha = np.array([
        0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.4, 2.8,
        3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0
    ])
    alpha_l = {
        0: [0, 0, 0],
        1: [0.4, 0.5],
        2: [0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        3: [1.4, 2.4, 2.8, 3.2, 3.2, 3.6, 4.0, 4.4, 5.6]
    }
    alpha_h = {
        0: [0.3, 0.4, 0.5],
        1: [6, 6],
        2: [1.2, 2, 2.4, 2.8, 2.8, 3.2, 3.6, 4.0, 5.2, 6],
        3: [6, 6, 6, 6, 6, 6, 6, 6, 6]
    }
    eta_range = {
        0: [0, 0.2, 0.4],
        1: [0, 0.2],
        2: [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0],
        3: [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0]
    }

    alpha_dict, eta_dict = {}, {}
    ms = {0: "o", 1: "s", 2: "D", 3: ">"}
    color = {0: "tab:blue", 1: "tab:red", 2: "tab:green", 3: "tab:orange"}
    label = {0: "H", 1: "DP", 2: "SP", 3: "DP+SP"}

    for key in sorted(alpha_l.keys()):
        for j, eta_j in enumerate(eta_range[key]):
            for alpha_i in alpha:
                if alpha_l[key][j] <= alpha_i <= alpha_h[key][j]:
                    if key in alpha_dict:
                        alpha_dict[key].append(alpha_i)
                        eta_dict[key].append(eta_j)
                    else:
                        alpha_dict[key] = [alpha_i]
                        eta_dict[key] = [eta_j]
        plt.plot(alpha_dict[key],
                 eta_dict[key],
                 ms[key],
                 c=color[key],
                 label=label[key])
    x = np.linspace(0, 4.05, 100)
    y = x
    plt.plot(x, y, "r-", label="EP line")
    plt.xlabel(r"$\tilde{\alpha}$", fontsize="x-large")
    plt.ylabel(r"$-\tilde{\eta}$", fontsize="x-large")
    plt.title(
        r"$L_x=20,L_y=5,\phi_A=\phi_B=\rho_0=100,v_0=1,D_r=0.02,\kappa=1$")
    plt.legend(loc="upper left", fontsize="large")
    plt.xlim(0, 6.05)
    plt.ylim(0, 4.05)
    plt.tight_layout()
    plt.show()
    plt.close()


def Dr_vs_alpha(Lx=20, Ly=5, phiA=80, phiB=80, rho0=80, k=0.7, eta=0):
    fout = f"fig/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_e{eta:g}.png"
    fnpz = f"data/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_e{eta:g}.npz"

    if os.path.exists(fnpz):
        with np.load(fnpz, "r") as data:
            alpha = data["alpha"]
            Dr = data["Dr"]
            v_std = data["v_std"]
        print("load data from", fnpz)
    else:
        print("cal v_std from raw data...")
        if Lx == 20 and Ly == 5 and phiA == phiB == rho0 == 80:
            seed = 1001
            # sftp://yduan@sohrab001
            folder = "/scratch03.local/yduan/QS5/L20_5_k0.7_alpha"
            pat = f"{folder}/L{Lx}_{Ly}_*_k{k:.2f}_" \
                f"p{phiA:g}_{phiB:g}_r{rho0:g}_" \
                f"e{eta:.3f}_*_{seed}.gsd"
        files = glob.glob(pat)
        alpha, Dr = np.zeros((2, len(files)))
        v_std = np.zeros_like(alpha)

        for i, f in enumerate(files):
            paras = get_paras(["Dr", "a"], f)
            alpha[i] = paras['a']
            Dr[i] = paras['Dr']
            snap = read_one_frame(f, -1)
            v_std[i] = np.std(snap.particles.charge)
        np.savez_compressed(fnpz, alpha=alpha, Dr=Dr, v_std=v_std)
    
    # plt.plot(alpha, Dr, "o")
    fig, ax = plt.subplots(constrained_layout=True)
    sca = ax.scatter(alpha, Dr, c=v_std)
    ax.set_xscale("log")
    ax.set_yscale("log")
    cb = fig.colorbar(sca, ax=ax)
    cb.set_label(r"$\Delta v$", fontsize="large")
    ax.set_xlabel(r"$\alpha$", fontsize="large")
    ax.set_ylabel(r"$D_r$", fontsize="large")
    title = f"$L_x={Lx},L_y={Ly}," \
        f"\\phi_A=\\phi_B=\\rho_0={phiA:g}," \
        f"\\kappa={k:g},\\eta={eta:g}$"
    fig.suptitle(title)
    # plt.show()
    plt.savefig(fout)
    plt.close()


def get_v_std_mean(fgsd, beg=-100):
    """ Cal time-averaged v_std

    Args:
        fgsd (str): filename of the gsd file
        beg (int, optional): Begining frame for time average.
            Defaults to -100.
    """
    v_std_mean = 0
    count = 0
    frames = read_frames(fgsd, beg=beg)
    for frame in frames:
        v_std_mean += np.std(frame.particles.charge)
        count += 1
    v_std_mean /= count
    return v_std_mean


def eta_vs_alpha(Lx=20, Ly=5, phiA=80, phiB=80, rho0=80, k=0.7, Dr=3):
    fout = f"fig/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_Dr{Dr:g}.png"
    fnpz = f"data/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_Dr{Dr:g}.npz"
    fdat = f"data/PD/L{Lx}_{Ly}_p{phiA}_{phiB}_r{rho0}_k{k}_Dr{Dr:g}.dat"

    paras = read_txt(fdat, ["eta", "alpha", "state"])
    eta = paras["eta"]
    alpha = paras["alpha"]
    state = paras["state"]

    if os.path.exists(fnpz):
        with np.load(fnpz, "r") as data:
            # alpha = data["alpha"]
            # eta = data["eta"]
            v_std = data["v_std"]
        print("load data from", fnpz)
    else:
        print("cal v_std from raw data...")
        if Lx == 20 and Ly == 5 and phiA == phiB == rho0 == 80:
            seed = 1001
            # sftp://yduan@sohrab001
            folder = "/scratch03.local/yduan/QS5/L20_5_k0.7_alpha"
            pat = f"{folder}/L{Lx}_{Ly}_Dr{Dr:.3f}_k{k:.2f}_" \
                f"p{phiA:g}_{phiB:g}_r{rho0:g}_" \
                f"e%.3f_a%.3f_{seed}.gsd"
        v_std = np.zeros_like(alpha)
        for i in range(alpha.size):
            fgsd = pat % (eta[i], alpha[i])
            v_std[i] = get_v_std_mean(fgsd)
        np.savez_compressed(fnpz, alpha=alpha, eta=eta, v_std=v_std)

    fig, ax = plt.subplots(constrained_layout=True)

    # phase boundaries predicted by theory
    x = np.linspace(0, 4.05, 100)
    ax.plot(x, x, "r-", label="EP line")
    x = np.linspace(0, 1, 100)
    y = (1 + x**2) / 2
    ax.plot(x, y, "b-", label=r"$\lambda_+=0$")
    x = np.linspace(1, 4.)
    y = (1 + x**2) / 2
    ax.plot(x, y, "-", label=r"$\lambda_-=0$", c="tab:orange")
    ax.axhline(1, linestyle="--", c="k")

    mk = {0: "o", 1: "s", 2: "D", 3: ">"}
    label = {0: "H", 1: "DP", 2: "SP", 3: "DP+SP"}
    for state_id in mk.keys():
        mask = state == state_id
        sca = ax.scatter(alpha[mask], -eta[mask], c=v_std[mask], vmin=0,
            vmax=v_std.max(), marker=mk[state_id], label=label[state_id])

    ax.set_ylim(-0.1, 3.2)
    cb = fig.colorbar(sca, ax=ax)
    cb.set_label(r"$\Delta v$", fontsize="large")
    ax.set_xlabel(r"$\alpha$", fontsize="large")
    ax.set_ylabel(r"$\eta$", fontsize="large")
    ax.legend()
    title = f"$L_x={Lx},L_y={Ly}," \
        f"\\phi_A=\\phi_B=\\rho_0={phiA:g}," \
        f"\\kappa={k:g},D_r={Dr:g}$"
    fig.suptitle(title)
    # plt.show()
    plt.savefig(fout)
    plt.close()



if __name__ == "__main__":
    Dr_vs_alpha()

    # eta_vs_alpha()
