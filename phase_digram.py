from re import X
import numpy as np
import matplotlib.pyplot as plt


def Dr_alpha_eta0():
    alpha = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                     1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4])
    Dr = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.10,
                  0.2, 0.4, 0.6, 0.8, 1.0, 2, 3])
    alpha_c = np.array([0, 0.3, 0.3, 0.4, 0.4, 0.4, 0.6,
                       0.9, 1.0, 1.2, 1.4, 3, 4, 4])

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
        plt.plot(alpha_dict[key], Dr_dict[key], ms[key],
                 c=color[key], label=label[key])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\tilde{\alpha}$", fontsize="x-large")
    plt.ylabel(r"$D_r$", fontsize="x-large")
    plt.title(
        r"$L_x=20,L_y=5,\phi_A=\phi_B=\rho_0=100,v_0=1,\tilde{\eta}=0,\kappa=1$")
    plt.legend(loc="lower right", fontsize="x-large")
    plt.tight_layout()
    plt.show()
    plt.close()


if __name__ == "__main__":
    alpha = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4,
                     1.6, 1.8, 2, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0])
    eta = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2,
                   1.4, 1.6, 2.0, 3.0, 4.0])
    alpha_l = {0: [0, 0, 0], 1: [0.4, 0.5],
               2: [0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               3: [1.4, 2.4, 2.8, 3.2, 3.2, 3.6, 4.0, 4.4, 5.6]}
    alpha_h = {0: [0.3, 0.4, 0.5], 1: [6, 6],
               2: [1.2, 2, 2.4, 2.8, 2.8, 3.2, 3.6, 4.0, 5.2, 6], 
               3: [6, 6, 6, 6, 6, 6, 6, 6, 6]}
    eta_range = {0: [0, 0.2, 0.4], 1: [0, 0.2],
                 2: [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0, 4.0],
                 3: [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 3.0]}

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
        plt.plot(alpha_dict[key], eta_dict[key], ms[key],
                 c=color[key], label=label[key])
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
