import numpy as np
import matplotlib.pyplot as plt


def get_Dr_threshold(eta, eta_AB, eta_BA, v0=1, q_thresh=np.pi*2):
    mu = 2 + 2 * eta
    nu = (-eta_AB * eta_BA) ** 0.5
    a = nu**2 / mu - 1/8
    Dr_thresh = v0 * np.sqrt(2) / 16 * (a - np.sqrt(a**2 - 1/64)) ** (-0.5) * q_thresh
    return Dr_thresh

def cal_alpha_threshold(eta, Dr, qc=np.pi*2, v0=1):
    a = 2 * (v0/(16 * Dr))** 2 * qc ** 2
    c = (a + 1/8) ** 2 / (2 * a)
    return np.sqrt(eta ** 2 + 2 * c * (1+eta))

def cal_etaAB_threshold(eta, Dr, qc=np.pi*2, v0=1):
    a = 2 * (v0/(16 * Dr))** 2 * qc ** 2
    c = (a + 1/8) ** 2 / (2 * a)
    return np.sqrt(2 * c * (1+eta))

def cal_lambda_pp(eta, eta_AB, eta_BA, Dr, q, v0=1):
    qq = q ** 2
    D = Dr / v0 + v0 / (16 * Dr) * q ** 2
    mu = 2 + eta + eta
    nu = np.sqrt(-eta_AB * eta_BA)
    return 0.5 * (-D + np.sqrt(D*D - mu * qq + 2j * nu * qq))


def plot_lambda_varied_eta_AB():
    q = np.logspace(-2, 3, 1000)

    eta = 0
    Dr = 0.1
    eta_AB = 0.8
    eta_BA = - eta_AB
    Dt = 0.01

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), constrained_layout=True)
    

    for eta_AB in [0.75, 1, 1.5, 2, 3]:
        eta_BA = -eta_AB
        y0 = cal_lambda_pp(eta, eta_AB, eta_BA, Dr, q)
        y = y0.real - Dt * q**2  
        line, = ax1.plot(q, y, label=r"$%g$" % eta_AB)
        ax2.axvline(q[np.argmax(y)], linestyle=":", color=line.get_c())
        print(q[np.argmax(y)])
        ax2.plot(q, y0.imag, label=r"$%g$" % eta_AB)

    ax1.axhline(0, linestyle="dashed", c="k")

    # ax.axvline(np.pi * 2, linestyle="dotted", c="tab:red")
    ax1.set_xscale("log")
    ax2.set_xscale("log")

    ax1.set_xlim(q[0], 10)
    ax2.set_xlim(q[0], 10)

    ax1.set_ylim(-1, 0.5)
    ax1.set_xlabel(r"$q$", fontsize="x-large")
    ax2.set_xlabel(r"$q$", fontsize="x-large")

    ax1.set_ylabel(r"$\Re(\lambda_{++})$", fontsize="x-large")
    ax2.set_ylabel(r"$\Im(\lambda_{++})$", fontsize="x-large")
    ax1.legend(title=r"$\eta_{AB}=$", title_fontsize="xx-large")
    fig.suptitle("$D_r=%g,\eta_{AA}=\eta_{BB}=%g, D_t=%g, \eta_{AB}=-\eta_{BA}$" % (Dr, eta, Dt), fontsize="large")
    plt.show()
    # plt.savefig("Re_lambda_Dt.pdf")
    plt.close()


if __name__ == "__main__":
    plot_lambda_varied_eta_AB()
