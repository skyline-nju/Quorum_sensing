import numpy as np
import matplotlib.pyplot as plt

def cal_mu(eta_AA, eta_BB):
    return 1 + 0.5 * (eta_AA + eta_BB)

def cal_Delta(eta_AB, eta_BA, eta_AA, eta_BB):
    return 0.25 * (eta_AA - eta_BB)**2 + eta_AB * eta_BA

def cal_minus_Delta_over_mu(Dr, v0=1, qc=np.pi*2):
    x = v0 * qc / Dr
    if x >= 4:
        return 0.5
    else:
        return 2 * (x/16 + 1/x) ** 2
    
def cal_Dr_over_v0(mu, Delta, qc=np.pi*2):
    x = -Delta / mu
    if isinstance(x, np.ndarray):
        y = np.zeros_like(x)
        mask = x <= 0.5
        y[mask] = 0
        mask = x > 0.5
        y[mask] = 0.5 * qc * (np.sqrt(x[mask]/2) + np.sqrt(x[mask]/2 - 0.25))
        return y
    else:
        if x <= 0.5:
            return 0
        else:
            return 0.5 * qc * (np.sqrt(x/2) + np.sqrt(x/2 - 0.25))

def cal_Dr_over_v0_2(mu, Delta, qc=np.pi*2):
    x = -Delta / mu
    if isinstance(x, np.ndarray):
        y = np.zeros_like(x)
        mask = x <= 0.5
        y[mask] = 0
        mask = x > 0.5
        y[mask] = 0.5 * qc * (np.sqrt(x[mask]/2) - np.sqrt(x[mask]/2 - 0.25))
        return y
    else:
        if x <= 0.5:
            return 0
        else:
            return 0.5 * qc * (np.sqrt(x/2) - np.sqrt(x/2 - 0.25))


def PD_Delta_mu(Dr, v0=1, qc=np.pi*2, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False

    # Delta
    x = np.linspace(-5, 5)

    # mu
    y = np.linspace(-5, 5)

    ax.axvline(0, c="tab:red", label=r"$\Delta=0$")
    ax.axhline(0, c="k", label=r"$\mu=0$")

    # Delta = mu **2
    x1 = y ** 2
    ax.plot(x1, y, c="tab:blue", label=r"$\Delta=\mu^2$")

    # minus Delta over mu
    k = cal_minus_Delta_over_mu(Dr, v0, qc)
    y2 = np.linspace(0, y.max(), 100)  # mu > 0
    x2 = -k * y2
    # label = r"$\frac{-\Delta}{\mu}=\left(\frac{-\Delta}{\mu}\right)^*$"
    label = r"$-\Delta/\mu = (-\Delta/\mu)^*$"
    ax.plot(x2, y2, c="tab:green", label=label)  # q*=q_-
    # plt.plot(-0.5 * y2, y2)

    ax.set_xlim(-4, 4)
    ax.set_ylim(-5, 5)
    ax.set_xlabel(r"$\Delta$", fontsize="x-large")
    ax.set_ylabel(r"$\mu$", fontsize="x-large")
    ax.text(-2, -1, "DP " + u"\u2160", fontsize="x-large")
    ax.text(-3, 2, "DP " + u"\u2161", fontsize="x-large")
    ax.text(-1, 4, "H", fontsize="x-large")
    ax.text(2, 4, "H", fontsize="x-large")
    ax.text(2, 0.5, "SP", fontsize="x-large")
    ax.text(2, -1, "SP", fontsize="x-large")
    ax.text(2, -3, "SP", fontsize="x-large")
    ax.legend()

    if flag_show:
        plt.show()
        plt.close()

def PD_alpha_eta(Dr, v0=1, qc=np.pi*2, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False

    alpha = np.linspace(0, 5)
    eta = np.linspace(-5, 5)
    
    # Exceptional Lines
    ax.plot(alpha, alpha, c="tab:red", label=r"$|\alpha/\eta|=1$")
    ax.plot(alpha, -alpha, c="tab:red")

    ax.axhline(-1, c="k", label=r"$\eta=-1$")

    y = -(1 + alpha**2)/2
    ax.plot(alpha, y, c="tab:blue", label=r"$\eta=-(1+\alpha^2)/2$")
    # ax.fill_between(alpha, -3.0, y, alpha=0.5)


    # minus Delta over mu
    k = cal_minus_Delta_over_mu(Dr, v0, qc)
    y = np.linspace(-1, eta.max())  # eta < -1
    x = np.sqrt((1+y)* k + y**2)
    label = r"$-\Delta/\mu = (-\Delta/\mu)^*$"
    ax.plot(x, y, label=label, c="tab:green")

    ax.set_xlim(0, 3)
    ax.set_ylim(-3, 1.5)
    ax.set_xlabel(r"$|\alpha|$", fontsize="x-large")
    ax.set_ylabel(r"$\eta$", fontsize="x-large")

    ax.text(2.5, -1.5, "DP " + u"\u2160", fontsize="x-large")
    ax.text(2.5, 0.5, "DP " + u"\u2161", fontsize="x-large")
    # ax.text(2.5, -1.5, "DP", fontsize="x-large")
    ax.text(0.3, 1.0, "H", fontsize="x-large")
    ax.text(0.5, 0, "H", fontsize="x-large")
    ax.text(0.05, -0.45, "H", fontsize="x-large")
    ax.text(2.3, -2.9, "SP", fontsize="x-large")
    ax.text(1.0, -2.5, "SP", fontsize="x-large")
    ax.text(0.1, -0.9, "SP", fontsize="x-large")

    # ax.legend(loc="upper right")    
    if flag_show:
        plt.show()
        # plt.savefig("PD_rho_eq.pdf")
        plt.close()


def PD_alpha_Dr(eta=0, v0=1, qc=np.pi*2, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    
    alpha_c = np.sqrt(0.5 * (1+eta) + eta**2)
    alpha = np.linspace(alpha_c, 10, 1000)
    mu = cal_mu(eta, eta)
    Delta = cal_Delta(eta + alpha, eta - alpha, eta, eta)
    Dr = cal_Dr_over_v0(mu, Delta, qc=qc) * v0
    line, = ax.plot(alpha, Dr)
    ax.axvline(alpha_c, linestyle="dashed", c="tab:pink", lw=2)
    # ax.set_xlim(0, 3.5)
    # ax.set_ylim(0, 15)
    # ax.set_xlabel(r"$|\alpha|$", fontsize="x-large")
    # ax.set_ylabel(r"$D_r/v_0$", fontsize="x-large")

    # ax.text(2.5, 8, "DP " + u"\u2161", fontsize="x-large")
    # ax.text(0.5, 8, "H", fontsize="x-large")
    # ax.text(1.2, 8, "H", fontsize="x-large")

    # ax.legend()
    if flag_show:
        plt.show()
        plt.close()
    else:
        return line


def PD_etaAB_eta(Dr, v0=1, qc=np.pi*2, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4),  constrained_layout=True)
        flag_show = True
    else:
        flag_show = False

    etaAB = np.linspace(0, 5)
    eta = np.linspace(-5, 0)
    
    ax.axvline(0, c="tab:red", label=r"$\eta_{AB}=0$", lw=4)
    # ax.plot(alpha, alpha, c="tab:red", label=r"$|\alpha/\eta|=1$")
    ax.axhline(-1, c="k", label=r"$\eta=-1$")

    # minus Delta over mu
    k = cal_minus_Delta_over_mu(Dr, v0, qc)
    y = np.linspace(-1, 5)  # eta < -1
    x = np.sqrt((1+y)* k)
    label = r"$-\Delta/\mu = (-\Delta/\mu)^*$"
    ax.plot(x, y, label=label, c="tab:green")

    ax.set_xlim(0, 4)
    ax.set_ylim(-3, 3)
    ax.set_xlabel(r"$|\eta_{AB}|$", fontsize="x-large")
    ax.set_ylabel(r"$\eta$", fontsize="x-large")

    ax.text(1.5, -2, "DP " + u"\u2160", fontsize="x-large")
    ax.text(1.5, 0.5, "DP " + u"\u2161", fontsize="x-large")
    ax.text(0.2, 0.2, "H", fontsize="x-large")

    # ax.legend()    
    if flag_show:
        plt.show()
        plt.close()