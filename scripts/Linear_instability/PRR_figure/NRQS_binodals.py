import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize


def get_v_XY(rho_Y, eta_XY, kappa=0.7):
    return 1 + kappa * np.tanh(eta_XY/kappa * (rho_Y-1))

def get_v_A(rho_A, rho_B, eta_AA, eta_AB, kappa=0.7):
    v_AA = get_v_XY(rho_A, eta_AA, kappa)
    v_AB = get_v_XY(rho_B, eta_AB, kappa)
    return v_AA * v_AB

def get_v_B(rho_A, rho_B, eta_BA, eta_BB, kappa=0.7):
    v_BA = get_v_XY(rho_A, eta_BA, kappa)
    v_BB = get_v_XY(rho_B, eta_BB, kappa)
    return v_BA * v_BB

def get_dv_XY(v_XY, eta_XY, kappa=0.7):
    return eta_XY * (1 - ((v_XY-1)/kappa)**2)

def get_v_XY_prime(rho_Y, eta_XY, kappa=0.7):
    v_XY = get_v_XY(rho_Y, eta_XY, kappa)
    return eta_XY * (1 - ((v_XY-1)/kappa)**2)

def get_mu_A(rho_A, rho_B, eta_AA, eta_AB, kappa=0.7):
    v_A = get_v_A(rho_A, rho_B, eta_AA, eta_AB, kappa)
    return np.log(rho_A * v_A)

def get_mu_B(rho_A, rho_B, eta_BA, eta_BB, kappa=0.7):
    v_B = get_v_B(rho_A, rho_B, eta_BA, eta_BB, kappa)
    return np.log(rho_B * v_B)

def f_int_A(phi_A, eta_AA, eta_BA):
    v_AA = get_v_XY(phi_A, eta_AA)
    v_BA = get_v_XY(phi_A, eta_BA)
    v_AA_prime = get_dv_XY(v_AA, eta_AA)
    return (1/phi_A + v_AA_prime/v_AA) * np.log(v_BA)


def f_int_B(phi_B, eta_AB, eta_BB):
    v_AB = get_v_XY(phi_B, eta_AB)
    v_BB = get_v_XY(phi_B, eta_BB)
    v_BB_prime = get_dv_XY(v_BB, eta_BB)
    return (1/phi_B + v_BB_prime/v_BB) * np.log(v_AB)

def get_RA_RB(phi_A, phi_B, eta_BA, eta_AB):
    v_BA = get_v_XY(phi_A, eta_BA)
    v_AB = get_v_XY(phi_B, eta_AB)
    return np.log(v_BA) * np.log(v_AB)


def func(phi_arr, eta_AA, eta_AB, eta_BA, eta_BB):
    """ When phase equilibrium is reached, one should have [f1, f2, f3] = [0, 0, 0]

    Args:
        phi_arr (array): [phi_A_g, phi_B_g, phi_A_l, phi_B_l]
        eta_AA (float): _description_
        eta_AB (float): _description_
        eta_BA (float): _description_
        eta_BB (float): _description_
    """
    phi_A_g, phi_B_g, phi_A_l, phi_B_l = phi_arr
    f1 = get_mu_A(phi_A_l, phi_B_l, eta_AA, eta_AB) - get_mu_A(phi_A_g, phi_B_g, eta_AA, eta_AB)
    f2 = get_mu_B(phi_A_l, phi_B_l, eta_BA, eta_BB) - get_mu_B(phi_A_g, phi_B_g, eta_BA, eta_BB)

    int_A = integrate.quad(f_int_A, phi_A_g, phi_A_l, args=(eta_AA, eta_BA))
    int_B = integrate.quad(f_int_B, phi_B_g, phi_B_l, args=(eta_AB, eta_BB))
    f3 = get_RA_RB(phi_A_l, phi_B_l, eta_BA, eta_AB) - get_RA_RB(phi_A_g, phi_B_g, eta_BA, eta_AB) + int_A + int_B
    return [f1, f2, f3]

def func2(x_arr, eta_AA, eta_AB, eta_BA, eta_BB, rhoA_0, rhoB_0):
    """ When phase equilibrium is reached, one should have [f1, f2, f3, f4, f5] = [0, 0, 0, 0, 0]

    Args:
        x_arr (array): [phi_A_g, phi_B_g, phi_A_l, phi_B_l, s]
        eta_AA (float): _description_
        eta_AB (float): _description_
        eta_BA (float): _description_
        eta_BB (float): _description_
    """
    phi_A_g, phi_B_g, phi_A_l, phi_B_l, s = x_arr
    f1 = get_mu_A(phi_A_l, phi_B_l, eta_AA, eta_AB) - get_mu_A(phi_A_g, phi_B_g, eta_AA, eta_AB)
    f2 = get_mu_B(phi_A_l, phi_B_l, eta_BA, eta_BB) - get_mu_B(phi_A_g, phi_B_g, eta_BA, eta_BB)

    int_A, err = integrate.quad(f_int_A, phi_A_g, phi_A_l, args=(eta_AA, eta_BA))
    int_B, err = integrate.quad(f_int_B, phi_B_g, phi_B_l, args=(eta_AB, eta_BB))
    f3 = get_RA_RB(phi_A_l, phi_B_l, eta_BA, eta_AB) - get_RA_RB(phi_A_g, phi_B_g, eta_BA, eta_AB) + int_A + int_B

    f4 = (1-s) * phi_A_g + s * phi_A_l - rhoA_0
    f5 = (1-s) * phi_B_g + s * phi_B_l - rhoB_0
    return [f1, f2, f3, f4, f5]

def func3(x_arr, eta_AA, eta_AB, eta_BA, eta_BB, rhoA_0, rhoB_0):
    """ When phase equilibrium is reached, one should have [f1, f2, f3, f4, f5] = [0, 0, 0, 0, 0]

    Args:
        x_arr (array): [phi_A_1, phi_B_1, phi_A_2, phi_B_2, phi_A_3, phi_B_3, s, t]
        eta_AA (float): _description_
        eta_AB (float): _description_
        eta_BA (float): _description_
        eta_BB (float): _description_
    """
    phi_A_1, phi_B_1, phi_A_2, phi_B_2, phi_A_3, phi_B_3, s, t = x_arr
    
    mu_A_1 = get_mu_A(phi_A_1, phi_B_1, eta_AA, eta_AB)
    mu_A_2 = get_mu_A(phi_A_2, phi_B_2, eta_AA, eta_AB)
    mu_A_3 = get_mu_A(phi_A_3, phi_B_3, eta_AA, eta_AB)
    mu_B_1 = get_mu_B(phi_A_1, phi_B_1, eta_BA, eta_BB)
    mu_B_2 = get_mu_B(phi_A_2, phi_B_2, eta_BA, eta_BB)
    mu_B_3 = get_mu_B(phi_A_3, phi_B_3, eta_BA, eta_BB)

    f1 = mu_A_1 - mu_A_2
    f2 = mu_A_1 - mu_A_3
    f3 = mu_B_1 - mu_B_2
    f4 = mu_B_1 - mu_B_3

    RA_RB_1 = get_RA_RB(phi_A_1, phi_B_1, eta_BA, eta_AB)
    int_A, err = integrate.quad(f_int_A, phi_A_1, phi_A_2, args=(eta_AA, eta_BA))
    int_B, err = integrate.quad(f_int_B, phi_B_1, phi_B_2, args=(eta_AB, eta_BB))
    f5 = get_RA_RB(phi_A_2, phi_B_2, eta_BA, eta_AB) - RA_RB_1 + int_A + int_B

    int_A, err = integrate.quad(f_int_A, phi_A_1, phi_A_3, args=(eta_AA, eta_BA))
    int_B, err = integrate.quad(f_int_B, phi_B_1, phi_B_3, args=(eta_AB, eta_BB))
    f6 = get_RA_RB(phi_A_3, phi_B_3, eta_BA, eta_AB) - RA_RB_1 + int_A + int_B

    f7 = (1-s-t) * phi_A_1 + s * phi_A_2 + t * phi_A_3 - rhoA_0
    f8 = (1-s-t) * phi_B_1 + s * phi_B_2 + t * phi_B_3 - rhoB_0
    return [f1, f2, f3, f4, f5, f6, f7, f8]

def func2_eq(x_arr, eta_AA, eta_AB, eta_BA, eta_BB, rhoA_0, rhoB_0):
    """ When phase equilibrium is reached, one should have [f1, f2, f3, f4, f5] = [0, 0, 0, 0, 0]

    Args:
        phi_arr (array): [phi_A_g, phi_B_g, phi_A_l, phi_B_l, s]
        eta_AA (float): _description_
        eta_AB (float): _description_
        eta_BA (float): _description_
        eta_BB (float): _description_
    """

    phi_A_g, phi_B_g, phi_A_l, phi_B_l, s = x_arr
    f1 = get_mu_A(phi_A_l, phi_B_l, eta_AA, eta_AB) - get_mu_A(phi_A_g, phi_B_g, eta_AA, eta_AB)
    f2 = get_mu_B(phi_A_l, phi_B_l, eta_BA, eta_BB) - get_mu_B(phi_A_g, phi_B_g, eta_BA, eta_BB)

    # int_A, err = integrate.quad(f_int_A, phi_A_g, phi_A_l, args=(eta_AA, eta_BA))
    # int_B, err = integrate.quad(f_int_B, phi_B_g, phi_B_l, args=(eta_AB, eta_BB))
    int_A, err = integrate.quad(get_v_XY, phi_A_g, phi_A_l, args=(eta_AA, 0.7))
    int_B, err = integrate.quad(get_v_XY, phi_B_g, phi_B_l, args=(eta_BB, 0.7))
    vAA_g = get_v_XY(phi_A_g, eta_AA)
    vAA_l = get_v_XY(phi_A_l, eta_AA)
    vBB_g = get_v_XY(phi_B_g, eta_BB)
    vBB_l = get_v_XY(phi_B_l, eta_BB)
    vAB_prime_g = get_v_XY_prime(phi_B_g, eta_AB)
    vAB_prime_l = get_v_XY_prime(phi_B_l, eta_AB)
    vBA_prime_g = get_v_XY_prime(phi_A_g, eta_BA)
    vBA_prime_l = get_v_XY_prime(phi_A_l, eta_BA)

    f3_0 = phi_A_l * (1 + vAA_l) - phi_A_g * (1 + vAA_g) + phi_B_l * (1 + vBB_l) - phi_B_g * (1 + vBB_g)
    f3_1 = 0.5 * phi_A_l * phi_B_l * (vBA_prime_l + vAB_prime_l) - 0.5 * phi_A_g * phi_B_g * (vBA_prime_g + vAB_prime_g)
    f3 = f3_0 + f3_1 - int_A - int_B
    
    f4 = (1-s) * phi_A_g + s * phi_A_l - rhoA_0
    f5 = (1-s) * phi_B_g + s * phi_B_l - rhoB_0
    return [f1, f2, f3, f4, f5]

def func3_eq(x_arr, eta_AA, eta_AB, eta_BA, eta_BB, rhoA_0, rhoB_0):
    """ When phase equilibrium is reached, one should have [f1, f2, f3, f4, f5] = [0, 0, 0, 0, 0]

    Args:
        x_arr (array): [phi_A_1, phi_B_1, phi_A_2, phi_B_2, phi_A_3, phi_B_3, s, t]
        eta_AA (float): _description_
        eta_AB (float): _description_
        eta_BA (float): _description_
        eta_BB (float): _description_
    """
    def get_p1(phiA, phiB, etaAA, etaAB, etaBA, etaBB):
        p1A = phiA * (1 + get_v_XY(phiA, etaAA))
        p1B = phiB * (1 + get_v_XY(phiB, etaBB))
        p1AB = 0.5 * phiA * phiB * (get_v_XY_prime(phiA, etaBA) + get_v_XY_prime(phiB, etaAB))
        return p1A + p1B + p1AB
    
    def get_p_int(phi_A_1, phi_B_1, phi_A_2, phi_B_2, etaAA, etaBB):
        int_A, err = integrate.quad(get_v_XY, phi_A_1, phi_A_2, args=(etaAA, 0.7))
        int_B, err = integrate.quad(get_v_XY, phi_B_1, phi_B_2, args=(etaBB, 0.7))
        return int_A + int_B



    phi_A_1, phi_B_1, phi_A_2, phi_B_2, phi_A_3, phi_B_3, s, t = x_arr
    
    mu_A_1 = get_mu_A(phi_A_1, phi_B_1, eta_AA, eta_AB)
    mu_A_2 = get_mu_A(phi_A_2, phi_B_2, eta_AA, eta_AB)
    mu_A_3 = get_mu_A(phi_A_3, phi_B_3, eta_AA, eta_AB)
    mu_B_1 = get_mu_B(phi_A_1, phi_B_1, eta_BA, eta_BB)
    mu_B_2 = get_mu_B(phi_A_2, phi_B_2, eta_BA, eta_BB)
    mu_B_3 = get_mu_B(phi_A_3, phi_B_3, eta_BA, eta_BB)

    f1 = mu_A_1 - mu_A_2
    f2 = mu_A_1 - mu_A_3
    f3 = mu_B_1 - mu_B_2
    f4 = mu_B_1 - mu_B_3


    f5 = get_p1(phi_A_2, phi_B_2, eta_AA, eta_AB, eta_BA, eta_BB) \
        - get_p1(phi_A_1, phi_B_1, eta_AA, eta_AB, eta_BA, eta_BB) \
        - get_p_int(phi_A_1, phi_B_1, phi_A_2, phi_B_2, eta_AA, eta_BB)

    f6 = get_p1(phi_A_3, phi_B_3, eta_AA, eta_AB, eta_BA, eta_BB) \
        - get_p1(phi_A_1, phi_B_1, eta_AA, eta_AB, eta_BA, eta_BB) \
        - get_p_int(phi_A_1, phi_B_1, phi_A_3, phi_B_3, eta_AA, eta_BB)

    f7 = (1-s-t) * phi_A_1 + s * phi_A_2 + t * phi_A_3 - rhoA_0
    f8 = (1-s-t) * phi_B_1 + s * phi_B_2 + t * phi_B_3 - rhoB_0
    return [f1, f2, f3, f4, f5, f6, f7, f8]

def plot_tri(ax, x, y, c, linestyle="-", lw=None):
    ax.plot(x[1:], y[1:], c=c, linestyle=linestyle, lw=lw)
    ax.plot(x[:2], y[:2], c=c, linestyle=linestyle, lw=lw)
    ax.plot([x[0], x[2]], [y[0], y[2]], c=c, linestyle=linestyle, lw=lw)


def binodals_J025():
    etaAA = -2
    etaBB = -2
    etaAB = -0.25
    etaBA = -0.25
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))

    # G + LA + LAB
    rhoA_0 = 2.25
    rhoB_0 = 1
    x0 = [0.5, 0.25, 2.5, 0.4, 5, 3, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    line, = ax.plot(pA1, pB1, ">")
    plot_tri(ax, pA1, pB1, c=line.get_c())

    rhoA_0 = 2.25
    rhoB_0 = 1
    x0 = [0.5, 0.25, 2.5, 0.4, 5, 3, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, -0.1, -0.1, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA_tmp = [x[0], x[2], x[4]]
    pB_tmp = [x[1], x[3], x[5]]
    line, = ax.plot(pA_tmp, pB_tmp, ">")
    plot_tri(ax, pA_tmp, pB_tmp, c=line.get_c())


    # G + LB + LAB
    rhoA_0 = 1
    rhoB_0 = 2.25
    x0 = [0.25, 0.5, 0.4, 2.5, 3, 5, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA2 = [x[0], x[2], x[4]]
    pB2 = [x[1], x[3], x[5]]
    line, = ax.plot(pA2, pB2, "<")
    plot_tri(ax, pA2, pB2, c=line.get_c())

    # G + LA
    rhoA_0 = 1.25
    rhoB_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0]) 
    rhoB_0 = 0.01
    n = 40
    rhoB0_arr = np.linspace(0.01, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    line, = ax.plot(phi_A_g, phi_B_g)
    ax.plot(phi_A_l, phi_B_l, c=line.get_c())

    edges = np.array([
        [3.8241, 0.3792, 23.2144, 0.6518],
        [3.8257, 0.7560, 23.1665, 1.3016],
        [3.7998, 1.1292, 23.2603, 1.9610],
        [3.7808, 1.5073, 23.3153, 2.6195]
    ]) / 10

    ms = 3
    line, = ax.plot(edges[:, 0], edges[:, 1], "o", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)



    # G + LB
    rhoB_0 = 1.25
    rhoA_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0])
    rhoA_0_min = 0.01
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    line, = ax.plot(phi_A_g, phi_B_g)
    ax.plot(phi_A_l, phi_B_l, c=line.get_c())


    # G + LAB
    rhoA_0 = 3
    rhoB_0_min = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[2] - pB1[0]) / (pA1[2] - pA1[0]) 
    rhoB_0_max = pB2[0] + (rhoA_0 - pA2[0]) * (pB2[2] - pB2[0]) / (pA2[2] - pA2[0])
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[0], pB1[0], pA1[2], pB1[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    line, = ax.plot(phi_A_g, phi_B_g)
    ax.plot(phi_A_l, phi_B_l, c=line.get_c())
    # ax.plot(np.ones_like(rhoB0_arr) * rhoA_0, rhoB0_arr)

    # LB + LAB
    rhoA_0 = 1.5
    rhoB_0_min = pB2[1] + (rhoA_0 - pA2[1]) * (pB2[2] - pB2[1]) / (pA2[2] - pA2[1])
    rhoB_0_max = 8.5
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[1], pB2[1], pA2[2], pB2[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    line, = ax.plot(phi_A_g, phi_B_g)
    ax.plot(phi_A_l, phi_B_l, c=line.get_c())

    # LA + LAB
    rhoB_0 = 1.5
    rhoA_0_min = pA1[1] + (rhoA_0 - pB1[1]) * (pA1[2] - pA1[1]) / (pB1[2] - pB1[1])
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    line, = ax.plot(phi_A_g, phi_B_g)
    ax.plot(phi_A_l, phi_B_l, c=line.get_c())

 


    # ax.plot([0.376, 2, 4.293], [0.211, 1.5, 3.302], "o")
    # ax.plot([0.299, 2, 4.73], [0.186, 1.5, 3.63], "s")

    rhoA_0 = 2
    rhoB_0 = 1.5
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.376, 0.211, 4.293, 3.302, 0.5]
    # x0 = [0.376, 0.211, 3, 2, 0.5]

    # sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    sol = optimize.root(func2_eq, x0, args=(-2, -0.1, -0.1, -2, rhoA_0, rhoB_0))
    x = sol.x
    ax.plot([x[0], x[2]], [x[1], x[3]], "-*")

    # sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    sol = optimize.root(func2, x0, args=(-2, -0.1, -0.1, -2, rhoA_0, rhoB_0))
    x = sol.x
    ax.plot([x[0], x[2]], [x[1], x[3]], "-p")


    x0 = [3.7808/10, 1.5073/10, 23.3153/10, 2.6195/10, 0.5]
    sol = optimize.root(func2, x0, args=(-2, -0.1, -0.1, -2, 1.25, 0.2))
    x = sol.x
    ax.plot([x[0], x[2]], [x[1], x[3]], "-p")

    # ax.plot([0.4075, 1.25, 2.25], [0.193, 0.25, 0.3187], "-o", fillstyle="none")
    # ax.plot([7.984/20, 1.25, 45.11/20], [3.849/20, 0.25, 6.374/20], "s", fillstyle="none")


    # rhoA_0 = 1.25
    # rhoB_0 = 0.25
    # rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    # phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    # x0 = [0.4, 0.19, 2.25, 0.32, 0.5]
    # sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # x = sol.x
    # ax.plot([x[0], x[2]], [x[1], x[3]], "*")

    # sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # x = sol.x
    # ax.plot([x[0], x[2]], [x[1], x[3]], "p")


    # rhoA_0 = 4
    # rhoB_0 = 1.65
    # x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    # sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # x = sol.x
    # ax.plot([x[0], rhoA_0, x[2]], [x[1], rhoB_0, x[3]], "-o")
    
       
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    plt.show()
    plt.close()


def binodals_eq_Jm01(ax):
    etaAA = -2
    etaBB = -2
    etaAB = -0.1
    etaBA = -0.1
    bc = "tab:red"
    ls = "dashed"
    lw = 1

    # G + LA + LAB
    rhoA_0 = 2.25
    rhoB_0 = 1
    x0 = [0.438, 0.3424, 2.418, 0.4252, 3.00, 2.49, 0.546, 0.2843]
    sol = optimize.root(func3_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    plot_tri(ax, pA1, pB1, c=bc, linestyle=ls, lw=lw)

    # G + LB + LAB
    rhoA_0 = 1
    rhoB_0 = 2.25
    x0 = [0.342, 0.438, 0.425, 2.418,  2.495,  3.00, 0.546, 0.28]
    sol = optimize.root(func3_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA2 = [x[0], x[2], x[4]]
    pB2 = [x[1], x[3], x[5]]
    plot_tri(ax, pA2, pB2, c=bc, linestyle=ls, lw=lw)

    # G + LA
    rhoA_0 = 1.25
    rhoB_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0]) 
    rhoB_0 = 0.01
    n = 40
    rhoB0_arr = np.linspace(0.01, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # G + LB
    rhoB_0 = 1.25
    rhoA_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0])
    rhoA_0_min = 0.01
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)


    # G + LAB
    rhoA_0 = 3
    rhoB_0_min = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[2] - pB1[0]) / (pA1[2] - pA1[0]) 
    rhoB_0_max = pB2[0] + (rhoA_0 - pA2[0]) * (pB2[2] - pB2[0]) / (pA2[2] - pA2[0])
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[0], pB1[0], pA1[2], pB1[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)


    # LB + LAB
    rhoA_0 = 1.5
    rhoB_0_min = pB2[1] + (rhoA_0 - pA2[1]) * (pB2[2] - pB2[1]) / (pA2[2] - pA2[1])
    rhoB_0_max = 8.5
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[1], pB2[1], pA2[2], pB2[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # LA + LAB
    rhoB_0 = 1.5
    rhoA_0_min = pA1[1] + (rhoA_0 - pB1[1]) * (pA1[2] - pA1[1]) / (pB1[2] - pB1[1])
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)


def binodals_Jm01(ax=None):
    etaAA = -2
    etaBB = -2
    etaAB = -0.1
    etaBA = -0.1
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))
        flag_show = True
    else:
        flag_show = False

    bc = "k"
    bc_eq = "tab:grey"

    # G + LA + LAB
    rhoA_0 = 2.25
    rhoB_0 = 1
    x0 = [0.5, 0.25, 2.5, 0.4, 5, 3, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    plot_tri(ax, pA1, pB1, c=bc)

    ax.fill(pA1, pB1, c="tab:orange", alpha=0.2)

    xx = np.array([29.4862, 3.8145, 23.1582]) / 10
    yy = np.array([23.8838, 2.9116, 3.6422]) / 10
    ms = 5
    ax.plot(xx, yy, "ko", fillstyle="none", ms=ms)
    ax.plot(yy, xx, "ko", fillstyle="none", ms=ms)

    # G + LB + LAB
    rhoA_0 = 1
    rhoB_0 = 2.25
    x0 = [0.25, 0.5, 0.4, 2.5, 3, 5, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA2 = [x[0], x[2], x[4]]
    pB2 = [x[1], x[3], x[5]]
    plot_tri(ax, pA2, pB2, c=bc)

    ax.fill(pA2, pB2, c="tab:orange", alpha=0.2)

    # G + LA
    rhoA_0 = 1.25
    rhoB_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0]) 
    rhoB_0 = 0.001
    n = 40
    rhoB0_arr = np.linspace(rhoB_0, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    edges = np.array([
        [3.8287, 0.4497, 23.0837, 0.5599],
        [3.8271, 1.3544, 23.0802, 1.6798],
        [3.8240, 2.2620, 23.0299, 2.7907],
        [3.8104, 2.6964, 23.0823, 3.3723]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 1], edges[:, 0], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)

    # G + LB
    rhoB_0 = 1.25
    rhoA_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0])
    rhoA_0_min = 0.001
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)


    # G + LAB
    rhoA_0 = 3
    rhoB_0_min = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[2] - pB1[0]) / (pA1[2] - pA1[0]) 
    rhoB_0_max = pB2[0] + (rhoA_0 - pA2[0]) * (pB2[2] - pB2[0]) / (pA2[2] - pA2[0])
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[0], pB1[0], pA1[2], pB1[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)


    # LB + LAB
    rhoA_0 = 1.5
    rhoB_0_min = pB2[1] + (rhoA_0 - pA2[1]) * (pB2[2] - pB2[1]) / (pA2[2] - pA2[1])
    rhoB_0_max = 8.5
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[1], pB2[1], pA2[2], pB2[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # LA + LAB
    rhoB_0 = 1.5
    rhoA_0_min = pA1[1] + (rhoA_0 - pB1[1]) * (pA1[2] - pA1[1]) / (pB1[2] - pB1[1])
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    rhoA_0 = 2
    rhoB_0 = 1.5
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.376, 0.211, 4.293, 3.302, 0.5]
    # x0 = [0.376, 0.211, 3, 2, 0.5]

    # sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))

    # sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # sol = optimize.root(func2, x0, args=(-2, -0.1, -0.1, -2, rhoA_0, rhoB_0))
    # x = sol.x
    # ax.plot([x[0], x[2]], [x[1], x[3]], "-p")


    # x0 = [3.7808/10, 1.5073/10, 23.3153/10, 2.6195/10, 0.5]
    # sol = optimize.root(func2, x0, args=(-2, -0.1, -0.1, -2, 1.25, 0.2))
    # x = sol.x
    # ax.plot([x[0], x[2]], [x[1], x[3]], "-p")

    # ax.plot([0.4075, 1.25, 2.25], [0.193, 0.25, 0.3187], "-o", fillstyle="none")
    # ax.plot([7.984/20, 1.25, 45.11/20], [3.849/20, 0.25, 6.374/20], "s", fillstyle="none")


    rhoA_0 = 1.25
    rhoB_0 = 0.1
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))


    # rhoA_0 = 4
    # rhoB_0 = 1.65
    # x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    # sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # x = sol.x
    # ax.plot([x[0], rhoA_0, x[2]], [x[1], rhoB_0, x[3]], "-o")
    
    
    binodals_eq_Jm01(ax)
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 3.5)


    ### Simulation data

    ## G + LAB

    edges = np.array([
        [3.2862, 3.2937, 26.0470, 25.9279],
        [3.426, 3.158, 26.63, 25.12],
        [3.565, 3.088, 27.786, 24.709],
        [3.73, 2.97, 28.49, 24.03]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 1], edges[:, 0], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)

    ## LA + LAB
    edges = np.array([
        [23.8438, 3.6572, 30.1036, 23.9658],
        [24.76, 3.69, 31.1, 23.87],
        [25.72, 3.655, 32.24, 24.0],
        [26.65, 3.69, 33.22, 23.90]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 1], edges[:, 0], "o", fillstyle="none", c=line.get_c(), ms=ms)

    edges = np.array([
        [29.43, 23.9, 23.31, 3.74, 3.87, 2.91]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 4], edges[:, 5], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 5], edges[:, 4], "o", fillstyle="none", c=line.get_c(), ms=ms)
    
    x, y = edges[0].reshape(3, 2).T
    plot_tri(ax, x, y, c="tab:grey", linestyle=":")
    plot_tri(ax, y, x, c="tab:grey", linestyle=":")


    # ax.plot(1.5, 1.5, "ro")
    # ax.plot(1.55, 1.45, "ro")
    # ax.plot(1.6, 1.4, "ro")
    # ax.plot(1.65, 1.35, "ro")


    # ax.plot(2.7, 1.4, "ro")
    # ax.plot(2.8, 1.4, "ro")
    # ax.plot(2.9, 1.4, "ro")
    # ax.plot(3.0, 1.4, "ro")


    if flag_show:
        plt.show()
        plt.close()


def binodals_eq_Jp01(ax):
    etaAA = -2
    etaBB = -2
    etaAB = 0.1
    etaBA = 0.1
    bc = "tab:red"
    ls = "dashed"
    lw = 1

    # G + LA + LAB
    rhoA_0 = 1.3
    rhoB_0 = 1.3
    x0 = [0.4, 0.4, 2., 0.3, 0.2, 2.5, 0.33, 0.33]
    sol = optimize.root(func3_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    plot_tri(ax, pA1, pB1, c=bc, linestyle=ls, lw=lw)

    # G + LB + LAB
    rhoA_0 = 3
    rhoB_0 = 3
    x0 = [3.5, 3.5, 0.5, 4, 4, 0.5, 0.33, 0.33]

    sol = optimize.root(func3_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA2 = [x[0], x[2], x[4]]
    pB2 = [x[1], x[3], x[5]]
    plot_tri(ax, pA2, pB2, c=bc, linestyle=ls, lw=lw)

    # LA + LB
    n = 40
    rhoA0_arr = np.linspace(pA1[1]+pA1[2], pA2[1]+pA2[2], n) * 0.5
    rhoB0_arr = rhoA0_arr
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))
    x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoB_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # G + LA
    rhoA_0 = 1.25
    rhoB_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0]) 
    rhoB_0 = 0.001
    n = 40
    rhoB0_arr = np.linspace(rhoB_0, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]
    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # G + LB
    rhoB_0 = 1.25
    rhoA_0_max = pA1[0] + (rhoB_0 - pB1[0]) * (pA1[2] - pA1[0]) / (pB1[2] - pB1[0])
    rhoA_0_min = 0.001
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # LB + LAB
    rhoA_0 = 1.25
    rhoB_0_min = pB2[1] + (rhoA_0 - pA2[1]) * (pB2[0] - pB2[1]) / (pA2[0] - pA2[1])
    rhoB_0_max = 8.5
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[1], pB2[1], pA2[0], pB2[0], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)

    # LA + LAB
    rhoB_0 = 1.5
    rhoA_0_min = pA2[0] + (rhoB_0 - pB2[0]) * (pA2[2] - pA2[0]) / (pB2[2] - pB2[0])
    print(rhoA_0_min)
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[0], pB2[0], pA2[2], pB2[2], 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2_eq, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc, linestyle=ls, lw=lw)
    ax.plot(phi_A_l, phi_B_l, c=bc, linestyle=ls, lw=lw)


def binodals_Jp01(ax=None):
    etaAA = -2
    etaBB = -2
    etaAB = 0.1
    etaBA = 0.1
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))
        flag_show = True
    else:
        flag_show = False

    bc = "k"

    # G + LA + LAB
    rhoA_0 = 1.3
    rhoB_0 = 1.3
    x0 = [0.4, 0.4, 2., 0.3, 0.2, 2.5, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    plot_tri(ax, pA1, pB1, c=bc)
    ax.fill(pA1, pB1, c="tab:orange", alpha=0.2)

    # G + LB + LAB
    rhoA_0 = 3
    rhoB_0 = 3
    x0 = [3.5, 3.5, 0.5, 4, 4, 0.5, 0.33, 0.33]

    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA2 = [x[0], x[2], x[4]]
    pB2 = [x[1], x[3], x[5]]
    plot_tri(ax, pA2, pB2, c=bc)

    ax.fill(pA2, pB2, c="tab:orange", alpha=0.2)

    # LA + LB
    n = 40
    rhoA0_arr = np.linspace(pA1[1]+pA1[2], pA2[1]+pA2[2], n) * 0.5
    rhoB0_arr = rhoA0_arr
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))
    x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoB_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)
    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)


    # G + LA
    rhoA_0 = 1.25
    rhoB_0_max = pB1[0] + (rhoA_0 - pA1[0]) * (pB1[1] - pB1[0]) / (pA1[1] - pB1[0]) 
    rhoB_0 = 0.001
    n = 40
    rhoB0_arr = np.linspace(rhoB_0, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # G + LB
    rhoB_0 = 1.25
    rhoA_0_max = pA1[0] + (rhoB_0 - pB1[0]) * (pA1[2] - pA1[0]) / (pB1[2] - pB1[0])
    rhoA_0_min = 0.001
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # LB + LAB
    rhoA_0 = 1.25
    rhoB_0_min = pB2[1] + (rhoA_0 - pA2[1]) * (pB2[0] - pB2[1]) / (pA2[0] - pA2[1])
    rhoB_0_max = 8.5
    n = 50
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[1], pB2[1], pA2[0], pB2[0], 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # LA + LAB
    rhoB_0 = 1.5
    rhoA_0_min = pA2[0] + (rhoB_0 - pB2[0]) * (pA2[2] - pA2[0]) / (pB2[2] - pB2[0])
    print(rhoA_0_min)
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [pA2[0], pB2[0], pA2[2], pB2[2], 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    ax.fill(xx, yy, c="tab:blue", alpha=0.2)
    
    binodals_eq_Jp01(ax)
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 3.5)


    ## simulation data
    ## G + LA
    edges = np.array([
        [3.837, 0.569, 23.169, 0.45],
        [3.859, 1.65, 23.2, 1.35],
        [3.83, 2.72, 23.14, 2.23]
    ]) / 10
    ms = 5
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 1], edges[:, 0], "o", fillstyle="none", c=line.get_c(), ms=ms)


    ## LAB + LA
    edges = np.array([
        [24.3, 23.54, 29.95, 3.81],
        [25.26, 23.46, 30.97, 3.825],
        [26.21, 23.49, 32.11, 3.79],
        [27.13, 23.54, 33.2, 3.824]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 1], edges[:, 0], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 3], edges[:, 2], "o", fillstyle="none", c=line.get_c(), ms=ms)

    ## LA + LB
    edges = np.array([
        [24.01, 3.17, 3.17, 24.01],
        [24.85, 3.27, 3.27, 24.85],
        [25.67, 3.36, 3.36, 25.67],
        [26.67, 3.49, 3.49, 26.67],
        [27.36, 3.52, 3.52, 27.36],
        [28.16, 3.62, 3.62, 28.16]
    ]) / 10

    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)

    edges = np.array([
        [23.5, 23.5, 29.1, 3.73, 3.73, 29.1],
        [3.83, 3.83, 23.16, 3.10, 3.10, 23.16]
    ]) / 10
    line, = ax.plot(edges[:, 0], edges[:, 1], "ko", fillstyle="none", ms=ms)
    ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)
    ax.plot(edges[:, 4], edges[:, 5], "o", fillstyle="none", c=line.get_c(), ms=ms)
    
    x, y = edges[0].reshape(3, 2).T
    plot_tri(ax, x, y, c="tab:grey", linestyle=":")
    x, y = edges[1].reshape(3, 2).T
    plot_tri(ax, x, y, c="tab:grey", linestyle=":")

    # ax.plot(1.25, 0.05, "ro")
    # ax.plot(1.25, 0.15, "ro")
    # ax.plot(1.25, 0.25, "ro")
    # ax.plot(2.7, 1.4, "ro")
    # ax.plot(2.8, 1.4, "ro")
    # ax.plot(2.9, 1.4, "ro")
    # ax.plot(3.0, 1.4, "ro", fillstyle="none")
    # ax.plot(1.5, 1.5, "ro")
    # ax.plot(1.4, 1.4, "ro")
    # ax.plot(1.35, 1.35, "ro")
    # ax.plot(1.45, 1.45, "ro")
    # ax.plot(1.55, 1.55, "ro")
    # ax.plot(1.6, 1.6, "ro")


    if flag_show:
        plt.show()
        plt.close()


def binodals_Jp05(ax=None):
    etaAA = -2
    etaBB = -2
    etaAB = 0.5
    etaBA = -0.5
    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))
        flag_show = True
    else:
        flag_show = False

    bc = "k"

    # G + LA + LAB
    rhoA_0 = 1.3
    rhoB_0 = 1.2
    x0 = [0.4, 0.4, 1.9, 0.3, 0.2, 2.5, 0.33, 0.33]
    sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    pA1 = [x[0], x[2], x[4]]
    pB1 = [x[1], x[3], x[5]]
    plot_tri(ax, pA1, pB1, c=bc)
    # ax.fill(pA1, pB1, c="tab:orange", alpha=0.2)
    # print(pA1, pB1)

    # # G + LB + LAB
    # rhoA_0 = 3.2
    # rhoB_0 = 3
    # x0 = [3.5, 3.5, 0.5, 4.2, 4.2, 0.5, 0.33, 0.33]

    # sol = optimize.root(func3, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    # x = sol.x
    # pA2 = [x[0], x[2], x[4]]
    # pB2 = [x[1], x[3], x[5]]
    # plot_tri(ax, pA2, pB2, c=bc)

    # ax.fill(pA2, pB2, c="tab:orange", alpha=0.2)

    # # LA + LB
    # n = 40
    # rhoA0_arr = np.linspace(pA1[1]+pA1[2], pA2[1]+pA2[2], n) * 0.5
    # rhoB0_arr = rhoA0_arr
    # phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))
    # x0 = [pA1[1], pB1[1], pA1[2], pB1[2], 0.5]
    # for i, rhoB_0 in enumerate(rhoB0_arr):
    #     sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoB_0, rhoB_0))
    #     x = sol.x
    #     x0 = x
    #     phi_A_g[i] = x[0]
    #     phi_B_g[i] = x[1]
    #     phi_A_l[i] = x[2]
    #     phi_B_l[i] = x[3]

    # ax.plot(phi_A_g, phi_B_g, c=bc)
    # ax.plot(phi_A_l, phi_B_l, c=bc)
    # xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    # yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    # ax.fill(xx, yy, c="tab:blue", alpha=0.2)


    # G + LA
    rhoA_0 = 1
    # rhoB_0_max = 0.85
    rhoB_0_max = 0.185
    rhoB_0 = 0.001
    n = 100
    rhoB0_arr = np.linspace(rhoB_0, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.25, rhoB_0, 2.5, rhoB_0, 0.5]
    for i, rhoB_0 in enumerate(rhoB0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    # xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    # yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    # ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # edges = np.array([
    #     [3.8287, 0.4497, 23.0837, 0.5599],
    #     [3.8271, 1.3544, 23.0802, 1.6798],
    #     [3.8240, 2.2620, 23.0299, 2.7907],
    #     [3.8104, 2.6964, 23.0823, 3.3723]
    # ]) / 10
    # ms = 3
    # line, = ax.plot(edges[:, 0], edges[:, 1], "o", fillstyle="none", ms=ms)
    # ax.plot(edges[:, 2], edges[:, 3], "o", fillstyle="none", c=line.get_c(), ms=ms)

    # G + LB
    rhoB_0 = 1.5
    # rhoA_0_max = pA1[0] + (rhoB_0 - pB1[0]) * (pA1[2] - pA1[0]) / (pB1[2] - pB1[0])
    # rhoA_0_max = 0.5
    rhoA_0_max = 0.275
    rhoA_0_min = 0.001
    n = 40
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [rhoA_0_min, 0.25, rhoA_0_min, 2.5, 0.5]
    for i, rhoA_0 in enumerate(rhoA0_arr):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)


    # xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    # yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    # ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # LB + LAB
    rhoA_0 = 1.25
    # rhoB_0_min = 2.8
    rhoB_0_min = 4.35
    rhoB_0_max = 12
    n = 80
    rhoB0_arr = np.linspace(rhoB_0_min, rhoB_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [0.43310682, 7.01958637, 2.37336856, 18.84892495, 0.42102216]
    for i, rhoB_0 in enumerate(rhoB0_arr[::-1]):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)

    print(x0)
    # xx = np.hstack((phi_A_g, phi_A_l[::-1]))
    # yy = np.hstack((phi_B_g, phi_B_l[::-1]))
    # ax.fill(xx, yy, c="tab:blue", alpha=0.2)

    # LA + LAB
    rhoB_0 = 1.
    # rhoA_0_min = 1.5
    rhoA_0_min = 3.4
    rhoA_0_max = 8.5
    n = 50
    rhoA0_arr = np.linspace(rhoA_0_min, rhoA_0_max, n)
    phi_A_l, phi_B_l, phi_A_g, phi_B_g = np.zeros((4, n))

    x0 = [5.04446035, 2.57116437, 10.51456208, 0.47976788, 0.63171397]
    for i, rhoA_0 in enumerate(rhoA0_arr[::-1]):
        sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
        x = sol.x
        x0 = x
        phi_A_g[i] = x[0]
        phi_B_g[i] = x[1]
        phi_A_l[i] = x[2]
        phi_B_l[i] = x[3]

    ax.plot(phi_A_g, phi_B_g, c=bc)
    ax.plot(phi_A_l, phi_B_l, c=bc)
    

    ## CCB + G
    # rhoA_0 = 0.75
    rhoA_0 = 2.4
    rhoB_0 = 1.25

    x0 = [1.8, 2.1, 3, 0.45, 0.5]
    sol = optimize.root(func2, x0, args=(etaAA, etaAB, etaBA, etaBB, rhoA_0, rhoB_0))
    x = sol.x
    print(x)
    ax.plot([x[0], x[2]], [x[1], x[3]])
    ax.plot(rhoA_0, rhoB_0, "o")

    
    # binodals_eq_Jp01(ax)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(0, 6.5)

    if flag_show:
        plt.show()
        plt.close()


def plot_FIG2():
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True, constrained_layout=True)
    binodals_Jm01(ax1)
    binodals_Jp01(ax2)


    label_fs = "large"
    ax1.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize=label_fs)
    ax2.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize=label_fs)
    ax1.set_ylabel(r"$\bar{\rho}_B/\rho_0$", fontsize=label_fs)

    ax1.set_title("(a) Mutual attraction", fontsize=label_fs)
    ax2.set_title("(b) Mutual repulsion", fontsize=label_fs)

    fs="large"
    ax1.text(0.1, 0.1, "G", fontsize=fs)
    ax1.text(3, 0.1, "LA", fontsize=fs)
    ax1.text(0.05, 3.2, "LB", fontsize=fs)
    ax1.text(3, 3.2, "LAB", fontsize=fs)

    ax1.text(0.1, 1, "G+LB", fontsize=fs, rotation=90)
    ax1.text(1, 0.1, "G+LA", fontsize=fs, rotation=0)
    ax1.text(1, 3.2, "LB+LAB", fontsize=fs, rotation=0)
    ax1.text(3, 1, "LA+LAB", fontsize=fs, rotation=90)
    ax1.text(1, 1, "G+LAB", fontsize=fs, rotation=45)
    ax1.text(0.5, 1, "G+LB+LAB", fontsize=fs, rotation=60)
    ax1.text(1, 0.5, "G+LA+LAB", fontsize=fs, rotation=30)



    ax2.text(0.1, 0.1, "G", fontsize=fs)
    ax2.text(3, 0.1, "LA", fontsize=fs)
    ax2.text(0.05, 3.2, "LB", fontsize=fs)
    ax2.text(3, 3.2, "LAB", fontsize=fs)

    ax2.text(0.1, 1, "G+LB", fontsize=fs, rotation=90)
    ax2.text(1, 0.1, "G+LA", fontsize=fs, rotation=0)
    ax2.text(1, 3.2, "LB+LAB", fontsize=fs, rotation=0)
    ax2.text(3, 1, "LA+LAB", fontsize=fs, rotation=90)

    ax2.text(1.3, 1.3, "LA+LB", fontsize=fs, rotation=-45)
    ax2.text(1.7, 1.7, "LA+LB+LAB", fontsize=fs, rotation=-45)
    ax2.text(0.5, 0.5, "G+LA+LB", fontsize=fs, rotation=-45)

    plt.show()
    # plt.savefig("fig/PD_mutual_attr_rep.pdf")
    plt.close()

def get_psi(phi_A, phi_B, eta_AA, eta_AB, eta_BA, eta_BB):
    def f_int_RX(phi_X, eta_XX, eta_YX):
        v_XX = get_v_XY(phi_X, eta_XX)
        v_YX = get_v_XY(phi_X, eta_YX)
        v_YX_prime = get_dv_XY(v_YX, eta_YX)
        return v_YX_prime / v_YX * np.log(phi_X * v_XX)

    def f_int_A(phi_A, eta_AA, eta_BA):
        v_AA = get_v_XY(phi_A, eta_AA)
        v_BA = get_v_XY(phi_A, eta_BA)
        v_BA_prime = get_dv_XY(v_BA, eta_BA)
        return v_BA_prime / v_BA * np.log(phi_A * v_AA)

    RA_RB = get_RA_RB(phi_A, phi_B, eta_BA, eta_AB)
    intA, err = integrate.quad(f_int_A, rho_min, phi_A, args=(eta_AA, eta_BA))
    intB, err = integrate.quad(f_int_RX, rho_min, phi_B, args=(eta_BB, eta_AB))
    return intA + intB + RA_RB


if __name__ == "__main__":
    # plot_FIG2()

    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    from matplotlib import image as mpimg

    fig = plt.figure(figsize=(10, 7.6), layout="constrained")
    subfigs = fig.subfigures(2, 1, wspace=0.001, hspace=0.01, height_ratios=[5, 2.8])
    (ax1, ax2) = subfigs[0].subplots(1, 2, sharex=True, sharey=True)
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True, constrained_layout=True)
    binodals_Jm01(ax1)
    binodals_Jp01(ax2)

    label_fs = "xx-large"
    ax1.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize=label_fs)
    ax2.set_xlabel(r"$\bar{\rho}_A/\rho_0$", fontsize=label_fs)
    ax1.set_ylabel(r"$\bar{\rho}_B/\rho_0$", fontsize=label_fs)

    ax1.set_title("(a) Mutual motility inhibition", fontsize=label_fs)
    ax2.set_title("(b) Mutual motility activation", fontsize=label_fs)

    fs="xx-large"
    ax1.text(0.1, 0.1, "G", fontsize=fs)
    ax1.text(3, 0.1, "LA", fontsize=fs)
    ax1.text(0.05, 3.2, "LB", fontsize=fs)
    ax1.text(3, 3.2, "LAB", fontsize=fs)

    ax1.text(0.1, 1, "G+LB", fontsize=fs, rotation=90)
    ax1.text(1, 0.15, "G+LA", fontsize=fs, rotation=0)
    ax1.text(1, 3.2, "LB+LAB", fontsize=fs, rotation=0)
    ax1.text(3, 1, "LA+LAB", fontsize=fs, rotation=90)
    ax1.text(0.8, 0.8, "G+LAB", fontsize=fs, rotation=45)
    ax1.text(0.5, 1, "G+LB+LAB", fontsize=fs, rotation=60)
    ax1.text(1, 0.5, "G+LA+LAB", fontsize=fs, rotation=30)



    ax2.text(0.1, 0.1, "G", fontsize=fs)
    ax2.text(3, 0.1, "LA", fontsize=fs)
    ax2.text(0.05, 3.2, "LB", fontsize=fs)
    ax2.text(3, 3.2, "LAB", fontsize=fs)

    ax2.text(0.1, 1, "G+LB", fontsize=fs, rotation=90)
    ax2.text(1, 0.1, "G+LA", fontsize=fs, rotation=0)
    ax2.text(1, 3.2, "LB+LAB", fontsize=fs, rotation=0)
    ax2.text(3, 1, "LA+LAB", fontsize=fs, rotation=90)

    ax2.text(1.3, 1.3, "LA+LB", fontsize=fs, rotation=-45)
    ax2.text(1.65, 1.65, "LA+LB+LAB", fontsize=fs, rotation=-45)
    ax2.text(0.5, 0.5, "G+LA+LB", fontsize=fs, rotation=-45)


    axes = subfigs[1].subplots(2, 3, sharex=True, sharey=True)

    fnames = ["G_LA.png", "G_LAB.png", "G_LA_LAB.png", "G_LA_LB.png", "LA_LB.png", "LA_LB_LAB.png"]
    titles = ["(c) G+LA", "(d) G+LAB", "(e) G+LA+LAB", "(f) G+LA+LB",
              "(g) LA+LB", "(h) LA+LB+LAB"]

    mx = [1.25, 1.5, 2, 1, 1.5, 2]
    my = [0.05, 1.5, 1, 1, 1.5, 2]
    mk = ["s", "p", ">", "P", "D", "*"]
    c = ["tab:blue", "tab:blue", "tab:orange", "tab:orange", "tab:blue", "tab:orange"]
    for i, ax in enumerate(axes.flat):
        fname = f"fig/snap/{fnames[i]}"
        im = mpimg.imread(fname)
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(titles[i], fontsize="xx-large")

        dx = 0.08
        ax_in = ax.inset_axes([0, 1-dx*3, dx, dx*3])
        ax_in.set_xticklabels([])
        ax_in.set_yticklabels([])
        ax_in.set_xticks([])
        ax_in.set_yticks([])

        ms = 8
        if i < 3:
            line, = ax1.plot(mx[i], my[i], mk[i], c=c[i], ms=ms)
        else:
            line, = ax2.plot(mx[i], my[i], mk[i], c=c[i], ms=ms)
        ax_in.plot(mx[i], my[i], mk[i], c=line.get_c(), ms=ms)

    # plt.show()
    plt.savefig("fig/PD_mutual_attr_rep_snaps.pdf")
    plt.close()