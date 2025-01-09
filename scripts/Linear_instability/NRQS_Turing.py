import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import schur, eigvals


def get_tilde_v_XY(eta_XY, rho_Y, bar_rho_Y, kappa=0.7):
    return 1 + kappa * np.tanh(eta_XY/kappa * (rho_Y - bar_rho_Y))

def get_tilde_v(eta_AA, eta_AB, eta_BA, eta_BB, rho_A, rho_B, bar_rho_A, bar_rho_B, kappa=0.7):
    drho_A = rho_A - bar_rho_A
    drho_B = rho_B - bar_rho_B
    inv_kappa = 1. / kappa
    v_AA = 1 + kappa * np.tanh(eta_AA * inv_kappa * drho_A)
    v_AB = 1 + kappa * np.tanh(eta_AB * inv_kappa * drho_B)
    v_BA = 1 + kappa * np.tanh(eta_BA * inv_kappa * drho_A)
    v_BB = 1 + kappa * np.tanh(eta_BB * inv_kappa * drho_B)
    return v_AA, v_AB, v_BA, v_BB

def get_v_A(rho_A, rho_B, bar_rho_A, bar_rho_B, etaAA, etaAB, bar_v_A=1., kappa=0.7):
    v_AA = get_tilde_v_XY(etaAA, rho_A, bar_rho_A, kappa)
    v_AB = get_tilde_v_XY(etaAB, rho_B, bar_rho_B, kappa) 
    return bar_v_A * v_AA * v_AB

def get_v_B(rho_A, rho_B, bar_rho_A, bar_rho_B, etaBA, etaBB, bar_v_B=1., kappa=0.7):
    v_BA = get_tilde_v_XY(etaBA, rho_A, bar_rho_A, kappa)
    v_BB = get_tilde_v_XY(etaBB, rho_B, bar_rho_B, kappa)
    return bar_v_B * v_BA * v_BB

def get_tilde_v_XY_derive(eta_XY, tilde_v_XY, kappa=0.7):
    return eta_XY * (1 - ((tilde_v_XY - 1)/kappa)**2)

def get_v0_omega(etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rho_A, bar_rho_B, bar_vA=1., bar_vB=1., kappa=0.7):
    v_AA, v_AB, v_BA, v_BB = get_tilde_v(etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rho_A, bar_rho_B, kappa)
    vA_0 = bar_vA * v_AA * v_AB
    vB_0 = bar_vB * v_BA * v_BB
    v_AA_deriv = get_tilde_v_XY_derive(etaAA, v_AA)
    v_AB_deriv = get_tilde_v_XY_derive(etaAB, v_AB)
    v_BA_deriv = get_tilde_v_XY_derive(etaBA, v_BA)
    v_BB_deriv = get_tilde_v_XY_derive(etaBB, v_BB)
    omega_AA = phiA * v_AA_deriv / v_AA
    omega_AB = phiA * v_AB_deriv / v_AB
    omega_BA = phiB * v_BA_deriv / v_BA
    omega_BB = phiB * v_BB_deriv / v_BB
    return vA_0, vB_0, omega_AA, omega_AB, omega_BA, omega_BB


def cal_Turing_prob(eta_AA, eta_BB, eta_AB, eta_BA=None, bar_rho_A=1, bar_rho_B=1, Dr_A=1, q_extr_thresh=None, p_range=[0, 4, 0, 4]):
    if eta_BA is None:
        eta_BA = -eta_AB

    phi_A = np.linspace(p_range[0], p_range[1], 400, endpoint=True)
    phi_B = np.linspace(p_range[2], p_range[3], 400, endpoint=True)

    pA, pB = np.meshgrid(phi_A, phi_B)
    Dr_B = 1

    vA_0, vB_0, omega_AA, omega_AB, omega_BA, omega_BB = get_v0_omega(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, bar_rho_A, bar_rho_B)

    sigma = vB_0 ** 2 / vA_0 ** 2 * Dr_A / Dr_B
    chi = omega_AB * omega_BA
    w1 = omega_AA + 1
    w2 = omega_BB + 1
    mask = np.logical_and(w1 + sigma * w2 > 0, w1 * w2 > chi)
    mask = np.logical_and(mask, w1 + w2 < 0)
    mask = np.logical_and(mask, (w1-w2)**2 + 4 * chi > 0)
    # mask = np.logical_and(omega_AA + sigma * omega_BB + (1+sigma)>0, omega_AA + omega_BB + 2 < 0)
    # mask = np.logical_and(mask, omega_AA * omega_BB - chi > 1)
    # mask = np.logical_and(mask, omega_AA * omega_BB - chi + omega_AA + omega_BB + 1 > 0)
    # mask = np.logical_and(mask, (omega_AA + omega_BB)**2 > 4 * (omega_AA * omega_BB - chi))
    
    if q_extr_thresh is not None:
        a, b, k_extr = np.zeros_like(pA), np.zeros_like(pA), np.zeros_like(pA)
        a = omega_AA * omega_BB - chi
        b = omega_AA + omega_BB
        k_extr[mask] = -b[mask]/ (2 * a[mask])
        q_extr = np.sqrt(20/3 * (1-k_extr))
        mask = np.logical_and(mask, q_extr < q_extr_thresh)
    M = np.zeros_like(pA)
    M[mask] = 1
    return 1 - np.mean(M)



def opt(x0, q_plus_thresh=None, q_extr_thresh=None, p_range=[0, 4, 0, 4], equal_Dr=False):
    def func(x):
        return cal_Turing_prob(x[0], x[1], x[2], eta_BA=x[3], bar_rho_A=x[4], bar_rho_B=x[5], Dr_A=x[6], q_extr_thresh=q_extr_thresh, p_range=p_range)
    
    def func_equal_Dr(x):
        return cal_Turing_prob(x[0], x[1], x[2], eta_BA=x[3], bar_rho_A=x[4], bar_rho_B=x[5], Dr_A=1, q_extr_thresh=q_extr_thresh, p_range=p_range)

    if x0[2] > 0:
        eta_AB_bnds = (0, 5)
        eta_BA_bnds = (-5, 0)
    else:
        eta_AB_bnds = (-5, 0)
        eta_BA_bnds = (0, 5)

    if not equal_Dr:
        bnds = ((-5, 1), (-5, 1), eta_AB_bnds, eta_BA_bnds, (0.1, 3), (0.1, 3), (0.1, 10))
        res = minimize(func, x0, method='Nelder-Mead', bounds=bnds, tol=1e-8)
    else:
        bnds = ((-5, 1), (-5, 1), eta_AB_bnds, eta_BA_bnds, (0.1, 3), (0.1, 3))
        res = minimize(func_equal_Dr, x0, method='Nelder-Mead', bounds=bnds, tol=1e-8)

    print(res)
    PD_pA_pB(res.x)
    return res.x


def PD_pA_pB(para=None):
    if para is None:
        # eta_AA = -2.3
        # eta_BB = -0.365
        # eta_AB = 1.16
        # bar_rho_A = 1.588
        # bar_rho_B = 1.557
        # Dr_A = 1.418
        # eta_AA = -2.3
        # eta_BB = -0.4
        # eta_AB = 1.2
        # bar_rho_A = 1.6
        # bar_rho_B = 1.6
        # Dr_A = 1.5
        # eta_BA = -eta_AB

        # eta_AA = -1.606
        # eta_BB = -0.171
        # eta_AB = 0.380
        # eta_BA = -1.44
        # bar_rho_A = 1.942
        # bar_rho_B = 2.618
        # Dr_A = 1.92

        eta_AA = -1.526
        eta_BB = -0.0067
        eta_AB = 0.19328
        eta_BA = -1.605
        bar_rho_A = 1.8222
        bar_rho_B = 4.3573
        Dr_A = 2.306

        # eta_AA = -1.8
        # eta_BB = -0.7
        # eta_AB = 1.2
        # bar_rho_A = 1.5
        # bar_rho_B = 1
        # Dr_A = 1
        # eta_BA = -eta_AB


        # eta_AA = -1.975
        # eta_BB = -1.3716
        # eta_AB = 1.905
        # bar_rho_A = 1.153
        # bar_rho_B = 1.2946
        # Dr_A = 1
        # eta_BA = -eta_AB

        # eta_AA = -2.149
        # eta_BB = -1.620
        # eta_AB = 1.759
        # eta_BA = -0.563
        # bar_rho_A = 1.466
        # bar_rho_B = 1.293
        # Dr_A = 0.763 

        # eta_AA = -2.45
        # eta_BB = -1.24
        # eta_AB = 1.13
        # eta_BA = -1.186
        # bar_rho_A = 1.550
        # bar_rho_B = 1.808
        # Dr_A = 1
    else:
        if para.size == 7:
            eta_AA = para[0]
            eta_BB = para[1]
            eta_AB = para[2]
            eta_BA = para[3]
            bar_rho_A = para[4]
            bar_rho_B = para[5]
            Dr_A = para[6]
        elif para.size == 6:
            eta_AA = para[0]
            eta_BB = para[1]
            eta_AB = para[2]
            eta_BA = para[3]
            bar_rho_A = para[4]
            bar_rho_B = para[5]
            Dr_A = 1

    Dr_B = 1

    phi_A = np.linspace(0, 5, 800, endpoint=True)
    phi_B = np.linspace(0, 5, 800, endpoint=True)

    pA, pB = np.meshgrid(phi_A, phi_B)

    vA_0, vB_0, omega_AA, omega_AB, omega_BA, omega_BB = get_v0_omega(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, bar_rho_A, bar_rho_B)

    sigma = vB_0 ** 2 / vA_0 ** 2 * Dr_A / Dr_B
    chi = omega_AB * omega_BA

    extent = [0, phi_A.max(), 0, phi_B.max()]
    mask = np.logical_and(omega_AA + sigma * omega_BB + (1+sigma)>0, omega_AA + omega_BB + 2 < 0)
    mask = np.logical_and(mask, omega_AA * omega_BB - chi > 1)
    mask = np.logical_and(mask, omega_AA * omega_BB - chi + omega_AA + omega_BB + 1 > 0)
    mask = np.logical_and(mask, (omega_AA + omega_BB)**2 > 4 * (omega_AA * omega_BB - chi))
    # M[mask] = 1
    a = omega_AA[mask] * omega_BB[mask] - chi[mask]
    b = omega_AA[mask] + omega_BB[mask]
    # q_plus = (-b + np.sqrt(b**2 - 4 * a))/ (2 * a)
    k_extr = (-b) / (2*a)
    q_extr = np.sqrt(20/3 * (1-k_extr))
    print("q_extr=", q_extr.min())
    M1 = np.zeros_like(pA)
    M1[mask] = q_extr

    M2 = np.zeros_like(pA)
    trM = -(1+sigma[mask]) - k_extr * (omega_AA[mask] + sigma[mask] * omega_BB[mask])
    detM = sigma[mask] * (1 + k_extr * omega_AA[mask]) * (1 + k_extr * omega_BB[mask]) - sigma[mask] * k_extr**2 * chi[mask]
    M2[mask] = q_extr**2 * vA_0[mask] ** 2 / (2 * Dr_A) * (trM + np.sqrt(trM**2 - 4 * detM))/2

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), constrained_layout=True, sharex=True, sharey=True)
    im1 = ax1.imshow(M1, origin="lower", extent=extent)
    im2 = ax2.imshow(M2, origin="lower", extent=extent)

    plt.colorbar(im1, ax=ax1, orientation="horizontal")
    plt.colorbar(im2, ax=ax2, orientation="horizontal")

    plt.show()
    plt.close()
    print(np.mean(M1))


class DMatrix:
    def __init__(self, etaAA, etaAB, etaBA, etaBB, phiA, phiB, Dr_A, Dr_B=None, bar_rhoA=1., bar_rhoB=1., bar_vA=1., bar_vB=1., kappa=0.7):
        vA_0, vB_0, omega_AA, omega_AB, omega_BA, omega_BB = get_v0_omega(
            etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rhoA, bar_rhoB, bar_vA, bar_vB, kappa)
        self.Dr_A = Dr_A
        if Dr_B is None:
            self.Dr_B = Dr_A
        else:
            self.Dr_B = Dr_B
        self.vA2_over_16Dr = vA_0 **2 / (16 * self.Dr_A)
        self.vB2_over_16Dr = vB_0 ** 2 / (16 *self.Dr_B)
        self.M = np.zeros((4, 4), complex)
        self.M[0, 1] = -1j * vA_0
        self.M[1, 0] = -0.5j * vA_0 * (1 + omega_AA)
        # self.M[1, 1] = 
        self.M[1, 2] = -0.5j * vA_0 * omega_AB
        self.M[2, 3] = -1j * vB_0
        self.M[3, 0] = -0.5j * vB_0 * omega_BA
        self.M[3, 2] = -0.5j * vB_0 * (1 + omega_BB)
        # self.M[3, 3] =

    def get_M(self, q):
        M = self.M.copy()
        M[1, 1] = -(self.Dr_A / q + self.vA2_over_16Dr * q)
        M[3, 3] = -(self.Dr_B / q + self.vB2_over_16Dr * q)
        return M
    
    def get_lambda(self, q):
        T, Z = schur(self.get_M(q), output="complex")
        eigen_values = np.array([T[0, 0], T[1, 1], T[2, 2], T[3, 3]]) * q
        return eigen_values
    
    def get_lambda_q(self, q_arr):
        lambda_arr = np.array([self.get_lambda(q) for q in q_arr])
        return lambda_arr
    
    def get_max_lambda(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr)
        max_lambda = lambda_arr.max(axis=1)
        return max_lambda.max()
    
    def get_lambda_max_lambda_q0(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr)
        max_lambda = lambda_arr.max(axis=1)
        return max_lambda.max(), max_lambda[0].real

    def get_lambda_max_Re_Im_q0(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr).max(axis=1)
        max_Re_idx = lambda_arr.real.argmax()
        max_Im_idx = np.abs(lambda_arr.imag).argmax()
        max_Re = lambda_arr[max_Re_idx].real
        if lambda_arr[max_Im_idx].real > 0:
            max_Im = np.abs(lambda_arr[max_Im_idx].imag)
        else:
            max_Im = np.abs(lambda_arr[max_Re_idx].imag)
        return max_Re, max_Im, lambda_arr[0].real



class DMatrix_w_surface_tension:
    def __init__(self, etaAA, etaAB, etaBA, etaBB, phiA, phiB, Dr_A, Dr_B=None, bar_rhoA=1., bar_rhoB=1., bar_vA=1., bar_vB=1., kappa=0.7, w_v2_over_q=True):
        vA_0, vB_0, omega_AA, omega_AB, omega_BA, omega_BB = get_v0_omega(
            etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rhoA, bar_rhoB, bar_vA, bar_vB, kappa)
        self.Dr_A = Dr_A
        if Dr_B is None:
            self.Dr_B = Dr_A
        else:
            self.Dr_B = Dr_B
        if w_v2_over_q:
            self.vA2_over_16Dr = vA_0 **2 / (16 * self.Dr_A)
            self.vB2_over_16Dr = vB_0 ** 2 / (16 *self.Dr_B)
        else:
            self.vA2_over_16Dr = 0
            self.vB2_over_16Dr = 0
        self.M = np.zeros((4, 4), complex)
        self.M[0, 1] = -1j * vA_0
        # self.M[1, 0] = -0.5j * vA_0 * (1 + omega_AA)
        # self.M[1, 1] = 
        # self.M[1, 2] = -0.5j * vA_0 * omega_AB
        self.M[2, 3] = -1j * vB_0
        # self.M[3, 0] = -0.5j * vB_0 * omega_BA
        # self.M[3, 2] = -0.5j * vB_0 * (1 + omega_BB)
        # self.M[3, 3] =
        self.omega = [omega_AA, omega_AB, omega_BA, omega_BB]
        self.vA_0 = vA_0
        self.vB_0 = vB_0

    def get_M(self, q):
        M = self.M.copy()
        M[1, 1] = -(self.Dr_A / q + self.vA2_over_16Dr * q)
        M[3, 3] = -(self.Dr_B / q + self.vB2_over_16Dr * q)
    
        k = 1-3/20 * q**2
        self.M[1, 0] = -0.5j * self.vA_0 * (1 + self.omega[0] * k)
        self.M[1, 2] = -0.5j * self.vA_0 * self.omega[1] * k
        self.M[3, 0] = -0.5j * self.vB_0 * self.omega[2] * k
        self.M[3, 2] = -0.5j * self.vB_0 * (1 + self.omega[3] * k)

        return M
    
    def get_lambda(self, q):
        T, Z = schur(self.get_M(q), output="complex")
        eigen_values = np.array([T[0, 0], T[1, 1], T[2, 2], T[3, 3]]) * q
        return eigen_values
    
    def get_lambda_q(self, q_arr):
        lambda_arr = np.array([self.get_lambda(q) for q in q_arr])
        return lambda_arr
    
    def get_max_lambda(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr)
        max_lambda = lambda_arr.max(axis=1)
        return max_lambda.max()
    
    def get_lambda_max_lambda_q0(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr)
        max_lambda = lambda_arr.max(axis=1)
        return max_lambda.max(), max_lambda[0].real

    def get_lambda_max_Re_Im_q0(self, q_arr):
        lambda_arr = self.get_lambda_q(q_arr).max(axis=1)
        max_Re_idx = lambda_arr.real.argmax()
        max_Im_idx = np.abs(lambda_arr.imag).argmax()
        max_Re = lambda_arr[max_Re_idx].real
        if lambda_arr[max_Im_idx].real > 0:
            max_Im = np.abs(lambda_arr[max_Im_idx].imag)
        else:
            max_Im = np.abs(lambda_arr[max_Re_idx].imag)
        return max_Re, max_Im, lambda_arr[0].real
    

def plot_lambda(ax, q, lamb, label=None, linestyle="-"):
    ax[0].plot(q, lamb.real, linestyle=linestyle)
    ax[1].plot(q, np.abs(lamb.imag), linestyle=linestyle)
    # ax[1].set_yscale("log")


def plot_lambda_varied_Dr(pA, pB, para=None):
    if para is None:
        eta_AA = -1.8
        eta_BB = -0.7
        eta_AB = 1.2
        eta_BA = -eta_AB
        bar_rho_A = 1.5
        bar_rho_B = 1
        pA = 124/80
        pB = 46/80
    else:
        if para.size == 7:
            eta_AA = para[0]
            eta_BB = para[1]
            eta_AB = para[2]
            eta_BA = para[3]
            bar_rho_A = para[4]
            bar_rho_B = para[5]
            Dr_A_over_B = para[6]
        elif para.size == 6:
            eta_AA = para[0]
            eta_BB = para[1]
            eta_AB = para[2]
            eta_BA = para[3]
            bar_rho_A = para[4]
            bar_rho_B = para[5]
            Dr_A_over_B = 1

    q = np.logspace(-2, 0.5 * np.log10(20/3), 1000, endpoint=False)
    # q = np.logspace(-2, 2, 1000, endpoint=False)


    fig, axes = plt.subplots(2, 3, figsize=(6, 8), sharex=True, sharey="row", constrained_layout=True)
    for i, Dr_B in enumerate([0.1, 1, 3]):
        Dr_A = Dr_B * Dr_A_over_B
        vA_0, vB_0, w_AA, w_AB, w_BA, w_BB = get_v0_omega(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, bar_rho_A, bar_rho_B)
        nu = vB_0 ** 2 / vA_0 ** 2 * Dr_A / Dr_B
        chi = w_AB * w_BA

        T = -(1+nu) - (w_AA + nu * w_BB)
        D = nu * (1+w_AA)*(1+w_BB) - nu * w_AB * w_BA
        Delta = T ** 2 - 4 * D
        lambda0 = (T + np.sqrt(Delta.astype(complex))) / 2 * q ** 2 * vA_0 ** 2 / (2 * Dr_A)
        # plot_lambda(axes[:, i], q, lambda0 * Dr_A)

        k = 1-3/20 * q**2
        T = -(1+nu) - k * (w_AA + nu * w_BB)
        D = nu * (1+k*w_AA)*(1+k*w_BB) - nu * k**2 * w_AB * w_BA
        Delta = T ** 2 - 4 * D
        lambda1 = (T + np.sqrt(Delta.astype(complex))) / 2 * q ** 2 * vA_0 ** 2 / (2 * Dr_A)
        plot_lambda(axes[:, i], q, lambda1 * Dr_A, linestyle="-")

        M2 = DMatrix(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, Dr_A, Dr_B, bar_rhoA=bar_rho_A, bar_rhoB=bar_rho_B)
        lambda2 = M2.get_lambda_q(q).max(axis=1)
        # plot_lambda(axes[:, i], q, lambda2 * Dr_A, linestyle="-.")


        M3 = DMatrix_w_surface_tension(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, Dr_A, Dr_B, bar_rhoA=bar_rho_A, bar_rhoB=bar_rho_B)
        lambda3 = M3.get_lambda_q(q).max(axis=1)
        plot_lambda(axes[:, i], q, lambda3 * Dr_A, linestyle="-.")


        M4 = DMatrix_w_surface_tension(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, Dr_A, Dr_B, bar_rhoA=bar_rho_A, bar_rhoB=bar_rho_B, w_v2_over_q=False)
        lambda4 = M4.get_lambda_q(q).max(axis=1)
        plot_lambda(axes[:, i], q, lambda4 * Dr_A, linestyle="--")

        axes[0, i].set_xscale("log")
        axes[0, i].axhline(0, c="k", linestyle="--")
        # axes[0, i].set_ylim(-0.1, 0.06)
        axes[1, i].set_ylim(-0.005, 0.1)
    plt.show()
    plt.close()


def plot_lambda_varied_phi(para=None):
    if para is None:
        eta_AA = -2.3
        eta_BB = -0.4
        eta_AB = 1.2
        eta_BA = -eta_AB
        bar_rho_A = 1.6
        bar_rho_B = 1.6
        Dr_A = 1.5

        # eta_AA = -1.606
        # eta_BB = -0.171
        # eta_AB = 0.380
        # eta_BA = -1.44
        # bar_rho_A = 1.942
        # bar_rho_B = 2.618
        # Dr_A = 1.92

        # eta_AA = -1.8
        # eta_BB = -0.7
        # eta_AB = 1.2
        # bar_rho_A = 1.5
        # bar_rho_B = 1
        # Dr_A = 1
        # eta_BA = -eta_AB
    else:
        if para.size == 7:
            eta_AA = para[0]
            eta_BB = para[1]
            eta_AB = para[2]
            eta_BA = para[3]
            bar_rho_A = para[4]
            bar_rho_B = para[5]
            Dr_A = para[6]
    Dr_B = 1

    pA = 3.4
    fig, axes = plt.subplots(3, 1, figsize=(6, 6), sharex=True, constrained_layout=True)
    for pB in [3.8, 3.9, 4.0, 4.1, 4.5]:
        vA_0, vB_0, w_AA, w_AB, w_BA, w_BB = get_v0_omega(eta_AA, eta_AB, eta_BA, eta_BB, pA, pB, bar_rho_A, bar_rho_B)
        nu = vB_0 ** 2 / vA_0 ** 2 * Dr_A / Dr_B
        chi = w_AB * w_BA

        q = np.logspace(-2, 0.5 * np.log10(20/3), 1000, endpoint=False)
        k = 1-3/20 * q**2
        T = -(1+nu) - k * (w_AA + nu * w_BB)
        D = nu * (1+k*w_AA)*(1+k*w_BB) - nu * k**2 * w_AB * w_BA

        axes[0].plot(q, T, label=r"$%g$" % pB)
        axes[1].plot(q, D)
        axes[0].set_ylabel(r"tr$\,\mathbf{M}$")
        axes[1].set_ylabel(r"det$\,\mathbf{M}$")
        Delta = T ** 2 - 4 * D
        mask = Delta > 0
        lamb = T[mask] + np.sqrt(Delta[mask])
        line, = axes[2].plot(q[mask], lamb, "-", label=r"$%g$" % pB)
        # line, = axes[2].plot(q[mask], lamb, "-", label=r"$\phi_A=%g$" % pA)

        mask = Delta <= 0
        axes[2].plot(q[mask], T[mask], ":", c=line.get_c())
        for ax in axes:
            ax.axhline(0, c="k", linestyle="--")
            ax.set_xscale("log")
        # plt.ylim(-0.1, 0.1)
        axes[2].set_xlabel(r"$q$")
        axes[2].set_ylabel(r"$\Re(\lambda)$")
        # axes[0].legend(loc="lower left", title=r"$\phi_B/\rho_0=$")
        axes[0].legend(loc="lower left", title=r"$\phi_A/\rho_0=$")

    # plt.ylim()
    plt.show()
    plt.close()


if __name__ == "__main__":
    para = np.array([100, 100, -1, 1, 1.0, 1, 1])
    # PD_pA_pB(para)

    plot_lambda_varied_Dr(1., 1., para)




