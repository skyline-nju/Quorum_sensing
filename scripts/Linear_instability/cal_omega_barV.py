import numpy as np


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


def get_v_A(rho_A, rho_B, bar_rho_A, bar_rho_B, etaAA, etaAB, vA_0=1., kappa=0.7):
    v_AA = get_tilde_v_XY(etaAA, rho_A, bar_rho_A, kappa)
    v_AB = get_tilde_v_XY(etaAB, rho_B, bar_rho_B, kappa) 
    return vA_0 * v_AA * v_AB


def get_v_B(rho_A, rho_B, bar_rho_A, bar_rho_B, etaBA, etaBB, vB_0=1., kappa=0.7):
    v_BA = get_tilde_v_XY(etaBA, rho_A, bar_rho_A, kappa)
    v_BB = get_tilde_v_XY(etaBB, rho_B, bar_rho_B, kappa)
    return vB_0 * v_BA * v_BB


def get_tilde_v_XY_derive(eta_XY, tilde_v_XY, kappa=0.7):
    return eta_XY * (1 - ((tilde_v_XY - 1)/kappa)**2)


def get_bar_v_omega(etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rho_A, bar_rho_B, vA0=1., vB0=1., kappa=0.7):
    v_AA, v_AB, v_BA, v_BB = get_tilde_v(etaAA, etaAB, etaBA, etaBB, phiA, phiB, bar_rho_A, bar_rho_B, kappa)
    bar_vA = vA0 * v_AA * v_AB
    bar_vB = vB0 * v_BA * v_BB
    v_AA_deriv = get_tilde_v_XY_derive(etaAA, v_AA)
    v_AB_deriv = get_tilde_v_XY_derive(etaAB, v_AB)
    v_BA_deriv = get_tilde_v_XY_derive(etaBA, v_BA)
    v_BB_deriv = get_tilde_v_XY_derive(etaBB, v_BB)
    omega_AA = phiA * v_AA_deriv / v_AA
    omega_AB = phiA * v_AB_deriv / v_AB
    omega_BA = phiB * v_BA_deriv / v_BA
    omega_BB = phiB * v_BB_deriv / v_BB
    return bar_vA, bar_vB, omega_AA, omega_AB, omega_BA, omega_BB


if __name__ == "__main__":
    pass