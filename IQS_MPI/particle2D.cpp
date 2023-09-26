#include "particle2D.h"

double QSP_2::v0_;
double QSP_2::kappa_;
double QSP_2::rho0_;
double QSP_2::rhoA_;
double QSP_2::rhoB_;
double QSP_2::eta_[2][2];

double IQSP_2::v0_;
double IQSP_2::kappa_;
double IQSP_2::rho0_;
double IQSP_2::eta_;
double IQSP_2::gamma_;

void QSP_2::set_QS_params(double v0, double kappa,
                          double rho0, double bar_rho_A, double bar_rho_B,
                          double phiA, double phiB,
                          double eta_AA, double eta_AB,
                          double eta_BA, double eta_BB) {
  v0_ = v0;
  kappa_ = kappa;
  rho0_ = rho0;
  rhoA_ = bar_rho_A;
  rhoB_ = bar_rho_B;
  eta_[0][0] = eta_AA / (rho0 * kappa);
  eta_[0][1] = eta_AB / (rho0 * kappa);
  eta_[1][0] = eta_BA / (rho0 * kappa);
  eta_[1][1] = eta_BB / (rho0 * kappa);
}


void IQSP_2::set_QS_params(double v0, double kappa, double rho0, double eta, double gamma) {
  v0_ = v0;
  kappa_ = kappa;
  rho0_ = rho0;
  eta_ = eta / (rho0 * kappa);
  gamma_ = gamma;
}

void IQSP_2::cal_v(double h) {
  double d_rho = rho_local - rho0_;
  double v_QS = v0_ * (1 + kappa_ * tanh(eta_ * d_rho));
  //p.v += h * gamma * (v_QS - p.v);
  //p.v = v_QS;

  // 2nd-order Runge-Kutta, https://blog.csdn.net/SanyHo/article/details/107017076
  double yn = v;
  double k1 = gamma_ * (v_QS - yn);
  double k2 = gamma_ * (v_QS - (yn + 0.5 * h * k1));
  v = yn + h * k2;
}

void IQSP_2::cal_v_4RK(double h) {
  // 4th-order Runge-Kutta
  double d_rho = rho_local - rho0_;
  double v_QS = v0_ * (1 + kappa_ * tanh(eta_ * d_rho));

  double yn = v;
  double k1 = gamma_ * (v_QS - yn);
  double k2 = gamma_ * (v_QS - (yn + 0.5 * h * k1));
  double k3 = gamma_ * (v_QS - (yn + 0.5 * h * k2));
  double k4 = gamma_ * (v_QS - (yn + h * k3));
  v = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}