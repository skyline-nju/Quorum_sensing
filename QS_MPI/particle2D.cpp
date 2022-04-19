#include "particle2D.h"

double QSP_2::v0_;
double QSP_2::kappa_;
double QSP_2::rho0_;
double QSP_2::eta_[2][2];

void QSP_2::set_QS_params(double v0, double kappa, double rho0,
                          double eta_AA, double eta_AB,
                          double eta_BA, double eta_BB) {
  v0_ = v0;
  kappa_ = kappa;
  rho0_ = rho0;
  eta_[0][0] = eta_AA / (rho0 * kappa);
  eta_[0][1] = eta_AB / (rho0 * kappa);
  eta_[1][0] = eta_BA / (rho0 * kappa);
  eta_[1][1] = eta_BB / (rho0 * kappa);
}
