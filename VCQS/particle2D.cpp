#include "particle2D.h"

double QSP_2::v0_;
double QSP_2::kappa_;
double QSP_2::rho0_;
double QSP_2::eta_;

void QSP_2::set_QS_params(double v0, double kappa, double rho0, double eta) {
  v0_ = v0;
  kappa_ = kappa;
  rho0_ = rho0;
  eta_ = eta / (rho0 * kappa);
}