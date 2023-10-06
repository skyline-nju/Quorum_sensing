#include "particle2D.h"

void IQSP_2::update_rho(double h, double gamma) {
  // 4th-order Runge-Kutta
  double yn = rho_local[1];
  double rho_QS = rho_local[0];
  double k1 = gamma * (rho_QS - yn);
  double k2 = gamma * (rho_QS - (yn + 0.5 * h * k1));
  double k3 = gamma * (rho_QS - (yn + 0.5 * h * k2));
  double k4 = gamma * (rho_QS - (yn + h * k3));
  rho_local[1] = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  rho_local[0] = 0.;
}
