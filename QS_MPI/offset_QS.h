#pragma once
#ifdef RUN_OFFSET_QS
#include "QS.h"

void run(const Vec_2<double>& gl_l,
  double phi, double rho0, double eps,
  double eta, double kappa,
  double Dt, double Dr, double h, double v0,
  int n_step, int snap_dt, double snap_log_dt,
  const std::string& ini_mode,
  unsigned long long seed,
  MPI_Comm group_comm);

#endif