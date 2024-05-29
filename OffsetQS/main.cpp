#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx = 48;
  double Ly = 48;

  double rho0 = 20;
  double rhoA = 5;
  double eta = 3;
  double Dt = 0.;
  double v0 = 1;
  double kappa = 0.9;
  double Dr = 5;
  double eps = 0.2;

  double h = 0.02;
  int n_step = 10000000;
  //int snap_interval = int(round(200 / h * 0.1));
  double snap_log_sep = 0.1;
  int log_interval = 10000;
  int seed = 1000;
  std::string ini_mode = "rand";  // should be "rand", "ordered", "bi-ordered" or "resume"

  int snap_interval = 1;
  if (snap_log_sep > 1) {
    snap_interval = int(snap_log_sep);
    snap_log_sep = -1;
  }
  int n_par = int(rhoA * Lx * Ly);

  typedef BiNode<QSP_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_QS = 1;
  double r_cut = r_QS + eps;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;

  // ini integrator
  EM_QS_iso integrator(h, Dt, Dr);

  // cal force
  AsymKernal kernal(r_QS, eps);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
    };
  auto f2 = [&kernal, &pdm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, pdm);
    };

  // set output
  char basename[255];
  char log_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "D:/code/Quorum_sensing/OffsetQS/data/";
#else
  char folder[] = "data/";
#endif
  snprintf(basename, 255, "L%g_%g_Dr%.3f_Dt%.3f_r%g_p%g_e%.3f_E%.3f_h%.3f_%d",
    Lx, Ly, Dr, Dt, rho0, rhoA, eta, eps, h, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int start = 0;
  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, n_step, snap_interval, start, h, snap_log_sep, gl_l, ini_mode);

  // set QS parameters
  QSP_2::set_QS_params(v0, kappa, rho0, eta);
  // ini particles
  ini_particles(p_arr, myran, ini_mode, n_par, gl_l, gsd);
  cl.create(p_arr);

  float* v_arr = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      double v = p_arr[i].cal_v();
      v_arr[i] = v;
      integrator.update(p_arr[i], pdm, myran, v);
    }
    kernal.reset_local_density(p_arr);
    cl.recreate(p_arr);
    gsd.dump(t, p_arr, v_arr);
    log.record(t);
  }
  delete[] v_arr;
}
