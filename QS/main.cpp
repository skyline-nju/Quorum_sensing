#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx = 20;
  double Ly = 5;
  double phi = 200;
  double rho0 = 100;
  double h = 0.01;
  int n_step = 100000;

  double v0 = 1;
  double eta = 0;
  double alpha = 2;
  double Dr = 0.02;
  double Dt = 0.;
  std::string ini_mode = "restart";  // should be "new" or "restart"

  int n_par = int(round(Lx * Ly * phi));

  typedef BiNode<QSP_2> node_t;
  Ranq2 myran(2);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 1;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;
  
  // ini integrator
  EM_QS_iso integrator(h, Dt, Dr);

  // cal force
  LinearDensityKernal kernal(r_cut);
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
  snprintf(basename, 255, "%g_%g_%g_%g_%.2f_%.2f_%.1f_%g_%g", Lx, Ly, phi, rho0, eta, alpha, v0, Dr, Dt);
  snprintf(log_file, 255, "%s.log", basename);
  snprintf(gsd_file, 255, "%s.gsd", basename);

  int snap_interval = 100;
  int log_interval = 100000;
  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode == "new") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
      if (i * 2 >= n_par) {
        p_arr[i].type_id = 1;
      }
    }
  } else {
    gsd.read_last_frame(p_arr);
  }

  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      double v = cal_v(v0, rho0, eta, alpha, p_arr[i]);
      integrator.update(p_arr[i], pdm, myran, v);
    }
    kernal.reset_local_density(p_arr);
    cl.recreate(p_arr);
    gsd.dump(t, p_arr);
    log.record(t);
  }
}