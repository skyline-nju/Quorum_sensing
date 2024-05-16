#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);

  double rho0 = 20;
  double rhoA = 20;
  double eta = 3;
  double Dt = atof(argv[3]);
  double v0 = 1;
  double kappa = 0.9;
  double Dr = 10;

  double angle = 180;
  double h = atof(argv[4]);
  int n_step = 100000000;
  //int snap_interval = int(round(200 / h * 0.1));
  int snap_interval = 10000;
  int log_interval = 10000;
  int seed = 1000;

  std::string ini_mode = "rand";  // should be "rand", "ordered", "bi-ordered" or "resume"

  int n_par = int(rhoA * Lx * Ly);

  typedef BiNode<QSP_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 1;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;
  
  // ini integrator
  EM_QS_iso integrator(h, Dt, Dr);

  // cal force
  //VisionConeKernal kernal(r_cut, angle);
  // LinearDensityKernal kernal(r_cut);
  double eps = 0.1;
  AsymKernal kernal(r_cut - abs(eps), eps);
  //NRTorque_Rep_2 kernal(eta, k, r_cut);
  //NRTorque_QS_2 kernal(r_cut, q, rho0);

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
  char folder[] = "D:/code/Quorum_sensing/VCQS/data/";
#else
  char folder[] = "/mnt/d/code/Quorum_sensing/VCQS/data/";
#endif
  // char folder[] = "./";
  // snprintf(basename, 255, "L%g_%g_Dr%.6f_Dt%.3f_r%g_p%g_e%.3f_h%.6f_%d", 
  //          Lx, Ly, Dr, Dt, rho0, rhoA, eta, h, seed);
    snprintf(basename, 255, "L%g_%g_Dr%.6f_Dt%.3f_r%g_p%g_e%.3f_E%.3f_h%.6f_%d", 
           Lx, Ly, Dr, Dt, rho0, rhoA, eta, eps, h, seed);
  //snprintf(basename, 255, "L%g_%g_Dr%g_Dt%g_r%g_p%g_e%.3f_a%g_h%.6f_%d",
  //  Lx, Ly, Dr, Dt, rho0, rhoA, eta, angle, h, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode == "rand") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
    }
  } else if (ini_mode == "ordered") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, 0., gl_l, Vec_2<double>());
    }
  } else if (ini_mode == "bi-ordered") {
    for (int i = 0; i < n_par; i++) {
      if (i < n_par / 2) {
        p_arr.emplace_back(myran, 0., gl_l, Vec_2<double>());
      } else {
        p_arr.emplace_back(myran, M_PI, gl_l, Vec_2<double>());
      }
  }
  } else {
    gsd.read_last_frame(p_arr);
  }
  
  cl.create(p_arr);
  
  float* rho_local = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      double v = cal_v(v0, kappa, rho0, eta, p_arr[i]);
      rho_local[i] = p_arr[i].rho_local[0];
      integrator.update(p_arr[i], pdm, myran, v);
    }
    kernal.reset_local_density(p_arr);
    cl.recreate(p_arr);
    gsd.dump(t, p_arr, rho_local);
    log.record(t);
  }
  delete[] rho_local;
}
