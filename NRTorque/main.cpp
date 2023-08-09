#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx =20;
  double Ly = Lx;

  double rho0 = 10;
  double eta = 0.5;
  double k = 0.1;
  double k2 = 100;
  double Dt = 0.1;
  double Dr = 1;
  double v0 = 1;
  double sigma = 0.2;

  double h = 0.01;
  int n_step = 3000000;
  int snap_interval = 5000;
  int log_interval = 10000;
  int seed = 100;

  std::string ini_mode = "new";  // should be "new" or "restart"

  int n_par = int(rho0 * Lx * Ly);

  typedef BiNode<Par_NRTorque_Rep_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 1;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;
  
  // ini integrator
  EM_NRTorque_Rep integrator(h, Dt, Dr, v0);

  // cal force
  //NRTorque_Rep_2 kernal(eta, k, r_cut);
  NRTorque_Rep_hardcore_2 kernal(eta, k, r_cut, k2, sigma);

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
  char folder[] = "D:/code/Quorum_sensing/data/";
  // char folder[] = "./";
  snprintf(basename, 255, "L%g_v%g_k%.3f_%.3f_r%.3f_Dr%.3f_Dt%.3f_e%.3f_h%.4f_%d", 
           Lx, v0, k, k2, rho0, Dr, Dt, eta, h, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode == "new") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
    }
  } else {
    gsd.read_last_frame(p_arr);
  }
  
  cl.create(p_arr);
  
  float* rho_local = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    for (int i = 0; i < n_par; i++) {
      rho_local[i] = (p_arr[i].n_neighbor + 1) / M_PI;
      integrator.update(p_arr[i], pdm, myran);
    }
    cl.recreate(p_arr);
    gsd.dump(t, p_arr, rho_local);
    log.record(t);

    //if (t % 100 == 0) {
    //  double rho_m = 0;
    //  for (int j = 0; j < n_par; j++) {
    //    rho_m += rho_local[j];
    //  }
    //  std::cout << "t=" << t << ", rho=" << rho_m / n_par << std::endl;
    //}
  }
  delete[] rho_local;
}
