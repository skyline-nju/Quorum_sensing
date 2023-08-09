#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx = 20;
  double Ly = Lx;

  double rho_c = 50;
  double M0 = 1;
  double k = 0.1;

  double h = 0.1;
  int n_step = 100000;

  double rho_gl = 20;

  double Dr = 0.1;
  double Dt = 0;
  double v0 = 1;
  int seed = 123;
  std::string ini_mode = "new";  // should be "new" or "restart"

  int n_par_A = int(rho_gl * Lx * Ly);
  int n_par_B = 0;
  int n_par = n_par_A + n_par_B;

  typedef BiNode<QSP_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 1;
  Grid_2 grid(gl_l, r_cut);
  PeriodicDomain_2 pdm(gl_l);
  CellListNode_2<node_t> cl(pdm, grid);
  std::vector<node_t> p_arr;
  
  // ini integrator
  EM_QS_iso_harmonic integrator(h, Dt, Dr, M0, k, gl_l);

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
  char folder[] = "D:/code/Quorum_sensing/data/";
  // char folder[] = "./";
  snprintf(basename, 255, "L%g_v%g_M%g_r%.3f_%.3f_k%.3f_Dr%.4f_%d", 
           Lx, v0, M0, rho_c, rho_gl, k, Dr, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int snap_interval = 10000;
  int log_interval = 10000;
  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode == "new") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
      if (i >= n_par_A) {
        p_arr[i].type_id = 1;
      }
    }
  } else {
    gsd.read_last_frame(p_arr);
  }
  
  cl.create(p_arr);
  
  float *rho_local = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      //double v = cal_v(v0, kappa, rho0, eta, J_AB, J_BA, p_arr[i]);
      double v = cal_v(v0, rho_c, p_arr[i]);
      integrator.update(p_arr[i], pdm, myran, v);
      rho_local[i] = p_arr[i].rho_local[0];
    }
    kernal.reset_local_density(p_arr);
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
