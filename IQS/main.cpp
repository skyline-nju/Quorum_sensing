#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"

int main(int argc, char* argv[]) {
  double Lx =40;
  double Ly = 5;

  int log_interval = 10000;
  double Dt = 0.01;
  double v0 = 1;
  double kappa = 0.9;

#ifdef _MSC_VER
  double rho0 = 10;
  double rhoA = 10;
  double eta = 6;
  double Dr = 0.;
  double gamma = 1;
  double h = 0.01;
  int n_step = 10000000;
  int snap_interval = 10000;
  int seed = 3002;
  std::string ini_mode = "bi-ordered";  // should be "rand", "ordered", or "resume"
#else
  double rho0 = atof(argv[1]);
  double rhoA = atof(argv[2]);
  double eta = atof(argv[3]);
  double Dr = atof(argv[4]);
  double gamma = atof(argv[5]);
  double h = atof(argv[6]);
  int n_step = atoi(argv[7]);
  int snap_interval = atoi(argv[8]);
  int seed = atoi(argv[9]);
  std::string ini_mode = argv[10]; // should be "new" or "restart"
#endif
  int n_par = int(rhoA * Lx * Ly);

  typedef BiNode<IQSP_2> node_t;
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
  LinearDensityKernal kernal(r_cut);
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
  char folder[] = "D:/code/Quorum_sensing/IQS/data/";
  // char folder[] = "./";
  snprintf(basename, 255, "L%g_%g_Dr%.3f_Dt%.3f_r%g_p%g_e%.3f_g%.3f_h%.6f_%d", 
           Lx, Ly, Dr, Dt, rho0, rhoA, eta, gamma, h, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, gl_l, ini_mode);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode != "resume") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
    }
    if (ini_mode == "ordered") {
      for (auto& p : p_arr) {
        p.u = Vec_2<double>(1, 0);
      }
    } else if (ini_mode == "bi-ordered") {
      for (int i = 0; i < p_arr.size(); i++) {
        if (i < p_arr.size() / 2) {
          p_arr[i].u = Vec_2<double>(1, 0);
        } else {
          p_arr[i].u = Vec_2<double>(-1, 0);
        }
      }
    }
  } else {
    gsd.read_last_frame(p_arr);
  }
  cl.create(p_arr);
  
  float* v = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      cal_v_underdamped(v0, kappa, rho0, eta, gamma, h, p_arr[i]);
      v[i] = p_arr[i].v;
      integrator.update(p_arr[i], pdm, myran, p_arr[i].v);
    }
    kernal.reset_local_density(p_arr);
    cl.recreate(p_arr);
    gsd.dump(t, p_arr, v);
    log.record(t);

    //if (t % 100 == 0) {
    //  double rho_m = 0;
    //  for (int j = 0; j < n_par; j++) {
    //    rho_m += rho_local[j];
    //  }
    //  std::cout << "t=" << t << ", rho=" << rho_m / n_par << std::endl;
    //}
  }
  delete[] v;
}
