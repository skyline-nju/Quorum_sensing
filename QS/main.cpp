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
  double phiA = atof(argv[1]);
  double phiB = atof(argv[2]);
  double rho0 = atof(argv[3]);
  double h = 0.01;

  double kappa = 0.7; 
  double eta = atof(argv[4]);
  double alpha = atof(argv[5]);
  double J_AB = alpha;
  double J_BA = -alpha;
  double Dr = atof(argv[6]);
  double Dt = 0;
  double v0 = 1;
  int seed = atoi(argv[7]);
  std::string ini_mode = argv[8];  // should be "new" or "restart"
  int n_step = atoi(argv[9]);
  int snap_interval = atoi(argv[10]);
  
  int n_par_A = int(round(Lx * Ly * phiA));
  int n_par_B = int(round(Lx * Ly * phiB));
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
  char op_file[255];
  char folder[] = "/scratch03.local/yduan/QS5/zero_eta/";
  // char folder[] = "./";
  snprintf(basename, 255, "L%g_%g_Dr%.3f_k%.2f_p%g_%g_r%g_e%.3f_J%.3f_%.3f_%d", 
           Lx, Ly, Dr, kappa, phiA, phiB, rho0, eta, J_AB, J_BA, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);
  snprintf(op_file, 255, "%sop_%s.dat", folder, basename);

  int log_interval = 10000;
  int start = 0;
  exporter::Snap_GSD_2 gsd(gsd_file, start, n_step, snap_interval, gl_l, ini_mode);

  exporter::LogExporter log(log_file, start, n_step, log_interval, n_par);

  int op_interval = 100;
  int op_flush_dt = snap_interval;
  exporter::OrderParaExporter op(op_file, n_step, op_interval, start, op_flush_dt, n_par);

  // ini particles
  p_arr.reserve(n_par);
  if (ini_mode == "rand") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
      if (i >= n_par_A) {
        p_arr[i].type_id = 1;
      }
    }
  } else if (ini_mode == "resume") {
    gsd.read_last_frame(p_arr);
  } else {
    std::cout << "Error, ini_mode must be 'rand' or 'resume'!" << std::endl;
    exit(1);
  }
  cl.create(p_arr);
  
  float *v_buf = new float[n_par];
  for (int t = 1; t <= n_step; t++) {
    kernal.reset_local_density(p_arr);
    cl.for_each_pair(f1, f2);
    kernal.normalize(p_arr);
    for (int i = 0; i < n_par; i++) {
      double v = cal_v(v0, kappa, rho0, eta, J_AB, J_BA, p_arr[i]);
      integrator.update(p_arr[i], pdm, myran, v);
      v_buf[i] = v;
    }
    cl.recreate(p_arr);
    gsd.dump(t, p_arr, v_buf);
    log.record(t);
    op.dump(t, p_arr, v_buf);
  }
  delete[] v_buf;
}
