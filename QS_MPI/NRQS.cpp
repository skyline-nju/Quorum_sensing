#include "NRQS.h"

void run(const Vec_2<double>& gl_l,
         double phiA, double phiB,
         double rho0, double bar_rho_A, double bar_rho_B,
         double eta, double J_AB, double J_BA, double kappa,
         double Dt, double Dr, double h, double v0,
         int n_step, int snap_dt, double snap_log_sep,
         const std::string& ini_mode,
         unsigned long long seed,
         MPI_Comm group_comm) {
  typedef BiNode<QSP_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;

  Vec_2<int> proc_size = Vec_2<int>(2, 2);
  Grid_2 grid(gl_l, r_cut, proc_size, group_comm);
  PeriodicDomain_2 dm(gl_l, grid, proc_size, group_comm);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, phiA + phiB, 60.);

  // cal force
  LinearDensityKernal kernal(r_cut);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
  };
  auto f2 = [&kernal, &dm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, dm);
  };

  auto f3 = [&kernal](node_t* p1, node_t* p2, const Vec_2<double>& offset) {
    kernal(*p1, *p2, p2->pos - p1->pos + offset);
  };
  auto for_all_pair_force = [&cl, &f1, &f3]() {
    //cl.for_each_pair(f1, f2);
    cl.for_each_pair_fast(f1, f3);
  };

  // ini integrator
  ABP_EM integrator(h, Dt, Dr);
  node_t::set_QS_params(v0, kappa, rho0, bar_rho_A, bar_rho_B, phiA, phiB, eta, J_AB, J_BA, eta);
  auto one_par_move = [&integrator, &dm, &myran](node_t& p) {
    double v = p.cal_v();
    integrator.update(p, dm, myran, v);
  };

  // set output
  char basename[255];
  char log_file[255];
  char gsd_file[255];
  char mc_file[255];
  snprintf(basename, 255, "L%g_%g_Dr%.3f_k%.2f_p%g_%g_r%g_%g_%g_e%.3f_J%.3f_%.3f_h%.3f_%lld",
    gl_l.x, gl_l.y, Dr, kappa, phiA, phiB, rho0, bar_rho_A, bar_rho_B, eta, J_AB, J_BA, h, seed);
#ifndef _MSC_VER
  char folder[255];
  char prefix[] = "/scratch03.local/yduan/QS_MPI_noself";
  if (gl_l.x == gl_l.y) {
    snprintf(folder, 255, "%s/L%g_r%g_%g_%g_h%.3f/", prefix, gl_l.x, rho0, bar_rho_A, bar_rho_B, h);
  } else {
    snprintf(folder, 255, "%s/L%g_%g_r%g_%g_%g_h%.3f/", prefix, gl_l.x, gl_l.y, rho0, bar_rho_A, bar_rho_B, h);
  }

#else
  char folder[] = "D:/code/Quorum_sensing/data/";
#endif

  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);
  snprintf(mc_file, 255, "%smc_%s.dat", folder, basename);


  int log_dt = 10000;
  int gl_np = round(phiA * gl_l.x * gl_l.y) + round(phiB * gl_l.x * gl_l.y);
  int start = 0;
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, snap_log_sep, gl_l, gl_np, ini_mode, group_comm);
  io::LogExporter log(log_file, n_step, log_dt, start, gl_np, group_comm);
  io::MassCenterExporter mc(mc_file, n_step, 100, start, gl_l, snap_log_sep, group_comm);

  // initialize particles and celllist
  ini(p_arr, dm, cl, phiA, phiB, ini_mode, myran, gsd, 40.);

  // run
  for (int t = 1; t <= n_step; t++) {
    kernal.reset_local_density(p_arr);
    cal_force(p_arr, cl, comm, for_all_pair_force);
    kernal.normalize(p_arr);
    integrate(p_arr, cl, one_par_move, comm);
    gsd.dump(t, p_arr);
    log.record(t);
    mc.record(t, p_arr);
  }

  if (my_rank == 0) {
    std::cout << "Finish simulation!" << std::endl;
  }
  MPI_Barrier(group_comm);

}


int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = 20;
  double Ly = 20;
  double phiA = atof(argv[1]);
  double phiB = atof(argv[2]);
  double rho0 = 20;
  double bar_rho_A = 20;
  double bar_rho_B = 20;
  double h = 1/3.;

  int n_step = atoi(argv[3]);
  double snap_log_sep = atof(argv[4]);
  int snap_dt = 1;
  if (snap_log_sep > 1) {
    snap_dt = int(snap_log_sep);
    snap_log_sep = -1;
  }

  double kappa = 0.7;
  double eta = atof(argv[5]);
  double J_AB = atof(argv[6]);
  double J_BA = atof(argv[7]);
  double Dr = atof(argv[8]);
  double Dt = 0.;
  double v0 = 1.;
  int seed = atoi(argv[9]);
  std::string ini_mode = argv[10]; // should be 'rand' or 'resume'
  Vec_2<double> gl_l(Lx, Ly);

  run(gl_l, phiA, phiB,
      rho0, bar_rho_A, bar_rho_B,
      eta, J_AB, J_BA, kappa,
      Dt, Dr, h, v0, 
      n_step, snap_dt, snap_log_sep,
      ini_mode, seed,
      MPI_COMM_WORLD);

  MPI_Finalize();

}
