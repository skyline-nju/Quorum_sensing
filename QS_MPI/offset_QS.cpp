#include "offset_QS.h"

void run(const Vec_2<double>& gl_l,
         double phi, double rho0, double eps,
         double eta, double kappa,
         double Dt, double Dr, double h, double v0,
         int n_step, int snap_dt, double snap_log_dt,
         const std::string& ini_mode,
         unsigned long long seed,
         MPI_Comm group_comm) {
  typedef BiNode<QSP_1> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  Ranq1 myran(seed + my_rank);
  std::vector<node_t> p_arr;
  double r_QS = 1;
  double r_cut = r_QS + eps;

  Vec_2<int> proc_size = Vec_2<int>(2, 2);
  Grid_2 grid(gl_l, r_cut, proc_size, group_comm);
  PeriodicDomain_2 dm(gl_l, grid, proc_size, group_comm);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, phi, 60.);

  // cal force
  OffsetKernal kernal(r_QS, eps);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
    };
  auto f2 = [&kernal, &dm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, dm);
    };

  //auto f3 = [&kernal](node_t* p1, node_t* p2, const Vec_2<double>& offset) {
  //  kernal(*p1, *p2, p2->pos - p1->pos + offset);
  //  };
  auto for_all_pair_force = [&cl, &f1, &f2]() {
    cl.for_each_pair(f1, f2);
    //cl.for_each_pair_fast(f1, f3);
    };

  // ini integrator
  ABP_EM integrator(h, Dt, Dr);
  node_t::set_QS_params(v0, kappa, rho0, eta);
  auto one_par_move = [&integrator, &dm, &myran](node_t& p) {
    double v = p.cal_v();
    integrator.update(p, dm, myran, v);
    };

  // set output
  char basename[255];
  char log_file[255];
  char gsd_file[255];
  snprintf(basename, 255, "L%g_%g_Dr%.3f_Dt%.3f_r%g_p%g_e%.3f_E%.3f_h%.3f_%lld",
    gl_l.x, gl_l.y, Dr, Dt, rho0, phi, eta, eps, h, seed);
#ifndef _MSC_VER
  char folder[255] = "abc/";
#else
  char folder[] = "data/";
#endif

  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int log_dt = 10000;
  int gl_np = round(phi * gl_l.x * gl_l.y);
  int start = 0;
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, snap_log_dt, gl_l, gl_np, ini_mode, group_comm);
  io::LogExporter log(log_file, n_step, log_dt, start, gl_np, group_comm);

  // initialize particles and celllist
  ini(p_arr, dm, cl, phi, ini_mode, myran, gsd, 40.);

  // run
  for (int t = 1; t <= n_step; t++) {
    kernal.reset_local_density(p_arr);
    cal_force(p_arr, cl, comm, for_all_pair_force);
    kernal.normalize(p_arr);
    integrate2(p_arr, cl, one_par_move, comm);
    gsd.dump(t, p_arr);
    log.record(t);
  }

  if (my_rank == 0) {
    std::cout << "Finish simulation!" << std::endl;
  }
  MPI_Barrier(group_comm);
}


#ifdef RUN_OFFSET_QS
int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = 24;
  double Ly = 24;
  double rho0 = atof(argv[1]);
  double phi = atof(argv[2]);
  double eta = atof(argv[3]);
  double Dr = atof(argv[4]);
  double Dt = atof(argv[5]);

  double eps = atof(argv[6]);
  double h = atof(argv[7]);
  int n_step = atoi(argv[8]);
  double snap_log_sep = atof(argv[9]);
  int snap_dt = 1;
  if (snap_log_sep > 1) {
    snap_dt = int(snap_log_sep);
    snap_log_sep = -1;
  }

  double kappa = 0.9;
  double v0 = 1.;
  int seed = atoi(argv[10]);
  std::string ini_mode = argv[11]; // should be 'rand' or 'resume'
  Vec_2<double> gl_l(Lx, Ly);

  run(gl_l, phi, rho0,
    eps, eta, kappa,
    Dt, Dr, h, v0,
    n_step, snap_dt, snap_log_sep,
    ini_mode, seed,
    MPI_COMM_WORLD);

  MPI_Finalize();
}
#endif