#include "NRQS.h"

void run(const Vec_2<double>& gl_l,
         double phiA, double phiB, double rho0,
         double eta, double J_AB, double J_BA, double kappa,
         double Dt, double Dr, double h, double v0,
         int n_step, int snap_dt,
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

  Vec_2<int> proc_size = Vec_2<int>(4, tot_proc/4);
  Grid_2 grid(gl_l, r_cut, proc_size, group_comm);
  PeriodicDomain_2 dm(gl_l, grid, proc_size, group_comm);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, phiA + phiB, 50.);

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
  node_t::set_QS_params(v0, kappa, rho0, eta, J_AB, J_BA, eta);
  auto one_par_move = [&integrator, &dm, &myran](node_t& p) {
    double v = p.cal_v();
    integrator.update(p, dm, myran, v);
  };

  // set output
  char basename[255];
  char log_file[255];
  char order_para_file[255];
  char gsd_file[255];
  char folder[] = "D:/code/Quorum_sensing/data/";
  snprintf(basename, 255, "L%g_%g_Dr%.3f_k%.2f_p%g_%g_r%g_e%.3f_J%.3f_%.3f_%lld",
    gl_l.x, gl_l.y, Dr, kappa, phiA, phiB, rho0, eta, J_AB, J_BA, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(order_para_file, 255, "%sop_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int gl_np = round(phiA * gl_l.x * gl_l.y) + round(phiB * gl_l.x * gl_l.y);
  int start = 0;  // may be reset as the last time step of the gsd file
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, gl_l, gl_np, ini_mode, group_comm);

  int log_dt = 10000;
  io::LogExporter log(log_file, n_step, log_dt, start, gl_np, group_comm);

  int OP_dt = 100;
  int flush_dt = snap_dt;
  io::OrderParaExporter op(order_para_file, n_step, OP_dt, start, flush_dt, gl_np, group_comm);

  // initialize particles and celllist
  ini(p_arr, dm, cl, phiA, phiB, ini_mode, myran, gsd, 50.);

  // run
  for (int t = 1; t <= n_step; t++) {
    kernal.reset_local_density(p_arr);
    cal_force(p_arr, cl, comm, for_all_pair_force);
    kernal.normalize(p_arr);
    integrate(p_arr, cl, one_par_move, comm);
    gsd.dump(t, p_arr);
    op.dump(t, p_arr);
    log.record(t);
  }

  if (my_rank == 0) {
    std::cout << "Finish simulation!" << std::endl;
  }
  MPI_Barrier(group_comm);

}


int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = 40;
  double Ly = 40;
  double phiA = 10;
  double phiB = 10;
  double rho0 = 10;
  double h = 0.01;
  int n_step = 10000;
  int snap_dt = 1000;

  double kappa = 0.7;
  double eta = 0;
  double J_AB = 0.1;
  double J_BA = 0.1;
  double Dr = 0.01;
  double Dt = 0.;
  double v0 = 1.;
  int seed = 124;

  std::string ini_mode = "rand"; // should be 'rand' or 'resume'

  Vec_2<double> gl_l(Lx, Ly);

  run(gl_l, phiA, phiB, rho0,
      eta, J_AB, J_BA, kappa,
      Dt, Dr, h, v0, 
      n_step, snap_dt,
      ini_mode, seed,
      MPI_COMM_WORLD);

  MPI_Finalize();

}
