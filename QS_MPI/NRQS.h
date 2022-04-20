#pragma once

#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "io2D.h"
#include "config.h"
#include "communicator2D.h"

#ifdef USE_MPI
#include "mpi.h"
#endif


template <typename TDomain, typename TPar, typename TRan>
void ini(std::vector<TPar>& p_arr, const TDomain& dm, CellListNode_2<TPar>& cl,
         double phiA, double phiB, const std::string& ini_mode,
         TRan& myran,
         io::Snap_GSD_2& gsd,
         double amplification) {
  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);
  double area_gl = dm.gl_l().x * dm.gl_l().y;
  double area = dm.l().x * dm.l().y;
  size_t n_A_gl = size_t(round(phiA * area_gl));
  size_t n_B_gl = size_t(round(phiB * area_gl));
  int n_gl = n_A_gl + n_B_gl;

  size_t n_max_per_core = size_t((phiA + phiB) * area * amplification);
  if (n_max_per_core > n_gl) {
    n_max_per_core = n_gl;
  }
  p_arr.reserve(n_max_per_core);

  double* x = new double[n_gl]{};
  double* y = new double[n_gl]{};
  double* theta = new double[n_gl]{};
  uint8_t* type_id = new uint8_t[n_gl]{};

  if (my_rank == 0) {
    if (ini_mode == "rand") {
      for (size_t i = 0; i < n_gl; i++) {
        x[i] = myran.doub() * dm.gl_l().x;
        y[i] = myran.doub() * dm.gl_l().y;
        theta[i] = M_PI * myran.doub() * 2;
        if (i < n_A_gl) {
          type_id[i] = 0;
        } else {
          type_id[i] = 1;
        }
      }
      std::cout << "Create " << n_gl << " particles with random pos and ori" << std::endl;
    } else if (ini_mode == "resume") {
      gsd.read_last_frame(x, y, theta, type_id);
    } else {
      std::cout << "ini_mode need be one of rand or resume" << std::endl;
      exit(1);
    }

  }
  MPI_Barrier(dm.comm());
  MPI_Bcast(x, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(y, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(theta, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(type_id, n_gl, MPI_UINT8_T, 0, dm.comm());

  for (size_t i = 0; i < n_gl; i++) {
    if (dm.contain_particle(x[i], y[i])) {
      size_t j = p_arr.size();
      p_arr.push_back(TPar());
      p_arr[j].pos.x = x[i];
      p_arr[j].pos.y = y[i];
      p_arr[j].u.x = std::cos(theta[i]);
      p_arr[j].u.y = std::sin(theta[i]);
      p_arr[j].type_id = type_id[i];
    }
  }

  int my_np = p_arr.size();
  int tot_np;
  MPI_Reduce(&my_np, &tot_np, 1, MPI_INT, MPI_SUM, 0, dm.comm());
  if (my_rank == 0) {
    if (tot_np != n_gl) {
      std::cout << "Error when creating particles with tot_np=" << tot_np
        << " while n_gl=" << n_gl << std::endl;
      exit(1);
    }
  }
  cl.create(p_arr);

  MPI_Barrier(dm.comm());
  delete[] x;
  delete[] y;
  delete[] theta;
  delete[] type_id;
}

void run(const Vec_2<double>& gl_l,
         double phiA, double phiB, double rho0,
         double eta, double J_AB, double J_BA, double kappa,
         double Dt, double Dr, double h, double v0,
         int n_step, int snap_dt,
         const std::string& ini_mode,
         unsigned long long seed,
         MPI_Comm group_comm);