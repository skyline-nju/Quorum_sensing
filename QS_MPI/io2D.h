#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "config.h"
#include "particle2D.h"
#include "gsd.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace io {

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

/**
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 */
class ExporterBase {
public:
  ExporterBase(int n_step, int sep, int start, MPI_Comm group_comm);

  bool need_export(const int i_step) {
    return (i_step - start_) % sep_ == 0;
  }

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
  MPI_Comm comm_;
  int my_rank_ = 0;
  int tot_proc_ = 1;
};

/**
 * @brief Exporter to output log
 *
 * Output the parameters after the initialization.
 * Output the beginning and endding time of the simulation.
 * Record time every certain time steps.
 */
class LogExporter : public ExporterBase {
public:
  LogExporter(const std::string& outfile,
              int n_step, int sep, int start,
              int np_gl, MPI_Comm group_comm);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  int step_count_ = 0;
};


class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
             int n_step, int sep, int start,
             const Vec_2<double>& gl_l,
             size_t n_par_gl,
             const std::string& open_flag,
             MPI_Comm group_comm);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr,
                         float* pos, uint32_t* type_id, float* v);

  uint64_t get_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, float* v);

  template <typename TFloat, typename TInt>
  void read(int i_frame, TFloat* x, TFloat* y, TFloat* theta, TInt* type_id);

  template <typename TFloat, typename TInt>
  void read_last_frame(TFloat* x, TFloat* y, TFloat* theta, TInt* type_id);

private:
  gsd_handle* handle_ = nullptr;
  Vec_2<double> gl_l_;
  size_t gl_np_;
};


template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr,
                                   float* pos,
                                   uint32_t* type_id,
                                   float* v) {
  int n_par = p_arr.size();
  int* n_par_arr = new int[tot_proc_]{};
  int* displs = new int[tot_proc_]{};
  MPI_Gather(&n_par, 1, MPI_INT, n_par_arr, 1, MPI_INT, 0, comm_);
  if (my_rank_ == 0) {
    size_t count = 0;
    for (int i = 0; i < tot_proc_; i++) {
      displs[i] += count;
      count += n_par_arr[i];
    }
    if (count != gl_np_) {
      std::cout << "Error, total particles = " << count << ", unequal to "
                << gl_np_ << std::endl;
      exit(1);
    }
  }

  float* my_pos = new float[n_par * 3];
  uint32_t* my_type_id = new uint32_t[n_par];
  float* my_v = new float[n_par];


  double half_Lx = gl_l_.x * 0.5;
  double half_Ly = gl_l_.y * 0.5;
  for (int j = 0; j < n_par; j++) {
    int j3 = j * 3;
    my_pos[j3    ] = p_arr[j].pos.x - half_Lx;
    my_pos[j3 + 1] = p_arr[j].pos.y - half_Ly;
    my_pos[j3 + 2] = p_arr[j].get_theta();
    my_type_id[j] = p_arr[j].type_id;
    my_v[j] = p_arr[j].cal_v();
  }

  MPI_Gatherv(my_type_id, n_par, MPI_UINT32_T,
              type_id, n_par_arr, displs, MPI_UINT32_T,
              0, comm_);
  MPI_Gatherv(my_v, n_par, MPI_FLOAT,
              v, n_par_arr, displs, MPI_FLOAT,
              0, comm_);

  for (int i = 0; i < tot_proc_; i++) {
    n_par_arr[i] *= 3;
    displs[i] *= 3;
  }

  MPI_Gatherv(my_pos, n_par * 3, MPI_FLOAT,
              pos, n_par_arr, displs, MPI_FLOAT,
              0, comm_);

  delete[] n_par_arr;
  delete[] displs;
  delete[] my_pos;
  delete[] my_type_id;
  delete[] my_v;
}

template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    float* pos = nullptr;
    uint32_t* type_id = nullptr;
    float* v = nullptr;
    
    if (my_rank_ == 0) {
      pos = new float[gl_np_ * 3];
      type_id = new uint32_t[gl_np_];
      v = new float[gl_np_];
    }

    get_data_from_par(p_arr, pos, type_id, v);
    uint32_t n_par = gl_np_;
    if (my_rank_ == 0) {
      uint64_t step = get_time_step();
      std::cout << "dump frame " << step << std::endl;
      gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
      gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
      gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
      gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);
      gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_par, 1, 0, v);
      gsd_end_frame(handle_);
    }
    delete[] pos;
    delete[] type_id;
    delete[] v;
  }
}

template <typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr, float* v) {
  if (need_export(i_step)) {
    uint32_t n_par = p_arr.size();
    float* pos = new float[n_par * 3];
    uint32_t *type_id = new uint32_t[n_par];
    get_data_from_par(p_arr, pos, type_id);
    uint64_t step = get_time_step();
  
    std::cout << "dump frame " << step << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);
    gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_par, 1, 0, v);
    gsd_end_frame(handle_);
    delete[] pos;
    delete []type_id;
  }
}

template<typename TFloat, typename TInt>
void io::Snap_GSD_2::read(int i_frame,
                          TFloat* x, TFloat* y, TFloat* theta,
                          TInt* type_id) {
  if (my_rank_ == 0) {
    uint32_t n_par;
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
    gsd_read_chunk(handle_, &n_par, chunk);
    if (n_par == gl_np_) {
      std::cout << "frame " << i_frame << ": find " << n_par << " particles" << std::endl;
    } else {
      std::cout << "Error for loading gsd file, wrong partiles number" << std::endl;
      exit(1);
    }
    float* pos = new float[size_t(n_par) * 3];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
    gsd_read_chunk(handle_, pos, chunk);
    uint32_t* type_id0 = new uint32_t[n_par];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/typeid");
    gsd_read_chunk(handle_, type_id0, chunk);

    double half_Lx = gl_l_.x * 0.5;
    double half_Ly = gl_l_.y * 0.5;
    for (int i = 0; i < n_par; i++) {
      x[i] = pos[i * 3 + 0] + half_Lx;
      y[i] = pos[i * 3 + 1] + half_Ly;
      theta[i] = pos[i * 3 + 2];
      type_id[i] = type_id0[i];
      tangle_1D(x[i], gl_l_.x);
      tangle_1D(y[i], gl_l_.y);
    }

    delete[] pos;
    delete[] type_id0;
  }
}

template<typename TFloat, typename TInt>
void io::Snap_GSD_2::read_last_frame(TFloat* x, TFloat* y,
                                     TFloat* theta, TInt* type_id) {
  if (my_rank_ == 0) {
    int nframes = gsd_get_nframes(handle_);
    if (nframes < 1) {
      std::cout << "Error, nframes=" << nframes << std::endl;
      exit(1);
    } else {
      read(nframes - 1, x, y, theta, type_id);
    }
  }
}

}

