#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
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
    if (log_sep_ < 0) {
      return i_step % sep_ == 0;
    } else {
      return need_export_log_scale(i_step);
    }
  }

  bool need_export_log_scale(const int i_step) {
    bool res = false;
    if (i_step + start_ == frames_[cur_frame_]) {
      cur_frame_++;
      res = true; 
    }
    return res;
  }


  void set_log_scale_frames(double h, double log_sep=0.1);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
  MPI_Comm comm_;
  int my_rank_ = 0;
  int tot_proc_ = 1;
  std::vector<int> frames_;
  int cur_frame_ = 0;
  double log_sep_ = -1;
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


class MassCenterExporter : public ExporterBase {
public:
  MassCenterExporter(const std::string& outfile,
    int n_step, int sep, int start,
    const Vec_2<double> &gl_l,
    int flush_sep, MPI_Comm group_comm);

  ~MassCenterExporter();

  template <typename TPar>
  void record(int i_step, const std::vector<TPar>& p_arr);
private:
  std::ofstream fout_;
  double Lx_;
  double Ly_;
  int flush_sep_;

};

template <typename TPar>
void MassCenterExporter::record(int i_step, const std::vector<TPar>& p_arr) {
  int cur_step = i_step + start_;
  if (need_export(cur_step)) {
    double sum_pos_vec[8]{}; // sum_A_sin_x, sum_A_cos_x, sum_A_sin_y, sum_A_cos_y,
                             // sum_B_sin_x, sum_B_cos_x, sum_B_sin_y, sum_B_cos_y
    for (const auto& p : p_arr) {
      double x = p.pos.x / Lx_ * M_PI * 2;
      double y = p.pos.y / Ly_ * M_PI * 2;
      if (p.type_id == 0) { // A particles
        sum_pos_vec[0] += sin(x);
        sum_pos_vec[1] += cos(x);
        sum_pos_vec[2] += sin(y);
        sum_pos_vec[3] += cos(y);
      } else { // B particles
        sum_pos_vec[4] += sin(x);
        sum_pos_vec[5] += cos(x);
        sum_pos_vec[6] += sin(y);
        sum_pos_vec[7] += cos(y);
      }
    }

    double gl_sum_pos_vec[8]{};
    MPI_Reduce(sum_pos_vec, gl_sum_pos_vec, 8, MPI_DOUBLE, MPI_SUM, 0, comm_);
    if (my_rank_ == 0) {
      double xm_A = atan2(gl_sum_pos_vec[0], gl_sum_pos_vec[1]) / (2 * M_PI) * Lx_;
      double ym_A = atan2(gl_sum_pos_vec[2], gl_sum_pos_vec[3]) / (2 * M_PI) * Ly_;
      double xm_B = atan2(gl_sum_pos_vec[4], gl_sum_pos_vec[5]) / (2 * M_PI) * Lx_;
      double ym_B = atan2(gl_sum_pos_vec[6], gl_sum_pos_vec[7]) / (2 * M_PI) * Ly_;
      fout_ << cur_step << "\t" << xm_A << "\t" << ym_A << "\t" << xm_B << "\t" << ym_B << "\n";
      if (cur_step % flush_sep_ == 0) {
        fout_ << std::flush;
      }
    }
  }
}

class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
             int n_step, int sep, int &start,
             double h, double log_sep,
             const Vec_2<double>& gl_l,
             size_t n_par_gl,
             const std::string& open_flag,
             MPI_Comm group_comm);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr,
                         float* pos, uint32_t* type_id, float* v);


  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr,
                         float* pos, float* v);

  uint64_t get_time_step();

  int reset_start_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TFloat, typename TInt>
  void read(int i_frame, TFloat* x, TFloat* y, TFloat* theta, TInt* type_id);

  template <typename TFloat, typename TInt>
  void read_last_frame(TFloat* x, TFloat* y, TFloat* theta, TInt* type_id);

  void ini(const std::string& filename, const std::string& open_flag, const Vec_2<double>& gl_l);

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

template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos, float* v) {
  int n_par = p_arr.size();
  int* n_par_arr = new int[tot_proc_] {};
  int* displs = new int[tot_proc_] {};
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
  float* my_v = new float[n_par];


  double half_Lx = gl_l_.x * 0.5;
  double half_Ly = gl_l_.y * 0.5;
  for (int j = 0; j < n_par; j++) {
    int j3 = j * 3;
    my_pos[j3] = p_arr[j].pos.x - half_Lx;
    my_pos[j3 + 1] = p_arr[j].pos.y - half_Ly;
    my_pos[j3 + 2] = p_arr[j].get_theta();
    my_v[j] = p_arr[j].cal_v();
  }


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
  delete[] my_v;
}


/********************************************************************
 * Caution: the particle speed output by this function is not correct,
 * because the local densities are lost when one particle moves from
 * one core to another core.
*********************************************************************/
template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    float* pos = nullptr;
#ifdef MULT_SPECIES
    uint32_t* type_id = nullptr;
#endif
    float* v = nullptr;
    
    if (my_rank_ == 0) {
      pos = new float[gl_np_ * 3];
#ifdef MULT_SPECIES
      type_id = new uint32_t[gl_np_];
#endif
      v = new float[gl_np_];
    }
#ifdef MULT_SPECIES
    get_data_from_par(p_arr, pos, type_id, v);
#else
    get_data_from_par(p_arr, pos, v);
#endif

    uint32_t n_par = gl_np_;
    if (my_rank_ == 0) {
      //uint64_t step = get_time_step();
      uint64_t step = start_ + i_step;
      //std::cout << "dump frame " << step << std::endl;
      gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
      gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
      gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
#ifdef MULT_SPECIES
      gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);
#endif
      gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_par, 1, 0, v);
      gsd_end_frame(handle_);
    }
    delete[] pos;
    delete[] v;
#ifdef MULT_SPECIES
    delete[] type_id;
#endif
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
#ifdef MULT_SPECIES
    uint32_t* type_id0 = new uint32_t[n_par];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/typeid");
    gsd_read_chunk(handle_, type_id0, chunk);
#endif
    double half_Lx = gl_l_.x * 0.5;
    double half_Ly = gl_l_.y * 0.5;
    for (int i = 0; i < n_par; i++) {
      x[i] = pos[i * 3 + 0] + half_Lx;
      y[i] = pos[i * 3 + 1] + half_Ly;
      theta[i] = pos[i * 3 + 2];
#ifdef MULT_SPECIES
      type_id[i] = type_id0[i];
#endif
      tangle_1D(x[i], gl_l_.x);
      tangle_1D(y[i], gl_l_.y);
    }

    delete[] pos;
#ifdef MULT_SPECIES
    delete[] type_id0;
#endif
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


/**
 * @brief Exporter for coarse-grained density and momentum field in 2D.
 *
 * Output coarse-grained density fields for A and B particles over boxes with linear
 * size "box_size" every "sep" time steps.
 *
 * Caution: the number of grids in both x and y directions must be multiple
 * of linear size of boxes for coarse graining.
 */
class FeildExporter : public ExporterBase {
public:
  template <typename TGrid, typename TDomain>
  FeildExporter(const std::string& outfile, int start, int n_step, int sep,
    const TGrid& grid, const TDomain& dm, int box_size = 4);

  template <typename TPar>
  void coarse_grain(const std::vector<TPar>& p_arr, unsigned short* rho_A, unsigned short* rho_B);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  void write_data(const unsigned short* rho_A, const unsigned short* rho_B);

  ~FeildExporter() { delete[] offset_; MPI_File_close(&fh_); }

protected:
  Vec_2<double> origin_;
  Vec_2<int> gl_n_;
  Vec_2<int> n_;
  Vec_2<int> o_;
  Vec_2<double> inv_lc_;

  MPI_File fh_{};
  MPI_Offset frame_size_;
  int n_grids_;
  MPI_Offset* offset_;
  int idx_frame_;
};

template <typename TGrid, typename TDomain>
FeildExporter::FeildExporter(const std::string& outfile, int start, int n_step, int sep,
  const TGrid& grid, const TDomain& dm, int box_size)
  : ExporterBase(n_step, sep, start, dm.comm()), origin_(dm.origin()) {
  gl_n_.x = grid.gl_n().x / box_size;
  gl_n_.y = grid.gl_n().y / box_size;
  n_.x = grid.n().x / box_size;
  n_.y = grid.n().y / box_size;
  o_.x = grid.origin().x / box_size;
  o_.y = grid.origin().y / box_size;
  inv_lc_.x = grid.lc().x / box_size;
  inv_lc_.y = grid.lc().y / box_size;
  frame_size_ = sizeof(unsigned short) * 2 * gl_n_.x * gl_n_.y;
  n_grids_ = n_.x * n_.y;

  offset_ = new MPI_Offset[n_.y * 2]{};

  for (int i_field = 0; i_field < 2; i_field++) {
    MPI_Offset start = sizeof(unsigned short) * gl_n_.x * gl_n_.y  * i_field;
    for (int row = 0; row < n_.y; row++) {
      int idx = row + i_field * n_.y;
      offset_[idx] = start + sizeof(unsigned short) * (o_.x + (o_.y + row) * gl_n_.x);
    }
  }

  MPI_File_open(dm.comm(), outfile.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh_);
  MPI_File_set_size(fh_, 0);
  idx_frame_ = 0;
}

template<typename TPar>
void FeildExporter::coarse_grain(const std::vector<TPar>& p_arr, unsigned short* rho_A, unsigned short* rho_B) {
  const auto end = p_arr.end();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    const TPar& p = *it;
    int idx = int((p.pos.x - origin_.x) * inv_lc_.x) + int((p.pos.y - origin_.y) * inv_lc_.y) * n_.x;
    if (p.type_id == 0) {
      rho_A[idx]++;
    } else {
      rho_B[idx]++;
    }
  }
}

template<typename TPar>
void FeildExporter::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    unsigned short *rho_A = new unsigned short[n_grids_]();
    unsigned short *rho_B = new unsigned short[n_grids_]();

    coarse_grain(p_arr, rho_A, rho_B);
    write_data(rho_A, rho_B);

    delete[] rho_A;
    delete[] rho_B;
  }
}


}

