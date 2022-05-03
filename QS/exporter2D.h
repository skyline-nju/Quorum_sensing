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

namespace exporter {

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
  ExporterBase() : n_step_(0) {}

  ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step), sep_(sep) {}

  bool need_export(const int i_step) const {return i_step % sep_ == 0;}

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step
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
#ifdef USE_MPI
  LogExporter(const std::string& outfile, int start, int n_step, int sep,
    int np, MPI_Comm group_comm);
#else
  LogExporter(const std::string& outfile, int start, int n_step, int sep,
    int np);
#endif

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
#ifdef USE_MPI
  MPI_Comm comm_;
#endif
  int step_count_ = 0;
};


class OrderParaExporter : public ExporterBase {
public:
  OrderParaExporter(const std::string& outfile, int n_step, int sep, int start,
                    int flush_dt, int npar_gl);

  ~OrderParaExporter();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, const float* v_arr);

private:
  std::ofstream fout_;
  int flush_dt_;
  int npar_gl_;
};

template<typename TPar>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& p_arr, const float* v_arr){
  if (need_export(i_step)) {
    int n = p_arr.size();
    double sum_rho[2]{};
    double sum_rho_square[2]{};
    double sum_rho_cross = 0;
    double sum_v = 0;
    double sum_v_square = 0;

    for (const auto& p : p_arr) {
      sum_rho[0] += p.rho_local[0];
      sum_rho_square[0] += p.rho_local[0] * p.rho_local[0];
      sum_rho[1] += p.rho_local[1];
      sum_rho_square[1] += p.rho_local[1] * p.rho_local[1];
      sum_rho_cross += p.rho_local[0] * p.rho_local[1];
    }

    for (int i = 0; i < n; i++) {
      sum_v += v_arr[i];
      sum_v_square += v_arr[i] * v_arr[i];
    }

    double sum[7] = {sum_rho[0], sum_rho[1],
                     sum_rho_square[0], sum_rho_square[1], sum_rho_cross,
                     sum_v, sum_v_square};
    
    double D_rho_A = sum[2] / npar_gl_ - (sum[0] / npar_gl_) * (sum[0] / npar_gl_);
    double D_rho_B = sum[3] / npar_gl_ - (sum[1] / npar_gl_) * (sum[1] / npar_gl_);
    double D_rho_c = sum[4] / npar_gl_ - (sum[0] / npar_gl_) * (sum[1] / npar_gl_);
    double D_v = sum[6] / npar_gl_ - (sum[5] / npar_gl_) * (sum[5] / npar_gl_);
    fout_ << std::fixed << std::setw(10) << std::setprecision(8) << start_ + i_step 
      << "\t" << D_v << "\t" << D_rho_A << "\t" << D_rho_B << "\t" << D_rho_c << "\n";
    if (i_step % flush_dt_ == 0) {
      fout_ << std::flush;
    }
  }
}

class XyzExporter_2 : public ExporterBase {
public:
  explicit XyzExporter_2(const std::string outfile, int start, int n_step, int sep,
    const Vec_2<double>& gl_l)
    : ExporterBase(start, n_step, sep), fout_(outfile), gl_l_(gl_l), sep_(sep) {}

  template <typename TPar>
  void dump_pos(int i_step, const std::vector<TPar>& par_arr);

  template <typename TPar>
  void dump_doub_pos(int i_step, const std::vector<TPar>& par_arr);

  template <typename TPar>
  void dump_multispecies(int i_step, const std::vector<TPar>& par_arr);

private:
  std::ofstream fout_;
  Vec_2<double> gl_l_;
  int sep_;
};

template <typename TPar>
void XyzExporter_2::dump_pos(int i_step, const std::vector<TPar>& par_arr) {
  //if (need_export(i_step)) {
  if (i_step % sep_ == 0) {
    int n_par = par_arr.size();
    fout_ << n_par << "\n";
    // comment line
    fout_ << "Lattice=\"" << gl_l_.x << " 0 0 0 " << gl_l_.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    for (int j = 0; j < n_par; j++) {
      fout_ << "\n" << "N\t"
        << par_arr[j].pos.x << "\t" << par_arr[j].pos.y;
    }
    fout_ << std::endl;
  }
}

template<typename TPar>
void XyzExporter_2::dump_doub_pos(int i_step, const std::vector<TPar>& par_arr) {
  if (i_step % sep_ == 0) {
    int n_par = par_arr.size();
    fout_ << n_par * 2 << "\n";
    // comment line
    fout_ << "Lattice=\"" << gl_l_.x << " 0 0 0 " << gl_l_.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    for (int j = 0; j < n_par; j++) {
      fout_ << "\n" << "N\t"
        << par_arr[j].pos.x << "\t" << par_arr[j].pos.y;
      double dx = 0.01 * par_arr[j].u.x;
      double dy = 0.01 * par_arr[j].u.y;
      fout_ << "\n" << "O\t"
        << par_arr[j].pos.x + dx << "\t" << par_arr[j].pos.y + dy;
    }
    fout_ << std::endl;
  }
}

template <typename TPar>
void XyzExporter_2::dump_multispecies(int i_step, const std::vector<TPar>& par_arr) {
  if (i_step % sep_ == 0) {
    int n_par = par_arr.size();
    fout_ << n_par << "\n";
    // comment line
    fout_ << "Lattice=\"" << gl_l_.x << " 0 0 0 " << gl_l_.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    for (int j = 0; j < n_par; j++) {
      fout_ << "\n";
      if (par_arr[j].type_id == 0) {
        fout_ << "N\t";
      } else {
        fout_ << "O\t";
      }
      fout_ << par_arr[j].pos.x << "\t" << par_arr[j].pos.y;
    }
    fout_ << std::endl;
  }
}

class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename, int& start, int n_step, int sep,
             const Vec_2<double>& gl_l, const std::string& open_flag);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos, uint32_t* type_id);

  uint64_t get_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, float* v);

  template <typename TPar>
  void read(int i_frame, std::vector<TPar>& p_arr);

  template <typename TPar>
  void read_last_frame(std::vector<TPar>& p_arr);

  int reset_start_time_step();

private:
  gsd_handle* handle_ = nullptr;
  Vec_2<double> half_l_;
};


template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos, uint32_t* type_id) {
  size_t n_par = p_arr.size();
  for (size_t j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    pos[j3    ] = p_arr[j].pos.x - half_l_.x;
    pos[j3 + 1] = p_arr[j].pos.y - half_l_.y;
    pos[j3 + 2] = p_arr[j].get_theta();
    type_id[j] = p_arr[j].type_id;
  }
}


template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
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
    gsd_end_frame(handle_);
    delete[] pos;
    delete []type_id;
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
template <typename TPar>
void Snap_GSD_2::read(int i_frame, std::vector<TPar>& p_arr) {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  std::cout << "frame " << i_frame  <<": find " << n_par << " particles" << std::endl;

  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);
  uint32_t* type_id = new uint32_t[n_par];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/typeid");
  gsd_read_chunk(handle_, type_id, chunk);

  p_arr.reserve(n_par);

  double Lx = half_l_.x * 2;
  double Ly = half_l_.y * 2;

  for (int j = 0; j < n_par; j++) {
    TPar p;
    size_t j3 = j * 3;
    double x = pos[j3] + half_l_.x;
    double y = pos[j3 + 1] + half_l_.y;
    double theta = pos[j3 + 2];
    
    if (x < 0) {
      x += Lx;
    } else if (x >= Lx) {
      x -= Lx;
    }
    if (y < 0) {
      y += Ly;
    } else if (y >= Ly) {
      y -= Ly;
    }

    p.pos = Vec_2<double>(x, y);
    p.u = Vec_2<double>(cos(theta), sin(theta));
    p.type_id = type_id[j];
    p_arr.push_back(p);
  }

  delete[] pos;
  delete[] type_id;
}


template <typename TPar>
void Snap_GSD_2::read_last_frame(std::vector<TPar>& p_arr) {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes-1, p_arr);
  }
}

}

