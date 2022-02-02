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

  ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step) {
    set_lin_frame(start, n_step, sep);
  }

  void set_lin_frame(int start, int n_step, int sep);

  bool need_export(const int i_step);

protected:
  int n_step_;    // total steps to run
  int start_ = 0; // The first step 
private:
  std::vector<int> frames_arr_; // frames that need to export
  std::vector<int>::iterator frame_iter_;
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
  Snap_GSD_2(const std::string& filename, int start, int n_step, int sep,
             const Vec_2<double>& gl_l, const std::string& open_flag);

  ~Snap_GSD_2();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TPar>
  void read(int i_frame, std::vector<TPar>& p_arr);

  template <typename TPar>
  void read_last_frame(std::vector<TPar>& p_arr) {
    read(gsd_get_nframes(handle_) - 1, p_arr);
  }

private:
  gsd_handle* handle_ = nullptr;
  Vec_2<double> half_l_;
};

template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    uint32_t n_par = p_arr.size();
    float* pos = new float[n_par * 3];
    uint32_t *type_id = new uint32_t[n_par];

    for (size_t j = 0; j < n_par; j++) {
      size_t j3 = j * 3;
      pos[j3    ] = p_arr[j].pos.x - half_l_.x;
      pos[j3 + 1] = p_arr[j].pos.y - half_l_.y;
      pos[j3 + 2] = p_arr[j].get_theta();
      type_id[j] = p_arr[j].type_id;
    }

    size_t n_frame = gsd_get_nframes(handle_);
    size_t i_frame;
    if (n_frame == 0) {
      i_frame = 0;
    } else {
      const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
      gsd_read_chunk(handle_, &i_frame, chunk);
      i_frame++;
    }
    //std::cout << "dump frame " << i_frame << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &i_frame);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);
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
  std::cout << "find " << n_par << " particles" << std::endl;

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

}

