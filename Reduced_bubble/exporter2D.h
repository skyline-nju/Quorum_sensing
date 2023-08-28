#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "vect.h"

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

    ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step), sep_(sep) {
      set_lin_frame(start, n_step, sep);
    }

    void set_lin_frame(int start, int n_step, int sep);

    bool need_export(const int i_step);

  protected:
    int n_step_;    // total steps to run
    int sep_;
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

  /**
 * @brief Output order parameters.
 *
 * If set use_sub_boxes=1, then calculate local orderparameters over boxes with
 * varies lienar size, otherwise just calcalute globa order parameters.
 */
  class OrderParaExporter_2 : public ExporterBase {
  public:
    OrderParaExporter_2(const std::string& outfile,
      int start, int n_step, int sep);

    void dump(int i_step, int n, double r) {
      if (need_export(i_step)) {
        fout_ << std::fixed << std::setw(16) << std::setprecision(10)
          << i_step << "\t" << n << "\t" << r << std::endl;
      }
    }

  private:
    std::ofstream fout_;
  };


}