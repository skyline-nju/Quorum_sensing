#include "exporter2D.h"

#ifdef USE_MPI
exporter::LogExporter::LogExporter(const std::string& outfile, 
                                   int start, int n_step, int sep, 
                                   int np, MPI_Comm group_comm)
  : ExporterBase(start, n_step, sep), n_par_(np), comm_(group_comm) {
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  if (my_rank == 0) {
    fout.open(outfile);
    t_start_ = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(t_start_);
    char str[100];
    tm now_time;
#ifdef _MSC_VER

    localtime_s(&now_time, &start_time);
#else
    localtime_r(&start_time, &now_time);
#endif
    std::strftime(str, 100, "%c", &now_time);
    fout << "Started simulation at " << str << "\n";
  }
}
#else
exporter::LogExporter::LogExporter(const std::string& outfile,
                                   int start, int n_step, int sep, int np)
  : ExporterBase(start, n_step, sep), n_par_(np) {
  fout.open(outfile, std::ios::app);
  t_start_ = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(t_start_);
  char str[100];
  tm now_time;
#ifdef _MSC_VER

  localtime_s(&now_time, &start_time);
#else
  localtime_r(&start_time, &now_time);
#endif
  std::strftime(str, 100, "%c", &now_time);
  fout << "Started simulation at " << str << "\n";
}
#endif

exporter::LogExporter::~LogExporter() {
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  if (my_rank == 0) {
#endif
    const auto t_now = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::to_time_t(t_now);
    char str[100];
    tm now_time;
#ifdef _MSC_VER
    localtime_s(&now_time, &end_time);
#else
    localtime_r(&end_time, &now_time);
#endif
    std::strftime(str, 100, "%c", &now_time);
    fout << "Finished simulation at " << str << "\n";
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    fout << "speed=" << std::scientific << step_count_ * double(n_par_) / elapsed_seconds.count()
      << " particle time step per second per core\n";
    fout.close();
#ifdef USE_MPI
}
#endif
}

void exporter::LogExporter::record(int i_step) {
  bool flag;
#ifdef USE_MPI
  int my_rank;
  MPI_Comm_rank(comm_, &my_rank);
  flag = my_rank == 0 && need_export(i_step);
#else
  flag = need_export(i_step);
#endif
  if (flag) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = int(dt / 3600);
    const auto min = int((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout << i_step << "\t" << hour << ":" << min << ":" << sec << std::endl;
  }
  step_count_++;
}

exporter::OrderParaExporter::OrderParaExporter(const std::string& outfile,
                                               int n_step, int sep, int start,
                                               int flush_dt, int npar_gl)
  : ExporterBase(start, n_step, sep), flush_dt_(flush_dt), npar_gl_(npar_gl) {
  fout_.open(outfile, std::ios::app);
}

exporter::OrderParaExporter::~OrderParaExporter() {
  fout_.close();
}

exporter::Snap_GSD_2::Snap_GSD_2(const std::string& filename,
                                 int& start, int n_step, int sep,
                                 const Vec_2<double>& gl_l,
                                 const std::string& open_flag)
  : ExporterBase(start, n_step, sep) {
  unsigned int version = gsd_make_version(1, 4);
  handle_ = new gsd_handle;
  if (open_flag == "rand") {
    int flag = gsd_create(filename.c_str(), "cpp", "hoomd", version);
    if (flag != 0) {
      std::cout << "Error when create " << filename << "; state=" << flag << std::endl;
      exit(1);
    }
    flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    }

    float box[6] = { gl_l.x, gl_l.y, 1, 0, 0, 0 };
    gsd_write_chunk(handle_, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, box);
    
    char types[] = {'A', 'B'};
    gsd_write_chunk(handle_, "particles/types", GSD_TYPE_INT8, 2, 1, 0, types);
  } else if (open_flag == "resume") {
    int flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    } else {
      std::cout << "open " << filename << std::endl;
    } 
  } else {
    std::cout << "Wrong open flag, which must be one of 'new' or 'append'!" << std::endl;
    exit(1);
  }
  half_l_ = Vec_2<double>(gl_l.x / 2, gl_l.y / 2);
  start = reset_start_time_step();
}

exporter::Snap_GSD_2::~Snap_GSD_2() {
  gsd_close(handle_);
  delete handle_;
}


uint64_t exporter::Snap_GSD_2::get_time_step() {
  uint64_t step;
  size_t n_frame = gsd_get_nframes(handle_);
  if (n_frame == 0) {
    step = sep_;
  } else {
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
    if (chunk) {
      gsd_read_chunk(handle_, &step, chunk);
      step += sep_;
    } else {
      step = sep_;
    }
  }
  return step;
}

int exporter::Snap_GSD_2::reset_start_time_step() {
  size_t n_frame = gsd_get_nframes(handle_);
  if (n_frame == 0) {
    start_ = 0;
  } else {
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
    if (chunk) {
      uint64_t last_step;
      gsd_read_chunk(handle_, &last_step, chunk);
      start_ = last_step;
    } else {
      std::cout << "Error, failed to read the time step of the last frame" << std::endl;
      exit(1);
    }
  }
  return start_;
}