#include "io2D.h"

io::ExporterBase::ExporterBase(int n_step, int sep, int start,
                               MPI_Comm group_comm)
  : n_step_(n_step), sep_(sep), start_(start), comm_(group_comm) {
  MPI_Comm_rank(comm_, &my_rank_);
  MPI_Comm_size(comm_, &tot_proc_);
}

void io::ExporterBase::set_log_scale_frames(double h, double log_sep) {
  log_sep_ = log_sep;
  int frame_size = int(20 / log_sep);
  frames_.reserve(frame_size);
  std::vector<double> log10_t;
  log10_t.reserve(frame_size);
  log10_t.push_back(1);
  for (int i = 1; i < frame_size; i++) {
    log10_t.push_back(log10_t[i-1]+log_sep);
  }
  for (int i = 0; i < frame_size; i++) {
    double t = pow(10, log10_t[i]);
    int i_step = round(t/h);
    frames_.push_back(i_step);
  }
  int cur_step = start_ + 1;
  if (cur_step <= frames_[0]) {
    cur_frame_ = 0;
  } else {
    for (int i = 1; i < frames_.size(); i++) {
      if (cur_step <= frames_[i] && cur_step > frames_[i-1]) {
        cur_frame_ = i;
        break;
      }
    }
  }
}

io::LogExporter::LogExporter(const std::string& outfile, 
                             int n_step, int sep, int start,
                             int np_gl, MPI_Comm group_comm)
  : ExporterBase(n_step, sep, start, group_comm), n_par_(np_gl){
  if (my_rank_ == 0) {
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
#ifdef WITH_SELF
    fout << "With particle self for QS" << "\n";
#else
    fout << "Without particle self for QS" << "\n";
#endif
  }
  MPI_Barrier(group_comm);
}


io::LogExporter::~LogExporter() {
  if (my_rank_ == 0) {
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
  }
}

void io::LogExporter::record(int i_step) {
  if (need_export(i_step + start_) && my_rank_ == 0) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = floor(dt / 3600);
    const auto min = floor((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout << i_step + start_ << "\t" << hour << ":" << min << ":" << sec << std::endl;
  }
  step_count_++;
}

io::MassCenterExporter::MassCenterExporter(const std::string& outfile,
  int n_step, int sep, int start, const Vec_2<double>& gl_l, int flush_sep, MPI_Comm group_comm)
  : ExporterBase(n_step, sep, start, group_comm), Lx_(gl_l.x), Ly_(gl_l.y), flush_sep_(flush_sep) {
  if (my_rank_ == 0) {
    fout_.open(outfile, std::ios::app);
  }
  MPI_Barrier(group_comm);
}

io::MassCenterExporter::~MassCenterExporter() {
  if (my_rank_ == 0) {
    fout_.close();
  }
}


io::Snap_GSD_2::Snap_GSD_2(const std::string& filename,
                           int n_step, int sep, int &start,
                           double h, double log_sep,
                           const Vec_2<double>& gl_l,
                           size_t n_par_gl,
                           const std::string& open_flag,
                           MPI_Comm group_comm)
  : ExporterBase(n_step, sep, start, group_comm), gl_l_(gl_l), gl_np_(n_par_gl) {
  handle_ = new gsd_handle;
  ini(filename, open_flag, gl_l);
  start = reset_start_time_step();
  if (log_sep > 0) {
    set_log_scale_frames(h, log_sep);
  }
  MPI_Barrier(group_comm);
}

void io::Snap_GSD_2::ini(const std::string& filename, const std::string& open_flag, const Vec_2<double>& gl_l) {
  unsigned int version = gsd_make_version(1, 4);
  if (my_rank_ == 0) {
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
      std::cout << "Wrong open flag, which must be one of 'rand' or 'resume'!" << std::endl;
      exit(1);
    }
  }
}

io::Snap_GSD_2::~Snap_GSD_2() {
  if (my_rank_ == 0) {
    gsd_close(handle_);
  }
  delete handle_;
}


uint64_t io::Snap_GSD_2::get_time_step() {
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

int io::Snap_GSD_2::reset_start_time_step() {
  if (my_rank_ == 0) {
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
        std::cout << "Warning, failed to read the time step of the last frame" << std::endl;
        start_ = 0;
      }
    }
  }
  MPI_Bcast(&start_, 1, MPI_INT, 0, comm_);
  return start_;
}

void io::FeildExporter::write_data(const unsigned short* rho_A,
                                   const unsigned short* rho_B) {
  const MPI_Offset frame_start = idx_frame_ * frame_size_;
  MPI_Offset offset1, offset2;
  int pos = 0;
  for (int row = 0; row < n_.y; row++) {
    offset1 = frame_start + offset_[row];
    offset2 = frame_start + offset_[row + n_.y];
    MPI_File_write_at(fh_, offset1, &rho_A[pos], n_.x, MPI_UINT16_T, MPI_STATUSES_IGNORE);
    MPI_File_write_at(fh_, offset2, &rho_B[pos], n_.x, MPI_UINT16_T, MPI_STATUSES_IGNORE);
    pos += n_.x;
  }
  idx_frame_++;
}