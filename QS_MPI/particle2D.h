#pragma once
#include "vect.h"
#include <iostream>
#include <cstring>
#define M_PI 3.14159265358979323846

template <typename T>
void tangle_1D(T& x, T L) {
  if (x < 0) {
    x += L;
  } else if (x >= L) {
    x -= L;
  }
}

template <typename TPar, typename TDomain>
bool check_overlap_2(const TPar& p_new, const std::vector<TPar>& p_arr,
                     double sigma_square, const TDomain& dm) {
  bool flag_overlap = false;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    Vec_2<double> r12(p_new.pos - (*it).pos);
    dm.untangle(r12);
    if (r12.square() < sigma_square) {
      flag_overlap = true;
      break;
    }
  }
  return flag_overlap;
}

template <typename TPar, typename TDomain, typename TRan>
void create_rand_2(std::vector<TPar>& p_arr, int n, TRan& myran,
                   const TDomain& dm, double sigma = 1.) {
  const double sigma_square = sigma * sigma;
  p_arr.reserve(n);
  for (int i = 0; i < n; i++) {
    int count = 0;
    int max_count = 1000;
    while (count < max_count) {
      TPar p_new(myran, dm.l(), dm.origin());
      bool flag_overlap = check_overlap_2(p_new, p_arr, sigma_square, dm);
      if (!flag_overlap) {
        p_arr.push_back(p_new);
        break;
      }
      count++;
    }
    if (count >= max_count) {
      std::cout << "count = " << count << std::endl;
      exit(1);
    }
  }

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      Vec_2<double> r12_vec = p_arr[j].pos - p_arr[i].pos;
      dm.untangle(r12_vec);
      if (r12_vec.square() < sigma_square) {
        std::cout << r12_vec << std::endl;
        exit(1);
      }
    }
  }
}

class BP_2 {
public:
  BP_2() : pos(), f() {}
  BP_2(const Vec_2<double>& pos0) : pos(pos0), f() {}
  BP_2(const Vec_2<double>& pos0, const Vec_2<double>& f0)
    : pos(pos0), f(f0) {}
  template <typename TRan>
  BP_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o)
    : pos(o.x + myran.doub() * l.x, o.y + myran.doub() * l.y), f() {}

  Vec_2<double> pos;
  Vec_2<double> f;
};

class BP_u_2: public BP_2 {
public:
  BP_u_2() :BP_2(), u() {}
  BP_u_2(const Vec_2<double>& pos0, const Vec_2<double>& u0) : BP_2(pos0), u(u0) {}
  BP_u_2(const Vec_2<double>& pos0, const Vec_2<double>& u0, const Vec_2<double>& f0) 
    : BP_2(pos0, f0), u(u0) {}
  template <typename TRan>
  BP_u_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) : BP_2(myran, l, o) {
    double theta = myran.doub() * M_PI * 2.;
    u.x = cos(theta);
    u.y = sin(theta);
  }

  const Vec_2<double>& update_ori(double d_theta) { 
    u.rotate(d_theta); return u; }

  Vec_2<double> u;
#ifdef HAS_CELL_INDEX
  int ic = 0;
#endif
};

class BP_u_tau_2 : public BP_u_2 {
public:
  BP_u_tau_2() :BP_u_2(), tau(0.) {}
  BP_u_tau_2(const Vec_2<double>& pos0, const Vec_2<double>& u0) 
    : BP_u_2(pos0, u0), tau(0.)  {}
  BP_u_tau_2(const Vec_2<double>& pos0, const Vec_2<double>& u0, const Vec_2<double>& f0)
    : BP_u_2(pos0, u0, f0), tau(0.) {}
  template <typename TRan>
  BP_u_tau_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) 
    : BP_u_2(myran, l, o), tau(0.) {}
  double tau;
};

class BP_theta_2 : public BP_2 {
public:
  BP_theta_2() : BP_2(), theta(0.) {}
  BP_theta_2(const Vec_2<double>& pos0, double theta0) : BP_2(pos0), theta(theta0) {}
  BP_theta_2(const Vec_2<double>& pos0, double theta0, const Vec_2<double>& f0)
    : BP_2(pos0, f0), theta(theta0) {}
  template <typename TRan>
  BP_theta_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o)
    : BP_2(myran, l, o), theta(myran.doub()* M_PI * 2.) {}

  const Vec_2<double> update_ori(double d_theta) {
    theta += d_theta;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return Vec_2<double>(c, s);
  }
  double theta;
};


class BP_theta_tau_2 : public BP_theta_2 {
public:
  BP_theta_tau_2() : BP_theta_2(), tau(0.) {}
  BP_theta_tau_2(const Vec_2<double>& pos0, double theta0) : BP_theta_2(pos0, theta0), tau(0.) {}
  BP_theta_tau_2(const Vec_2<double>& pos0, double theta0, const Vec_2<double>& f0)
    : BP_theta_2(pos0, theta0, f0), tau(0.) {}

  template <typename TRan>
  BP_theta_tau_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o)
    : BP_theta_2(myran, l, o), tau(0.) {}

  double tau;
};

class QSP_2 {
public:
  QSP_2() : u(), pos(), type_id(0) {}
  QSP_2(const Vec_2<double>& pos0) : u(), pos(), type_id(0) {}
  QSP_2(const Vec_2<double>& pos0, const Vec_2<double>& u0)
    :u(u0), pos(pos0), type_id(0) {}
  template <typename TRan>
  QSP_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  double get_theta() const { return atan2(u.y, u.x); }

  const Vec_2<double>& update_ori(double d_theta) { u.rotate(d_theta); return u; }

  template <typename TInt>
  void copy_pos_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_pos_from(const char* src, TInt& idx);

  template <typename TInt>
  void copy_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_from(const char* dest, TInt& idx);

  double cal_v() const {
    return v0_
      * (1 + kappa_ * tanh(eta_[type_id][0] * (rho_local[0] - rho0_)))
      * (1 + kappa_ * tanh(eta_[type_id][1] * (rho_local[1] - rho0_)));
  };

  static void set_QS_params(double v0, double kappa, double rho0,
                            double eta_AA, double eta_AB,
                            double eta_BA, double eta_BB);

  Vec_2<double> u;
  Vec_2<double> pos;
  uint8_t type_id;
  float rho_local[2]{};

private:
  static double v0_;
  static double kappa_;
  static double rho0_;
  static double eta_[2][2];   // eta / (rho0 * kappa)
};

template <typename TRan>
QSP_2::QSP_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) {
  double theta = myran.doub() * M_PI * 2.;
  u.x = cos(theta);
  u.y = sin(theta);
  pos.x = o.x + myran.doub() * l.x;
  pos.y = o.y + myran.doub() * l.y;
  type_id = 0;
}

template <typename TInt>
void QSP_2::copy_pos_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &pos.x, 17);
  idx += 17;
}

template <typename TInt>
void QSP_2::copy_pos_from(const char* src, TInt& idx) {
  std::memcpy(&pos.x, src + idx, 17);
  idx += 17;
}

template <typename TInt>
void QSP_2::copy_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &u.x, 33);
  idx += 33;
}

template <typename TInt>
void QSP_2::copy_from(const char* src, TInt& idx) {
  std::memcpy(&u.x, src + idx, 33);
  idx += 33;
}