#pragma once
#include "vect.h"
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class SpringForce_2 {
public:
  SpringForce_2() = default;
  explicit SpringForce_2(double sigma, double k) 
    : spring_const_(k), sigma_(sigma), r_cut_square_(sigma* sigma) {}

  void set_spring_const(double spring_const) { spring_const_ = spring_const; }

  void set_sigma(double sigma) { sigma_ = sigma; r_cut_square_ = sigma * sigma; }

  void eval(const Vec_2<double> &r12_vec, Vec_2<double>& f12_vec) const {
    const double r12_square = r12_vec.square();
    if (r12_square < r_cut_square_) {
      const double r12 = std::sqrt(r12_square);
      const double f = (sigma_ - r12) * spring_const_;
      const double tmp = f / r12;
      f12_vec.x = tmp * r12_vec.x;
      f12_vec.y = tmp * r12_vec.y;
    }
  }
private:
  double spring_const_;
  double sigma_;
  double r_cut_square_;
};

class WCAForce_2 {
public:
  WCAForce_2() = default;
  WCAForce_2(double eps, double sigma = 1) {
    eps24_ = eps * 24;
    double r_cut = std::pow(2.0, 1. / 6) * sigma;
    r_cut_square_ = r_cut * r_cut;
    sigma_square_ = sigma * sigma;
  }

  void set_sigma(double sigma) {
    double r_cut = std::pow(2.0, 1. / 6) * sigma;
    r_cut_square_ = r_cut * r_cut;
    sigma_square_ = sigma * sigma;
  }

  double get_r_cut() { return std::sqrt(r_cut_square_); }

  void eval(const Vec_2<double>& r12_vec, Vec_2<double>& f12_vec) const {
    const double r12_square = r12_vec.square();
    if (r12_square < r_cut_square_) {
      double r_2 = sigma_square_ / r12_square;
      double r_6 = r_2 * r_2 * r_2;
      double tmp = eps24_ * (2 * r_6 * r_6 - r_6) * r_2;
      f12_vec.x = tmp * r12_vec.x;
      f12_vec.y = tmp * r12_vec.y;
    }
  }
private:
  double eps24_;
  double r_cut_square_;
  double sigma_square_;
};

class LinearDensityKernal {
public:
  LinearDensityKernal(double r_cut) : r_cut_(r_cut), r_cut_square_(r_cut* r_cut) {}

  template <typename TPar>
  void accum_density(double r12, TPar& p1, TPar& p2) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

  template <typename TPar>
  void normalize(std::vector<TPar>& p_arr) const;

  template <typename TPar>
  void reset_local_density(std::vector<TPar>& p_arr) const;

private:
  double r_cut_;
  double r_cut_square_;
};

template<typename TPar>
void LinearDensityKernal::accum_density(double r12, TPar& p1, TPar& p2) const {
  double density_weighted = r_cut_ - r12;
  p1.rho_local += density_weighted;
  p2.rho_local += density_weighted;

  //!TODO need divided by a nomalized factor
}

template <typename TPar>
void LinearDensityKernal::operator ()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  const double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
}

template <typename TPar, typename BdyCondi>
void LinearDensityKernal::operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  const double r12_square = r12.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
}

template <typename TPar>
void LinearDensityKernal::normalize(std::vector<TPar>& p_arr) const {
  double norm = 3. / M_PI;
  for (auto& p : p_arr) {
    //p.rho_local[p.type_id] += 1;
    p.rho_local *= norm;

  }
}

template <typename TPar>
void LinearDensityKernal::reset_local_density(std::vector<TPar>& p_arr) const {
  for (auto& p : p_arr) {
    p.rho_local = 0.;
  }
}
