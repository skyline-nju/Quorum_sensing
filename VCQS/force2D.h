#pragma once
#include "vect.h"
#include <cmath>

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

// U(r12_vec, q1, q2) = C * exp(-lambda(r-1.))/r12**2 (q1 - q2) * r21_vec
class AmphiphilicWCA_2 {
public:
  AmphiphilicWCA_2() = default;
  AmphiphilicWCA_2(double eps, double lambda, double c, double r_cut)
    : C_(c), lambda_(lambda) {
    eps24_ = eps * 24;
    double r_cut_WCA = std::pow(2.0, 1. / 6);
    rcut_square_WCA_ = r_cut_WCA * r_cut_WCA;
    rcut_square_AN_ = r_cut * r_cut;
  }

  void eval(const Vec_2<double>& r12_vec, const Vec_2<double>& q1, const Vec_2<double>& q2,
    Vec_2<double>& f12_vec, double& tau1, double& tau2) {
    const double r12_square = r12_vec.square();
    if (r12_square < rcut_square_AN_) {
      double r = sqrt(r12_square);
      double r_2 = 1. / r12_square;
      double grad_V_pre = (lambda_ * r + 2.) * r_2;
      double V = C_ * exp(-lambda_ * (r - 1.)) * r_2;
      Vec_2<double> Vq1 = V * q1;
      Vec_2<double> Vq2 = V * q2;
      Vec_2<double> Vq1_m_q2 = Vq1 - Vq2;
      f12_vec = -(grad_V_pre * (Vq1_m_q2.dot(r12_vec))) * r12_vec + Vq1_m_q2;
      tau1 = Vq1.cross(r12_vec);
      tau2 = -Vq2.cross(r12_vec);
      if (r12_square < rcut_square_WCA_) {
        double r_6 = r_2 * r_2 * r_2;
        double tmp = eps24_ * (2 * r_6 * r_6 - r_6) * r_2;
        f12_vec.x += tmp * r12_vec.x;
        f12_vec.y += tmp * r12_vec.y;
      }
    }
  }

private:
  double C_;
  double lambda_;
  double eps24_;
  double rcut_square_WCA_;
  double rcut_square_AN_;
};

class AsymKernal {
public:
  AsymKernal(double r_cut, double epsilon=0)
    : r_cut_(r_cut), r_cut_square_(r_cut* r_cut), eps_(epsilon) {}

  template <typename TPar>
  void accum_density(double r12, TPar& p1, const TPar& p2) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

  template <typename TPar>
  void normalize(std::vector<TPar>& p_arr) const;

  template <typename TPar>
  void reset_local_density(std::vector<TPar>& p_arr) const;

protected:
  double r_cut_;
  double r_cut_square_;
  double eps_;
};

template<typename TPar>
void AsymKernal::accum_density(double r12, TPar& p1, const TPar& p2) const {
  double density_weighted = 1 - r12 / r_cut_;
  p1.rho_local[p2.type_id] += density_weighted;
  //p2.rho_local[p1.type_id] += density_weighted;

  //!TODO need divided by a nomalized factor
}

template <typename TPar>
void AsymKernal::operator ()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - (p1.pos + eps_ * p1.u);
  Vec_2<double> r21 = p1.pos - (p2.pos + eps_ * p2.u);
  const double r12_square = r12.square();
  const double r21_square = r21.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
  if (r21_square < r_cut_square_) {
    double r21_module = sqrt(r21_square);
    accum_density(r21_module, p2, p1);
  }
}

template <typename TPar, typename BdyCondi>
void AsymKernal::operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - (p1.pos + eps_ * p1.u);
  Vec_2<double> r21 = p1.pos - (p2.pos + eps_ * p2.u);
  bc.untangle(r12);
  bc.untangle(r21);
  const double r12_square = r12.square();
  const double r21_square = r21.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
  if (r21_square < r_cut_square_) {
    double r21_module = sqrt(r21_square);
    accum_density(r21_module, p2, p1);
  }
}

template <typename TPar>
void AsymKernal::normalize(std::vector<TPar>& p_arr) const {
  double norm = 3. / M_PI / r_cut_square_;
  for (auto& p : p_arr) {
    //p.rho_local[p.type_id] += 1;
    p.rho_local[0] *= norm;
    p.rho_local[1] *= norm;
  }
  //std::cout << p_arr[0].rho_local[0] << std::endl;
}

template <typename TPar>
void AsymKernal::reset_local_density(std::vector<TPar>& p_arr) const {
  for (auto& p : p_arr) {
    p.rho_local[0] = 0.;
    p.rho_local[1] = 0.;
  }
}


class VisionConeKernal {
public:
  VisionConeKernal(double r_cut, double visual_angle) : r_cut_(r_cut), r_cut_square_(r_cut* r_cut) {
    cross_dot_ = cos(visual_angle / 2 / 180 * M_PI);
    norm_ = 3 / M_PI / (visual_angle / 360);
  }

  template <typename TPar>
  void accum_density(double r12, TPar& p1, TPar& p2) const;

  template <typename TPar>
  void accum_density(const Vec_2<double> r12_vec, TPar& p1, TPar& p2) const;

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
  double cross_dot_;
  double norm_;

};

template<typename TPar>
void VisionConeKernal::accum_density(double r12, TPar& p1, TPar& p2) const {
  double density_weighted = r_cut_ - r12;
  p1.rho_local[p2.type_id] += density_weighted;
  p2.rho_local[p1.type_id] += density_weighted;

  //!TODO need divided by a nomalized factor
}

template<typename TPar>
void VisionConeKernal::accum_density(const Vec_2<double> r12_vec, TPar& p1, TPar& p2) const {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    double density_weighted = r_cut_ - r12_module;
    double tmp = cross_dot_ * r12_module;
    if (r12_vec.dot(p1.u) > tmp) {
      p1.rho_local[p2.type_id] += density_weighted;
    }
    if (r12_vec.dot(p2.u) < -tmp) {
      p2.rho_local[p1.type_id] += density_weighted;
    }

  }
}

template <typename TPar>
void VisionConeKernal::operator ()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  accum_density(r12, p1, p2);
  //const double r12_square = r12.square();
  //if (r12_square < r_cut_square_) {
  //  double r12_module = sqrt(r12_square);
  //  accum_density(r12_module, p1, p2);
  //}
}

template <typename TPar, typename BdyCondi>
void VisionConeKernal::operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  accum_density(r12, p1, p2);
  //const double r12_square = r12.square();
  //if (r12_square < r_cut_square_) {
  //  double r12_module = sqrt(r12_square);
  //  accum_density(r12_module, p1, p2);
  //}
}

template <typename TPar>
void VisionConeKernal::normalize(std::vector<TPar>& p_arr) const {
  //double norm = 3. / M_PI;
  for (auto& p : p_arr) {
    //p.rho_local[p.type_id] += 1;
    p.rho_local[0] *= norm_;
    p.rho_local[1] *= norm_;
  }
}

template <typename TPar>
void VisionConeKernal::reset_local_density(std::vector<TPar>& p_arr) const {
  for (auto& p : p_arr) {
    p.rho_local[0] = 0.;
    p.rho_local[1] = 0.;
  }
}

class NRTorque_2 {
public:
  NRTorque_2() = default;
  NRTorque_2(double eta, double r_cut) : Eta_(eta), r_cut_square_(r_cut* r_cut) {}

  template <typename TPar>
  void eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const;

  template<typename TPar>
  void operator()(TPar& p1, TPar& p2) const;

  template<typename TPar, typename BdyCondi>
  void operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const;
protected:
  double Eta_;
  double r_cut_square_;
};

template <typename TPar>
void NRTorque_2::eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r = sqrt(r12_square);
    Vec_2<double> r12_hat = r12_vec / r;
    p1.tau += Eta_ * p1.u.cross(r12_hat);
    p2.tau -= Eta_ * p2.u.cross(r12_hat);
    p1.n_neighbor += 1;
    p2.n_neighbor += 1;
  }
}

template<typename TPar>
void NRTorque_2::operator()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  eval(r12, p1, p2);
}

template<typename TPar, typename BdyCondi>
void NRTorque_2::operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  eval(r12, p1, p2);
}


class NRTorque_Rep_2 {
public:
  NRTorque_Rep_2() = default;
  NRTorque_Rep_2(double eta, double k, double r_cut) : Eta_(eta), k_(k),  r_cut_(r_cut), r_cut_square_(r_cut * r_cut) {}

  template <typename TPar>
  void eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const;

  template<typename TPar>
  void operator()(TPar& p1, TPar& p2) const;

  template<typename TPar, typename BdyCondi>
  void operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const;
protected:
  double Eta_;
  double k_;
  double r_cut_;
  double r_cut_square_;
};

template<typename TPar>
void NRTorque_Rep_2::eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r = sqrt(r12_square);
    double dr = r_cut_ - r;
    Vec_2<double> r12_hat = r12_vec / r;
    p1.tau += Eta_ * dr * p1.u.cross(r12_hat);
    p2.tau -= Eta_ * dr * p2.u.cross(r12_hat);
    Vec_2<double> f_rep = k_ * dr * r12_vec;
    p1.f -= f_rep;
    p2.f += f_rep;
    p1.n_neighbor += 1;
    p2.n_neighbor += 1;
  }
}

template<typename TPar>
void NRTorque_Rep_2::operator()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  eval(r12, p1, p2);
}

template<typename TPar, typename BdyCondi>
void NRTorque_Rep_2::operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  eval(r12, p1, p2);
}

class NRTorque_Rep_hardcore_2 : public NRTorque_Rep_2 {
public:
  NRTorque_Rep_hardcore_2() = default;
  NRTorque_Rep_hardcore_2(double eta, double k, double r_cut, double k2, double sigma)
    : NRTorque_Rep_2(eta, k, r_cut), k2_(k2), sigma_(sigma){}

  template <typename TPar>
  void eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const;

  template<typename TPar>
  void operator()(TPar& p1, TPar& p2) const {
    Vec_2<double> r12 = p2.pos - p1.pos;
    eval(r12, p1, p2);
  }

  template<typename TPar, typename BdyCondi>
  void operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
    Vec_2<double> r12 = p2.pos - p1.pos;
    bc.untangle(r12);
    eval(r12, p1, p2);
  }

  double k2_;
  double sigma_;
};

template<typename TPar>
void NRTorque_Rep_hardcore_2::eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r = sqrt(r12_square);
    double dr = r_cut_ - r;
    Vec_2<double> r12_hat = r12_vec / r;
    p1.tau += Eta_ * dr * p1.u.cross(r12_hat);
    p2.tau -= Eta_ * dr * p2.u.cross(r12_hat);
    Vec_2<double> f_rep = k_ * dr * r12_vec;
    if (r < sigma_) {
      f_rep += k2_ * (sigma_ - r) * r12_vec;
    }
    p1.f -= f_rep;
    p2.f += f_rep;
    p1.n_neighbor += 1;
    p2.n_neighbor += 1;
  }
}

class NRTorque_QS_2 {
public:
  NRTorque_QS_2(double r_cut, double torque_strength, double rho0)
    : r_cut_(r_cut), r_cut_square_(r_cut* r_cut), torque_strength_(torque_strength/rho0) {}

  template <typename TPar>
  void eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2)  const;

  template<typename TPar>
  void operator()(TPar& p1, TPar& p2) const {
    Vec_2<double> r12 = p2.pos - p1.pos;
    eval(r12, p1, p2);
  }

  template<typename TPar, typename BdyCondi>
  void operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
    Vec_2<double> r12 = p2.pos - p1.pos;
    bc.untangle(r12);
    eval(r12, p1, p2);
  }

  template <typename TPar>
  void normalize(std::vector<TPar>& p_arr) const;

  template <typename TPar>
  void reset_local_density(std::vector<TPar>& p_arr) const {
    for (auto& p : p_arr) {
      p.rho_local = 0.;
    }
  }
private:
  double r_cut_;
  double r_cut_square_;
  double torque_strength_;
};

template <typename TPar>
void NRTorque_QS_2::eval(const Vec_2<double>& r12_vec, TPar& p1, TPar& p2) const {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r = sqrt(r12_square);
    Vec_2<double> r12_hat = r12_vec / r;
    double dr = r_cut_ - r;
    p1.tau += dr * p1.u.cross(r12_hat);
    p2.tau -= dr * p2.u.cross(r12_hat);
    p1.rho_local += dr;
    p2.rho_local += dr;
  }
}

template<typename TPar>
void NRTorque_QS_2::normalize(std::vector<TPar>& p_arr) const {
  double norm = 3. / M_PI;
  for (auto& p : p_arr) {
    p.rho_local += 1;
    p.rho_local *= norm;
    p.tau *= torque_strength_;
  }
}