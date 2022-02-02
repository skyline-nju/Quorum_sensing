#pragma once
#include <cmath>
#include "particle2D.h"

class EM_BD_iso {
public:
  EM_BD_iso(double h) : h_(h), Dt_(std::sqrt(24. * h)) {}

  template <class TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran) const;

  void set_h(double h) { h_ = h; }

  void set_Dt(double Dt) { Dt_ = Dt; }
protected:
  double h_;
  double Dt_; // sqrt(24 * h) by default
};

template<class TPar, class BdyCondi, class TRan>
void EM_BD_iso::update(TPar& p, const BdyCondi& bc, TRan& myran) const {
  bc.wall_force(p);
  p.pos.x += p.f.x * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += p.f.y * h_ + (myran.doub() - 0.5) * Dt_;
  bc.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
}

class EM_ABD_iso: public EM_BD_iso {
public:
  EM_ABD_iso(double h, double Pe) : EM_BD_iso(h), Dr_(3), sqrt_Dr_(std::sqrt(72. * h)), Pe_(Pe) {}

  template <class TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran) const;

  void set_Pe(double Pe) { Pe_ = Pe; }
  void set_Dr(double Dr) { Dr_ = Dr; sqrt_Dr_ = std::sqrt(24. * h_ * Dr); }
protected:
  double Dr_;
  double sqrt_Dr_; // sqrt(72 * h) by default
  double Pe_;
};

template<class TPar, class BdyCondi, class TRan>
void EM_ABD_iso::update(TPar& p, const BdyCondi& bc, TRan& myran) const {
  const double d_theta = (myran.doub() - 0.5) * sqrt_Dr_;
  Vec_2<double> u = p.update_ori(d_theta);
  bc.wall_force(p);
  p.pos.x += (p.f.x + u.x * Pe_) * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += (p.f.y + u.y * Pe_) * h_ + (myran.doub() - 0.5) * Dt_;
  bc.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
}

class EM_ABD_aniso : public EM_ABD_iso {
public:
  EM_ABD_aniso(double h, double Pe) : EM_ABD_iso(h, Pe) {}

  template <typename TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran) const;
};

template<typename TPar, class BdyCondi, class TRan>
void EM_ABD_aniso::update(TPar& p, const BdyCondi& bc, TRan& myran) const {
  const double d_theta = p.tau * Dr_ * h_ + (myran.doub() - 0.5) * sqrt_Dr_;
  Vec_2<double> u = p.update_ori(d_theta);
  bc.wall_force(p);
  p.pos.x += (p.f.x + u.x * Pe_) * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += (p.f.y + u.y * Pe_) * h_ + (myran.doub() - 0.5) * Dt_;
  bc.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
  p.tau = 0.;
}

class EM_ABD_aniso_Dt0 : public EM_ABD_iso {
public:
  EM_ABD_aniso_Dt0(double h, double Pe) : EM_ABD_iso(h, Pe) {}

  template <typename TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran) const;
};

template<typename TPar, class BdyCondi, class TRan>
void EM_ABD_aniso_Dt0::update(TPar& p, const BdyCondi& bc, TRan& myran) const {
  //const double d_theta = p.tau * Dr_ * h_ + (myran.doub() - 0.5) * sqrt_Dr_;
  const double d_theta = p.tau * h_ + (myran.doub() - 0.5) * sqrt_Dr_;

  Vec_2<double> u = p.update_ori(d_theta);
  bc.wall_force(p);
  p.pos.x += (p.f.x + u.x * Pe_) * h_;
  p.pos.y += (p.f.y + u.y * Pe_) * h_;
  bc.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
  p.tau = 0.;
}

class EM_QS_iso {
public:
  EM_QS_iso(double h, double Dr) : h_(h), sqrt_Dr_(sqrt(24 * Dr * h)) {}

  template <typename TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran, double v) const;

protected:
  double h_;
  double sqrt_Dr_;
};

template<typename TPar, class BdyCondi, class TRan>
void EM_QS_iso::update(TPar& p, const BdyCondi& bc, TRan& myran, double v) const {
  const double d_theta = (myran.doub() - 0.5) * sqrt_Dr_;
  Vec_2<double> u = p.update_ori(d_theta);

  double displace = h_ * v;
  p.pos.x += u.x * displace;
  p.pos.y += u.y * displace;
  bc.tangle(p.pos);
}

  //!TODO Calculate velocity from local density