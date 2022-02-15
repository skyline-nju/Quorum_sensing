#pragma once
#include <cmath>
#include "particle2D.h"
#include "node.h"

class BD_EM {
public:
  BD_EM(double h)
    : h_(h), Dt_(std::sqrt(24. * h)), Dr_(std::sqrt(72. * h)) {}

  template <class TPar, class TDomain, class TRan>
  void update(TPar& p, const TDomain& dm, TRan& myran) const;

  template <class TPar, class TDomain, class TRan>
  void update(TPar& p, double Pe, const TDomain& dm, TRan& myran) const;

  template <class TDomain, class TRan>
  void update(BiNode<BP_u_tau_2>& p, double Pe, const TDomain& dm, TRan& myran) const;

protected:
  double h_;
  double Dt_; // sqrt(24 * h) by default
  double Dr_; // sqrt(72 * h) by default
};

template <class TPar, class TDomain, class TRan>
void BD_EM::update(TPar& p, const TDomain& dm, TRan& myran) const{
  p.pos.x += p.f.x * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += p.f.y * h_ + (myran.doub() - 0.5) * Dt_;
  dm.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
}

template<class TPar, class TDomain, class TRan>
void BD_EM::update(TPar& p, double Pe, const TDomain& dm, TRan& myran) const{
  const double d_theta = (myran.doub() - 0.5) * Dr_;
  Vec_2<double> u = p.update_ori(d_theta);
  p.pos.x += (p.f.x + u.x * Pe) * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += (p.f.y + u.y * Pe) * h_ + (myran.doub() - 0.5) * Dt_;
  dm.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
}

template<class TDomain, class TRan>
void BD_EM::update(BiNode<BP_u_tau_2>& p, double Pe, const TDomain& dm, TRan& myran) const {
  const double d_theta = p.tau * 3. * h_ + (myran.doub() - 0.5) * Dr_;
  Vec_2<double> u = p.update_ori(d_theta);
  p.pos.x += (p.f.x + u.x * Pe) * h_ + (myran.doub() - 0.5) * Dt_;
  p.pos.y += (p.f.y + u.y * Pe) * h_ + (myran.doub() - 0.5) * Dt_;
  dm.tangle(p.pos);
  p.f.x = 0.;
  p.f.y = 0.;
  p.tau = 0.;
}

class EM_QS_iso {
public:
  EM_QS_iso(double h, double Dt, double Dr) 
    : h_(h), Dt_(sqrt(24 * Dt * h)), Dr_(sqrt(24 * Dr * h)), trans_noise_on_(Dt > 0.) {}

  template <typename TPar, class BdyCondi, class TRan>
  void update(TPar& p, const BdyCondi& bc, TRan& myran, double v) const;

protected:
  double h_;
  double Dt_; // sqrt(24 * h) by default
  double Dr_; // sqrt(72 * h) by default
  bool trans_noise_on_;
};

template<typename TPar, class BdyCondi, class TRan>
void EM_QS_iso::update(TPar& p, const BdyCondi& bc, TRan& myran, double v) const {
  const double d_theta = (myran.doub() - 0.5) * Dr_;
  Vec_2<double> u = p.update_ori(d_theta);

  double displace = h_ * v;
  p.pos.x += u.x * displace;
  p.pos.y += u.y * displace;
  if (trans_noise_on_) {
    p.pos.x += (myran.doub() - 0.5) * Dt_;
    p.pos.y += (myran.doub() - 0.5) * Dt_;
  }
  bc.tangle(p.pos);
}

template<typename TPar>
double cal_v(double v0, double kappa, double rho0, double eta, double JAB, double JBA, TPar& p) {
  double v;
  double d_rhoA = (p.rho_local[0] - rho0) / (rho0 * kappa);
  double d_rhoB = (p.rho_local[1] - rho0) / (rho0 * kappa);

  if (p.type_id == 0) {
    v = v0 * (1 + kappa * tanh(eta * d_rhoA)) * (1 + kappa * tanh(JAB * d_rhoB));
  } else {
    v = v0 * (1 + kappa * tanh(eta * d_rhoB)) * (1 + kappa * tanh(JBA * d_rhoA));
  }
  return v;
}
