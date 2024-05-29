#pragma once
#include "vect.h"
#include <cmath>
#include "communicator2D.h"
#include "config.h"

class LinearDensityKernal {
public:
  LinearDensityKernal(double r_cut) : r_cut_(r_cut), r_cut_square_(r_cut* r_cut) {}

  template <typename TPar>
  void accum_density(double r12, TPar& p1, TPar& p2) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec);

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
  p1.update_local_density(p2, density_weighted);
}

template<typename TPar>
void LinearDensityKernal::operator()(TPar& p1, TPar& p2, const Vec_2<double>& r12_vec) {
  const double r12_square = r12_vec.square();
  if (r12_square < r_cut_square_) {
    double r12_module = sqrt(r12_square);
    accum_density(r12_module, p1, p2);
  }
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
    p.norm_local_density(norm);
  }
}

template <typename TPar>
void LinearDensityKernal::reset_local_density(std::vector<TPar>& p_arr) const {
  for (auto& p : p_arr) {
    p.set_local_density(0.);
  }
}


template <typename TNode, typename TFunc>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl,
               Communicator_2& comm, TFunc for_all_pair_force) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}
