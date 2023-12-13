#include <iostream>
#include "vect.h"
#include "rand.h"
#include "exporter2D.h"


class Bubble {
public:
  Bubble(Vec_2<double> pos0, double r0) : pos(pos0), radius(r0) {}

  template<typename TRan>
  Bubble(TRan& myran, double radius0, const Vec_2<double>& L)
    : pos(myran.doub()* L.x, myran.doub()* L.y), radius(radius0) {}

  int coalescence(Bubble& b2, const Vec_2<double>& L);

  template <typename TRan>
  void diffuse(TRan& myran, const Vec_2<double>& L, double h);

  Vec_2<double> pos;
  double radius;
  bool status = true;
};


int Bubble::coalescence(Bubble& b2, const Vec_2<double>& L) {
  int flag_colliding = 0;
  if (status && b2.status) {
    Vec_2<double> dr = b2.pos - pos;
    Vec_2<double> half_L = L * 0.5;
    if (dr.x < -half_L.x) {
      dr.x += L.x;
    } else if (dr.x > half_L.x) {
      dr.x -= L.x;
    }
    if (dr.y < -half_L.y) {
      dr.y += L.y;
    } else if (dr.y > half_L.y) {
      dr.y -= L.y;
    }

    double r2 = dr.square();
    double d = radius + b2.radius;
    if (r2 < d * d) {
      flag_colliding = 1;
      b2.status = false;
      pos += 0.5 * (b2.pos - pos);
      if (pos.x < 0) {
        pos.x += L.x;
      } else if (pos.x >= L.x) {
        pos.x -= L.x;
      }
      if (pos.y < 0) {
        pos.y += L.y;
      } else if (pos.y >= L.y) {
        pos.y -= L.y;
      }
      radius = sqrt(radius * radius + b2.radius * b2.radius);
    }
  }
  return flag_colliding;
}

template <typename TRan>
void Bubble::diffuse(TRan& myran, const Vec_2<double>& L, double h) {
  // D ~ 1/R
  // double D = sqrt(0.1/ radius * 24 * h);

  // D ~ 1/R^3
  //double D = sqrt(0.001 / (radius * radius * radius) * 24 * h);

  // D ~ 1/R^5.7
  double D = sqrt(pow(0.1/radius, 5.7) * 24 * h);

  pos.x += (myran.doub() - 0.5) * D;
  pos.y += (myran.doub() - 0.5) * D;
  if (pos.x < 0) {
    pos.x += L.x;
  } else if (pos.x >= L.x) {
    pos.x -= L.x;
  }
  if (pos.y < 0) {
    pos.y += L.y;
  } else if (pos.y >= L.y) {
    pos.y -= L.y;
  }

}


void remove_dead_bubbles(std::vector<Bubble>& p_arr, int n_dead) {
  int n = p_arr.size();
  int n_removed = 0;
  for (int i = n - 1; i > -1; i--) {
    if (!p_arr[i].status) {
      if (i == n - 1) {
        p_arr.pop_back();
      } else {
        p_arr[i] = p_arr.back();
        p_arr.pop_back();
      }
      n_removed++;
      if (n_removed == n_dead) {
        break;
      }
    }
  }
}

int main(int argc, char* argv[]) {
  Ranq2 myran(1);

  Vec_2<double> L(400, 400);
  double r0 = 0.1;
  double h = 0.1;
  int N0 = 16000;
  int n_step = 100000;
  std::vector<Bubble> p_arr;



  exporter::OrderParaExporter_2 out("output4.dat", 0, n_step, 100);

  p_arr.reserve(N0);
  for (int i = 0; i < N0; i++) {
    p_arr.emplace_back(myran, r0, L);
  }

  double packing_frac = 3.14159 * r0 * r0 * N0 / (L.x * L.y);
  double sum_r_square = r0 * r0 * N0;
  std::cout << "packing fraction = " << packing_frac << std::endl;

  for (int t = 1; t <= n_step; t++) {
    for (auto& p : p_arr) {
      p.diffuse(myran, L, h);
    }

    int coalescence_num = 0;
    int n = p_arr.size();
    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        coalescence_num += p_arr[i].coalescence(p_arr[j], L);
      }
    }
    remove_dead_bubbles(p_arr, coalescence_num);

    double r_new = sqrt(sum_r_square / p_arr.size());
    for (auto& p : p_arr) {
      p.radius = r_new;
    }

    out.dump(t, p_arr.size(), r_new);



  }
  




}