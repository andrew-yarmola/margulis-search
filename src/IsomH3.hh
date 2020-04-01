#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "AJ.h"
#include "Generators.hh"
#include "roundoff.h"
#include "types.hh"
#include "assert.h"

template<typename T>
const T four_cosh_re_length(const SL2<T>& w) {
  T tr = w.a + w.d;
  return abs_sqrd(tr) + abs(tr*tr - 4);
}

template<typename T>
const T norm_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  return (tr1 * tr1 - 4) * (tr2 * tr2 - 4);
}

// TODO: Check that we do need the product of the square roots and
// not the square root of the product.
template<typename T>
const T norm(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T sh1 = sqrt((tr1 * tr1 - 4));
  T sh2 = sqrt((tr2 * tr2 - 4));
  if (absUB(tr1 + sh1) < 2) { sh1 = -sh1; }
  if (absUB(tr2 + sh2) < 2) { sh2 = -sh2; }
  return sh1 * sh2;
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + (w1.b * w2.c + w1.c * w2.b) * 2;
}

template<typename T>
const T cosh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  return cosh_perp_normed(w1,w2)/norm(w1,w2);
}

template<typename T>
const T cosh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return z*z; 
}

template<typename T>
const T cosh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  return cosh_perp_sqrd_normed(w1,w2)/norm_sqrd(w1,w2);
}

template<typename T>
const T abs_cosh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return abs_sqrd(z); 
}

template<typename T>
const T abs_cosh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  return abs_cosh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2)); // TODO: check if this gives best error as both sqrt and division make error blow up
}

template<typename T>
const T sinh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp_normed(w1,w2); 
  return ch * ch - norm_sqrd(w1,w2);
}

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp_normed(w1,w2); 
  T n_sqrd = norm_sqrd(w1,w2);
  T sh = sqrt(ch * ch - n_sqrd);
  if (absUB(ch + sh) < absUB(norm(w1,w2))) { sh = -sh; }
  return sh;
}

template<typename T>
const T sinh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp(w1,w2);
  T sh = sqrt(ch * ch - 1);
  if (absUB(ch + sh) < 1) { sh = -sh; }
  return sh;
}

template<typename T>
const T abs_sinh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  return abs(sinh_perp_sqrd_normed(w1,w2));
}

template<typename T>
const T abs_sinh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  return abs_sinh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2)); // TODO: check if this gives best error
}

template<typename T>
const T cosh_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T abs_cp_sqrd = abs_cosh_perp_sqrd(w1,w2);
  T abs_sp_sqrd = abs_sinh_perp_sqrd(w1,w2);
  return abs_cp_sqrd + abs_sp_sqrd;
}

template<typename T>
const T cosh_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T abs_cp_sqrd_normed = abs_cosh_perp_sqrd_normed(w1,w2);
  T abs_sp_sqrd_normed = abs_sinh_perp_sqrd_normed(w1,w2);
  return abs_cp_sqrd_normed + abs_sp_sqrd_normed;
}

template<typename T>
const T sinh_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_2_re_perp(w1,w2);
  T sh = sqrt(ch * ch - 1);
  if (absUB(ch + sh) < 1) { sh = -sh; }
  return sh;
}

template<typename T>
const T sinh_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_2_re_perp_normed(w1,w2);
  T n = norm_sqrd(w1,w2);
  T abs_n_sqrd = abs_sqrd(n);
  T sh = sqrt(ch * ch - abs_n_sqrd);
  if (absUB(ch + sh) < absUB(n)) { sh = -sh; }
  return sh;
}

template<typename T>
const T e_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp(w1,w2);
  T sp = sinh_perp(w1,w2);
  return cp + sp;
}

template<typename T>
const T e_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  return abs_sqrd(e_perp(w1,w2));
}

template<typename T>
const T e_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp_normed(w1,w2);
  T sp = sinh_perp_normed(w1,w2);
  return abs_sqrd(cp + sp);
}

template<typename T>
const T e_m_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp(w1,w2);
  T sp = sinh_perp(w1,w2);
  return cp - sp;
}

template<typename T>
const T e_m_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp_normed(w1,w2);
  T sp = sinh_perp_normed(w1,w2);
  return abs_sqrd(cp - sp);
}

template<typename T>
const double e_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  return absUB(e_perp(w1,w2));
}

template<typename T>
const double e_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  return absLB(e_perp(w1,w2));
}

template<typename T>
const std::pair<T,T> four_cosh_margulis_simple(const SL2<T>& w1, const SL2<T>& w2) {
  /*
  print_center("w1.a:", w1.a);
  print_center("w1.b:", w1.b);
  print_center("w1.c:", w1.c);
  print_center("w1.d:", w1.d);
  print_center("w2.a:", w2.a);
  print_center("w2.b:", w2.b);
  print_center("w2.c:", w2.c);
  print_center("w2.d:", w2.d);
  */
  // retuns 4 cosh( margulis ) and exp(2t) for w1,w2
  const T tr1 = w1.a + w1.d;
  const T tr2 = w2.a + w2.d;
  const T x1 = abs_sqrd(tr1);
  const T x2 = abs(tr1*tr1 - 4);
  const T x2_sqrd = abs_sqrd(tr1*tr1 - 4);
  const T y1 = abs_sqrd(tr2);
  const T y2 = abs(tr2*tr2 - 4);
  const T y2_sqrd = abs_sqrd(tr2*tr2 - 4);
  // exp(2re(P))x2y2
  const T e_2_re_perp_n = e_2_re_perp_normed(w1,w2);
  const T e_m_2_re_perp_n = e_m_2_re_perp_normed(w1,w2);
  const T e_2_re_p = e_2_re_perp(w1,w2);
  // cosh(2(re(P))x2y2
  const T ch_2_re_perp_n = cosh_2_re_perp_normed(w1,w2); 
  // sinh(2(re(P))x2y2
  const T sh_2_re_perp_n = sinh_2_re_perp_normed(w1,w2); 

  /*
  printf("***********************************\n");
  print_center("tr(w1):", tr1);
  print_center("tr(w2):", tr2);
  print_center("x1:", x1);
  print_center("x2:", x2);
  print_center("y1:", y1);
  print_center("y2:", y2);
  print_center("cosh(perp)sqrt((tr1^2-4)(tr2^2-4)):", cosh_perp_normed(w1,w2));
  print_center("sinh(perp)sqrt((tr1^2-4)(tr2^2-4)):", sinh_perp_normed(w1,w2));
  print_center("exp(2re(perp))x2y2:", e_2_re_perp_n);
  print_center("exp(-2re(perp))x2y2:", e_m_2_re_perp_n);
  print_center("cosh(2re(perp))x2y2:", ch_2_re_perp_n);
  print_center("sinh(2re(perp))x2y2:", sh_2_re_perp_n);
  printf("***********************************\n");
  print_center("e_2_re_perp_n*y1 - x1*y2_sqrd:",e_2_re_perp_n*y1 - x1*y2_sqrd);
  print_center("e_2_re_perp_n - y2_sqrd:", e_2_re_perp_n - y2_sqrd);
  print_center("a/b:", (e_2_re_perp_n*y1 - x1*y2_sqrd)/(e_2_re_perp_n - y2_sqrd));
  print_center("sinh(2re(perp))x2y2:", sh_2_re_perp_n);
  print_center("(y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)):", (y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));
  print_center("a/b:", sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));
  print_center("x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))):", x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));
  print_center("x2_sqrd - e_m_2_re_perp_n:", x2_sqrd - e_m_2_re_perp_n);
  print_center("a/b:", (x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))))/(x2_sqrd - e_m_2_re_perp_n));
  printf("***********************************\n");
  */

  std::vector<T> versions;
  versions.push_back((e_2_re_perp_n*y1 - x1*y2_sqrd)/(e_2_re_perp_n - y2_sqrd) +
             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));

  // swap x and y to hope for better error (formula is symmetric, but not in an obious way)
  versions.push_back((e_2_re_perp_n*x1 - y1*x2_sqrd)/(e_2_re_perp_n - x2_sqrd) +
             sh_2_re_perp_n/((x1-y1)+sqrt((x1-y1)*(x1-y1) + ((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd))));
  
// Alternate version where first denominator is scaled by e^(-2re(P))y2/x2. The current version should keep the denominator larger
//  const T four_cosh_marg_v1_alt1 = (y1*x2_sqrd - x1*e_m_2_re_perp_n)/(x2_sqrd - e_m_2_re_perp_n) +
//             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));

// Alternate version where we split off x1. The current version should keep error slightly smaller (since fewer additions)
//  const T four_cosh_marg_v1_alt2 = (x1 + ((y1-x1)*x2_sqrd)/(x2_sqrd - e_m_2_re_perp_n)) +
//             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));

  // Collect everything as one fraction. Helps when re(P) is very close to zero 
  versions.push_back(x1 + ((y1-x1)*(e_2_re_perp_n-x2_sqrd) - sh_2_re_perp_n*((y1-x1) -
              sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))))/((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd));

  // swap x and y to hope for better error (formula is symmetric, but not in an obious way)
  versions.push_back(y1 + ((x1-y1)*(e_2_re_perp_n-y2_sqrd) - sh_2_re_perp_n*((x1-y1) -
              sqrt((x1-y1)*(x1-y1) + ((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd))))/((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd));

  /* 
  print_type("4cosh(margulis) v1:", versions[0]);
  print_type("4cosh(margulis) v2:", versions[1]);
  print_type("4cosh(margulis) v3:", versions[2]);
  print_type("4cosh(margulis) v4:", versions[3]);
  */

  std::sort(versions.begin(), versions.end(), sort_comp<T>);

  T four_cosh_marg = versions[0];

  T exp_2_t = (x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))))/(x2_sqrd - e_m_2_re_perp_n);
  /*
  print_center("exp(2re(perp)):", e_2_re_p);
  print_center("exp(2t):", (const T) exp_2_t);
  print_center("exp(2t) - exp(2re(perp)):", exp_2_t - e_2_re_p);
  print_center("exp(2re(perp)) exp(2t):", e_2_re_p - exp_2_t);
  */

  // if we realize the margulis number outside the ortho segement, it means it is on one of the end point
  const T one = T(1);
  const T zero = T(0);
  T fc_marg_lb = max_local(x2 + x1, y2 + y1);
  // print_type("lower bound for four cosh margulis:", fc_marg_lb);
  if (absUB(exp_2_t) < 1.0) {
    four_cosh_marg = x2 + x1;
    exp_2_t = one;
  } else if (absLB(exp_2_t - 1.0) == 0) {
    four_cosh_marg = max_local(four_cosh_marg, fc_marg_lb);
    exp_2_t = max_local(exp_2_t, one);
  } else if (absLB(exp_2_t) > absUB(e_2_re_p)) {
      four_cosh_marg = y2 + y1;
      exp_2_t = e_2_re_p;  
  } else if (absLB(exp_2_t - e_2_re_p) == 0) {
    four_cosh_marg = max_local(four_cosh_marg, fc_marg_lb);
    exp_2_t = min_local(exp_2_t, e_2_re_p);
  }
  if (absLB(y2 + y1) > absUB((e_2_re_p*x2)+x1)) { 
    four_cosh_marg = y2+y1;
    exp_2_t = e_2_re_p; 
  }
  if (absLB(x2 + x1) > absUB((e_2_re_p*y2)+y1)) { 
    four_cosh_marg = x2+x1;
    exp_2_t = zero; 
  }

  //  print_type("4cosh(margulis) final:", four_cosh_marg);
  //  print_type("exp_2_t final:", exp_2_t);

  std::pair<T,T> result(four_cosh_marg, exp_2_t);

  return result;
}


template<typename T>
const T cosh_move_j(const SL2<T>& w) {
    T q = abs_sqrd(w.c) + abs_sqrd(w.d);
    T z = w.a * conj(w.c) + w.b * conj(w.d);
    return (abs_sqrd(z) + (q - 1) * (q - 1))/(q * 2) + 1; 
}

// We compute |tr(w1)^2 - 4| + |tr(w1 w2 W1 W2) - 2|
// with optimzation for x and y specifically
template<typename T>
const T jorgensen(const SL2<T>& w1, const SL2<T>& w2) {
  SL2<T> W1 = inverse(w1); 
  SL2<T> W2 = inverse(w2);
  SL2<T> C = w1*w2*W1*W2;
  T tr1 = w1.a + w1.d; 
  T tr2 = C.a + C.d; 
  return abs(tr1*tr1 - 4) + abs(tr2 - 2);
}

// Eliminate bad boxes that can't generate non-elementay groups
template<typename T>
const T jorgensen_xy(const Params<T>& p) {
  T z = p.sinhLy2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhLx2) * 4;
}

// Eliminate bad boxes that can't generate non-elementay groups
template<typename T>
const T jorgensen_yx(const Params<T>& p) {
  T z = p.sinhLx2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhLy2) * 4;
}

template<typename T>
const T jorgensen_xw(const SL2<T>& w, const Params<T>& p) {
  T shLx2 = p.sinhLx2;
  T td = w.a - w.d;
  T z = w.c * p.expmdx - w.b * p.expdx;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLx2);
}

template<typename T>
const T jorgensen_wx(const SL2<T>& w, const Params<T>& p) {
  T shLx2 = p.sinhLx2;
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T z = w.c * p.expmdx - w.b * p.expdx;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLx2);
}

template<typename T>
const T jorgensen_yw(const SL2<T>& w, const Params<T>& p) {
  T shLy2 = p.sinhLy2;
  T td = w.a - w.d;
  T z = w.c * p.expdyf - w.b * p.expmdyf;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLy2);
}

template<typename T>
const T jorgensen_wy(const SL2<T>& w, const Params<T>& p) {
  T shLy2 = p.sinhLy2;
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T z = w.c * p.expdyf - w.b * p.expmdyf;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLy2);
}

// Distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ax_wax(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expmdx - w.b * p.expdx;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Distance between axis(y) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ay_way(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expdyf - w.b * p.expmdyf;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ax_way(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
        ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
  return  abs(z - 2) + abs(z + 2);
}

// Distance between axis(y) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ay_wax(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expmdx - (w.b * w.b) * p.expdx) * p.expmdyf +
        ((w.d * w.d) * p.expdx  - (w.c * w.c) * p.expmdx) * p.expdyf;
  return  abs(z - 2) + abs(z + 2);
}

#endif // __IsomH3_h

//  printf("###################################\n");
//  print_type("4 cosh(margulis) :", four_cosh_marg);
//  print_type("exp(2t) :", exp_2_t);
//  printf("absLB(exp(2t)) = %f and absUB(exp(2t)) = %f :\n", absLB(exp_2_t), absUB(exp_2_t));
//  printf("###################################\n");

//  float_pair result;
//
//  if (upper_margulis) {
//    result.first = absUB(four_cosh_marg);
//  } else {
//    result.first = absLB(four_cosh_marg);
//  }
//  if (upper_t) {
//    result.second = absUB(exp_2_t);
//  } else {
//    result.second = absLB(exp_2_t);
//  }
//  return result;

/*
#define MAX_LOOPS 10000

template<typename T>
const float_pair four_cosh_margulis(const SL2<T>& w1, const SL2<T>& w2, bool upper_margulis, bool upper_t) {
  // TODO Optimize for x and y as w1 (or w2)
  T margulis = T() + pow(2,50);
  T exp_2_t;
  int n = 1;
  SL2<T> A(w1);
  // print_SL2(A);
  printf("LB 4cosh(re(A)) vs UB margulis: %f < %f\n", absLB(four_cosh_re_length(A)), absUB(margulis));
  int loops = 0;
  while (absLB(four_cosh_re_length(A)) < absUB(margulis)) {
    int m = 1;
    SL2<T> B(w2);
    // print_SL2(B);
    printf("LB 4cosh(re(B)) vs UB margulis: %f < %f\n", absLB(four_cosh_re_length(B)), absUB(margulis));
    while (absLB(four_cosh_re_length(B)) < absUB(margulis)) {
      if (loops > MAX_LOOPS) { break; }
      std::pair<T,T> m_pair = four_cosh_margulis_simple(A,B);
      T margulis_new = m_pair.first;
      T exp_2_t_new = m_pair.second;
      printf("w1^%d and w2^%d give margulis = %f and exp(2t) = %f\n", n, m, absLB(margulis_new), absLB(exp_2_t_new)); 
      if (absLB(margulis_new) == 0) {
        fprintf(stderr, "Margulis LB is zero!\n");
        print_type("Got:", (const T) margulis_new);
      }
      // We use LB for a partial ordering because we assume that size(margulis) is about the same for each comp
      if (absLB(margulis_new) < absLB(margulis) || (absLB(margulis_new) == absLB(margulis) && absUB(margulis_new) < absUB(margulis))) {
        margulis = margulis_new;
        exp_2_t = exp_2_t_new;
      }
//      } else {
//        double margulis_f = upper_margulis ? infinity() : 0.0;
//        double exp_2_t_f = upper_t ? infinity() : 0.0;
//        return float_pair(margulis_f, exp_2_t_f); 
//      }
      m += 1;
      B = pow(w2, m); // reducing number of powers needed, might be better to just accumuate
      loops += 1;
      // print_SL2(B);
      printf("LB 4cosh(re(B)) vs UB margulis: %f < %f\n", absLB(four_cosh_re_length(B)), absUB(margulis));
    }
    if (loops > MAX_LOOPS) { break; }
    n += 1;
    A = pow(w1, n); // reducing the number of powers needed, might be better to just accumulate
    // print_SL2(A);
    printf("LB 4cosh(re(A)) vs UB margulis: %f < %f\n", absLB(four_cosh_re_length(A)), absUB(margulis));
  }
  double margulis_f = upper_margulis ? absUB(margulis) : absLB(margulis);
  double exp_2_t_f = upper_t ? absUB(exp_2_t) : absLB(exp_2_t);
  return float_pair(margulis_f, exp_2_t_f); 
}
*/

