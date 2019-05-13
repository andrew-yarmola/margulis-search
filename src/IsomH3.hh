#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "AJ.h"
#include "Params.hh"
#include "roundoff.h"
#include "types.hh"

template<typename T>
const T apply_mobius(const SL2<T>& w, const T& z) {
  // TODO: check if we care about division by zero
  return (w.a*z + w.b)/(w.c*z + w.d);
}

template<typename T>
const T four_cosh_re_length(const SL2<T>& w) {
  T tr = w.a + w.d;
  return abs_sqrd(tr) + abs(tr*tr - 4);
}

// TODO: Check that we do need the product of the square roots and
// note the square root of the product.
template<typename T>
const T norm_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  return (tr1*tr1-4)*(tr2*tr2-4);
}

template<typename T>
const T norm(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T sh1 = sqrt((tr1*tr1 - 4));
  T sh2 = sqrt((tr2*tr2 - 4));
  if (absUB(tr1 + sh1) < 2) { sh1 = -sh1; }
  if (absUB(tr2 + sh2) < 2) { sh2 = -sh2; }
  return sh1*sh2;
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2;
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
  return ch*ch - norm_sqrd(w1,w2);
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
  // cosh(2(re(P))x2y2
  const T ch_2_re_perp_n = cosh_2_re_perp_normed(w1,w2); 
  // sinh(2(re(P))x2y2
  const T sh_2_re_perp_n = sinh_2_re_perp_normed(w1,w2); 

  printf("***********************************\n");
  print_type("tr(w1):", tr1);
  print_type("tr(w2):", tr2);
  print_type("x1:", x1);
  print_type("x2:", x2);
  print_type("y1:", y1);
  print_type("y2:", y2);
  print_type("exp(2re(P))x2y2:", e_2_re_perp_n);
  print_type("exp(-2re(P))x2y2:", e_m_2_re_perp_n);
  print_type("cosh(2re(P))x2y2:", ch_2_re_perp_n);
  print_type("sinh(2re(P))x2y2:", sh_2_re_perp_n);
  printf("***********************************\n");

  const T four_cosh_marg = (e_2_re_perp_n*y1 - x1*y2_sqrd)/(e_2_re_perp_n - y2_sqrd) +
             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));

  const T exp_2_t = (x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) + ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))))/(x2_sqrd - e_m_2_re_perp_n);

  std::pair<T,T> result(four_cosh_marg,exp_2_t);

  return result;

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
} 

template<typename T>
const float_pair four_cosh_margulis(const SL2<T>& w1, const SL2<T>& w2, bool upper_margulis, bool upper_t) {
  // TODO Optimize for x and y as w1 (or w2)
  T margulis = T() + pow(2,50);
  T exp_2_t;
  int n = 1;
  int m = 1;
  SL2<T> A(w1);
  printf("%f < %f\n", absLB(four_cosh_re_length(A)), absUB(margulis));
  while (absLB(four_cosh_re_length(A)) < absUB(margulis)) {
    SL2<T> B(w2);
    printf("%f < %f\n", absLB(four_cosh_re_length(B)), absUB(margulis));
    while (absLB(four_cosh_re_length(B)) < absUB(margulis)) {
      std::pair<T,T> m_pair = four_cosh_margulis_simple(A,B);
      T margulis_new = m_pair.first;
      T exp_2_t_new = m_pair.second;
      printf("w1^%d and w2^%d give margulis = %f and exp(2t) = %f\n", n, m, absLB(margulis_new), absLB(exp_2_t_new)); 
      if (absLB(margulis_new) > 0) {
        // We use LB for a partial ordering because we assume that size(margulis) is aobut the same for each comp
        if (absLB(margulis_new) < absLB(margulis)) {
          margulis = margulis_new;
          exp_2_t = exp_2_t_new;
        }
      } else {
        fprintf(stderr, "Failed Margulis computation with error %f\n", absLB(margulis_new)); 
      }
      m += 1;
      B = pow(w2,m); // reducing number of powers needed, might be better to just accumuate
      printf("%f < %f\n", absLB(four_cosh_re_length(B)), absUB(margulis));
    }
    n += 1;
    A = pow(w2,n); // reducing the number of powers needed, might be better to just accumulate
    printf("%f < %f\n", absLB(four_cosh_re_length(A)), absUB(margulis));
  }
  double margulis_f = upper_margulis ? absUB(margulis) : absLB(margulis);
  double exp_2_t_f = upper_t ? absUB(exp_2_t) : absLB(exp_2_t);
  return float_pair(margulis_f, exp_2_t_f); 
}

//template<typename T>
//const bool margulis_is_inner(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
//  // Check is the point realizing the margulis constant always lies in the interior fo the ortholine
//  // TODO : Check signs ERROR ROUNDING
//  T tr1 = w1.a + w1.d;
//  T tr2 = w2.a + w2.d;
//  T x1 = tr1*tr1;
//  T x2 = tr1*tr1 - 4;
//  T y1 = tr2*tr2;
//  T y2 = tr2*tr2 - 4;
//  // cosh(2(re(P)) x2 y2
//  double ch_2_re_perp_LB_normed = cosh_2_re_perp_LB_normed(w1,w2); 
//  double ch_2_re_perp_UB_normed = cosh_2_re_perp_UB_normed(w1,w2); 
//
//  // TODO: For now, box must only contain interior points
//  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
//  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
//  fprintf(stderr, "u + w normed 12 LB : %f, u + w normed 21 LB : %f and coshh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, ch_2_re_perp_LB_normed);
//  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || ch_2_re_perp_LB_normed <= 1) {
//    fprintf(stderr, "Box contains non-interior points\n");
//    return false; // Box contains non-interior points
//  } else {
//    return true;
//  }

template<typename T>
const T cosh_perp_x(const SL2<T>& w) {
  // returns cosh(R) where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T sh = sqrt(tr*tr - 4);
  if (absLB(tr + sh) < 4) { sh = -sh; }
  return td/sh; 
}

template<typename T>
const T cosh_perp_x_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  return td*params.sinhL2*2; 
}

template<typename T>
const T sinh_perp_x_normed(const SL2<T>& w, const Params<T>& params) {
  // TODO: check if correct branch of sqrt
  return sqrt(-(w.b * w.c))*params.sinhL2*4; 
}

template<typename T>
const T cosh_perp_x_sq(const SL2<T>& w) {
  // return cosh(R)^2 where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return (td*td)/(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_x_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  return (td*td)*(params.coshL2*params.coshL2-1)*4; 
}

template<typename T>
const T sinh_perp_x_sq_normed(const SL2<T>& w, const Params<T>& params) {
  return (1 - w.a * w.d)*(params.coshL2*params.coshL2-1)*16; 
}

template<typename T>
const T cosh_perp_y(const SL2<T>& w, const Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)*params.sinhD2*2;
}

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  // TODO check if correct sqrt
  return sqrt(ch*ch-(tr*tr - 4))*params.sinhD2*2;
}

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch = td * params.coshP + dd * params.sinhP;
  return (ch*ch)/(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  return (ch*ch)*(params.coshD2-1)*4;
}

template<typename T>
const T sinh_perp_y_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  return ((ch*ch)-(tr*tr - 4))*(params.coshD2*params.coshD2-1)*4;
}

// We only care about upper bounds on the Jorgenesen inequality. That is,
// If over a box |tr(w1)^2 - 4| + |tr(w1 w2 W1 W2) - 2| < 1, then one of the
// words must be the identity or the box contains no discrete subgroups
template<typename T>
const double jorgensen_UB(const SL2<T>& w1, const SL2<T>& w2) {
  SL2<T> W1 = inverse(w1); 
  SL2<T> W2 = inverse(w2);
  SL2<T> C = w1*w2*W1*W2;
  T tr1 = w1.a + w1.d; 
  T tr2 = C.a + C.d; 
  return (1+EPS)*(absUB(tr1*tr1 - 4) + absUB(tr2 - 2));
}

template<typename T>
const double jorgensen_xw_UB(const SL2<T>& w, const Params<T>& params) {
  T shL2 = params.sinhL2;
  T sh2bc = shL2 * shL2 * w.b * w.c; 
  return 4*((1+EPS)*(absUB(shL2*shL2) + absUB(sh2bc)));
}

template<typename T>
const double jorgensen_wx_UB(const SL2<T>& w, const Params<T>& params) {
  T shL2 = params.sinhL2;
  T sh2bc = shL2 * shL2 * w.b * w.c; 
  T tr = w.a + w.d;
  return (1+EPS)*(absUB(tr*tr - 4) + 4*absUB(sh2bc));
}

template<typename T>
const double jorgensen_yw_UB(const SL2<T>& w, const Params<T>& params) {
  // Note: this seems too complecated
  T shD2 = params.sinhD2;
  T sinhP = params.sinhP;
  T coshP = params.coshP;
  T z = w.a * w.d;
  T f = z*4 - w.a*w.a - w.b*w.b - w.c*w.c - w.d*w.d - 2;
  T g = z*4 + f*sinhP*sinhP + (w.d - w.a)*(w.b - w.c)*sinhP*coshP*2 - 4;
  return (1+EPS)*(4*absUB(shD2*shD2) + absUB(shD2*shD2*g));
}

template<typename T>
const double jorgensen_wy_UB(const SL2<T>& w, const Params<T>& params) {
  // Note: this seems too complecated
  T shD2 = params.sinhD2;
  T sinhP = params.sinhP;
  T coshP = params.coshP;
  T z = w.a * w.d;
  T f = z*4 - w.a*w.a - w.b*w.b - w.c*w.c - w.d*w.d - 2;
  T g = z*4 + f*sinhP*sinhP + (w.d - w.a)*(w.b - w.c)*sinhP*coshP*2 - 4;
  T tr = w.a + w.d;
  return (1+EPS)*(absUB(tr*tr - 4) + absUB(shD2*shD2*g));
}

//template<typename T>
//const double four_cosh_margulis_simple(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
//  // retuns 4 cosh( margulis ) for w1,w2
//  const T tr1 = w1.a + w1.d;
//  const T tr2 = w2.a + w2.d;
//  const T x1 = abs_sqrd(tr1);
//  const T x2 = abs(tr1*tr1 - 4);
//  const T y1 = abs_sqrd(tr2);
//  const T y2 = abs(tr2*tr2 - 4);
//  // exp(2re(P)) 
//  const T ep = e_perp(w1,w2);
//  const T e_2_re_perp = abs_sqrd(ep);
//  // cosh(2(re(P))
//  const T ch_2_re_perp = cosh_2_re_perp(w1,w2); 
//  // sinh(2(re(P))
//  const T sh_2_re_perp = sinh_2_re_perp(w1,w2); 
//
//  print_type("tr(w1):", tr1);
//  print_type("tr(w2):", tr2);
//  print_type("x1:", x1);
//  print_type("x2:", x2);
//  print_type("y1:", y1);
//  print_type("y2:", y2);
//  print_type("exp(P):", ep);
//  print_type("exp(2re(P)):", e_2_re_perp);
//  print_type("cosh(2re(P)):", ch_2_re_perp);
//  print_type("sinh(2re(P)):", sh_2_re_perp);
//
//  T result = (e_2_re_perp*(y1*x2) - (x1*y2))/(e_2_re_perp*x2 - y2) +
//             (sh_2_re_perp*(x2*y2))/((y1-x1)+sqrt((y1-x1)*(y1-x1) + (ch_2_re_perp*2)*(x2*y2) - (x2*x2) - (y2*y2)));
//
//  if (upper) {
//    return absUB(result);
//  } else {
//    return absLB(result);
//  }
//} 


#endif // __IsomH3_h
