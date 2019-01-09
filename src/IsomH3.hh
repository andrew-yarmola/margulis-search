#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Params.hh"
#include "roundoff.h"
#include "types.hh"

template<typename T>
const T cosh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2)/sqrt((tr1*tr1 - 4)*(tr2*tr2 - 4));
};

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2;
};

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z  = td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2; 
  return sqrt(z*z - (tr1*tr1 - 4)*(tr2*tr2 - 4));
};

template<typename T>
const T cosh_perp_sq(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z   = td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2;
  return (z*z)/((tr1*tr1 - 4)*(tr2*tr2 - 4));
};

template<typename T>
const T cosh_perp_sq_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z   = td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2;
  return z*z; 
};

template<typename T>
const T sinh_perp_sq_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z  = td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2; 
  return z*z - (tr1*tr1 - 4)*(tr2*tr2 - 4);
};

template<typename T>
const T cosh_2_perp(const SL2<T>& w1, const SL2<T>& w2) {
  return 2 * cosh_perp_sq(w1,w2) - 1;
};

template<typename T>
const T cosh_perp_x(const SL2<T>& w) {
  // returns cosh(R) where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return td/sqrt(tr*tr - 4); 
};

template<typename T>
const T cosh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return td*params.sinhL2*2; 
};

template<typename T>
const T sinh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  return sqrt(1 - w.a * w.d)*params.sinhL2*4; 
};

template<typename T>
const T cosh_perp_x_sq(const SL2<T>& w) {
  // return cosh(R)^2 where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return (td*td)/(tr*tr - 4); 
};


template<typename T>
const T cosh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return (td*td)*(params.coshL2*params.coshL2-1)*4; 
};

template<typename T>
const T sinh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  return (1 - w.a * w.d)*(params.coshL2*params.coshL2-1)*16; 
};

template<typename T>
const T cosh_perp_y(const SL2<T>& w, Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr*tr - 4); 
};

template<typename T>
const T cosh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)*params.sinhD2*2;
};

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return sqrt(z*z-(tr*tr - 4))*params.sinhD2*2;
};

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return (z*z)/(tr*tr - 4); 
};

template<typename T>
const T cosh_perp_y_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return (z*z)*(params.coshD2-1)*4;
};

template<typename T>
const T sinh_perp_y_sq_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return ((z*z)-(tr*tr - 4))*(params.coshD2*params.coshD2-1)*4;
};

template<typename T>
const double cosh_2_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double sinh_2_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  double ch_LB = cosh_2_re_perp_LB(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(ch_LB*ch_LB)-1),0));
};

template<typename T>
const double cosh_2_re_perp_LB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double sinh_2_re_perp_LB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  double ch_LB = cosh_2_re_perp_LB_normed(w1,w2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4);
  double norm_UB = absUB(x2y2*x2y2);
  printf("norm_UB %f\n", norm_UB);
  printf("cosh_sq_LB %f\n", (1-EPS)*(ch_LB*ch_LB));
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(ch_LB*ch_LB) - norm_UB),0));
};

template<typename T>
const double cosh_2_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
};

template<typename T>
const double sinh_2_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  double ch_UB = cosh_2_re_perp_UB(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(ch_UB*ch_UB)-1),0));
};

template<typename T>
const double cosh_2_re_perp_UB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double sinh_2_re_perp_UB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  double ch_UB = cosh_2_re_perp_UB_normed(w1,w2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4);
  double norm_LB = absLB(x2y2*x2y2);
  printf("norm_LB %f\n", norm_LB);
  printf("cosh_sq_UB %f\n", (1+EPS)*(ch_UB*ch_UB));
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(ch_UB*ch_UB) - norm_LB),0));
};

template<typename T>
const double cosh_2_re_perp_x_LB(const SL2<T>& w) {
  T cp_sq = cosh_perp_x_sq(w);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_x_LB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_x_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_x_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_x_UB(const SL2<T>& w) {
  T cp_sq = cosh_perp_x_sq(w);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_x_UB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_x_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_x_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_y_LB(const SL2<T>& w, Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_y_LB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_y_UB(const SL2<T>& w, Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_y_UB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double four_cosh_margulis(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // retuns 2 cosh( margulis ) for w1,w2
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // Normed cosh and sihn values
  T cpnsq = cosh_perp_sq_normed(w1,w2);
  T spnsq = sinh_perp_sq_normed(w1,w2);
  T e2pn = cpnsq + spnsq + sqrt(spnsq*cpnsq)*2;
  // exp(2re(P)) x2 y2 
  double e_2_re_perp_LB_normed = absLB(e2pn);
  double e_2_re_perp_UB_normed = absUB(e2pn);
  // cosh(2(re(P)) x2 y2
  double ch_2_re_perp_LB_normed = cosh_2_re_perp_LB_normed(w1,w2); 
  double ch_2_re_perp_UB_normed = cosh_2_re_perp_UB_normed(w1,w2); 
  // sinh(2(re(P)) x2 y2
  double sh_2_re_perp_LB_normed = sinh_2_re_perp_LB_normed(w1,w2); 
  double sh_2_re_perp_UB_normed = sinh_2_re_perp_UB_normed(w1,w2);

  // Variables used in formulas paper 
  double al_LB = (1-EPS)*(absLB(y1)-absUB(x1));
  double al_UB = (1+EPS)*(absUB(y1)-absLB(x1));
  
  double beta_LB = (1-EPS)*((1-EPS)*(e_2_re_perp_LB_normed*absLB(y1)) - absUB(y2*y2*x1));
  double beta_UB = (1+EPS)*((1+EPS)*(e_2_re_perp_UB_normed*absUB(y1)) - absLB(y2*y2*x1));

  double kappa_LB = (1-EPS)*(e_2_re_perp_LB_normed - absUB(y2*y2));
  double kappa_UB = (1+EPS)*(e_2_re_perp_UB_normed - absLB(y2*y2));
 
  double eta_LB = (1-EPS)*((1-EPS)*(2*ch_2_re_perp_LB_normed - absUB(x2*x2)) - absUB(y2*y2));
  double eta_UB = (1+EPS)*((1+EPS)*(2*ch_2_re_perp_UB_normed - absLB(x2*x2)) - absLB(y2*y2));
  
  printf("al : %f, %f\n", al_LB, al_UB);
  printf("beta : %f, %f\n", beta_LB, beta_UB);
  printf("kappa : %f, %f\n", kappa_LB, kappa_UB);
  printf("eta : %f, %f\n", eta_LB, eta_UB);
  printf("s : %f, %f\n", sh_2_re_perp_LB_normed, sh_2_re_perp_UB_normed); 
  printf("cosh2rePnormed : %f, %f\n", ch_2_re_perp_LB_normed, ch_2_re_perp_UB_normed);

  // TODO: For now, box must only contain interior points
  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
  printf("u + w normed 12 LB : %f, u + w normed 21 LB : %f and sinh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, sh_2_re_perp_LB_normed);
  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || sh_2_re_perp_LB_normed <= 0) {
    printf("Box contains non-interior points\n");
    return -1; // Box contains non-interior points
  }

  if (al_LB >= 0 || kappa_LB > 0) {
    // TODO: verify that this increasing and decreasing behavior is correct within these bounds
    if (beta_LB < 0 || eta_LB < 0) {
      printf("Error: non-implemented state: beta_LB  < 0 or eta < 0\n");
      return -5;
    }
    if (upper) {
      // return upper bound
      double denom = (1-EPS)*(al_LB + (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(al_LB*al_LB) + eta_LB),0)));
      if (denom == 0) { return -5; }
      return (1+EPS)*((1+EPS)*(beta_UB / kappa_LB) + (1+EPS)*(sh_2_re_perp_UB_normed / denom));
    } else {
      // return lower bound
      double denom = (1+EPS)*(al_UB + (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(al_UB*al_UB) + eta_UB),0)));
      return (1-EPS)*((1-EPS)*(beta_LB / kappa_UB) + (1-EPS)*(sh_2_re_perp_LB_normed / denom));
    }     

  } else if (al_UB <= 0 || kappa_UB < 0) {
    return four_cosh_margulis(w2, w1, upper);
  } else { 
    return -2;
  }
}; 

template<typename T>
const double exp_2_t(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // retuns 2 cosh( margulis ) for w1,w2
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // Normed cosh and sihn values
  T cpnsq = cosh_perp_sq_normed(w1,w2);
  T spnsq = sinh_perp_sq_normed(w1,w2);
  T em2pn = cpnsq + spnsq - sqrt(spnsq*cpnsq)*2;
  // exp(2re(P)) x2 y2 
  double e_minus_2_re_perp_LB_normed = absLB(em2pn);
  double e_minus_2_re_perp_UB_normed = absUB(em2pn);
  // cosh(2(re(P)) x2 y2
  double ch_2_re_perp_LB_normed = cosh_2_re_perp_LB_normed(w1,w2); 
  double ch_2_re_perp_UB_normed = cosh_2_re_perp_UB_normed(w1,w2); 

  // Variables used in formulas paper 
  double delta_LB = (1-EPS)*(absLB(y1*x2) - absUB(x1*x2));
  double delta_UB = (1+EPS)*(absUB(y1*x2) - absLB(x1*x2));
  
  double zeta_LB = (1-EPS)*(absLB(x2*x2) - e_minus_2_re_perp_UB_normed);
  double zeta_UB = (1+EPS)*(absUB(x2*x2) - e_minus_2_re_perp_LB_normed);
 
  double omega_LB = (1-EPS)*((1-EPS)*(2*((1-EPS)*(ch_2_re_perp_LB_normed*absLB(x2*x2))) - absUB(x2*x2*x2*x2)) - absUB(y2*y2*x2*x2));
  double omega_UB = (1+EPS)*((1+EPS)*(2*((1+EPS)*(ch_2_re_perp_UB_normed*absUB(x2*x2))) - absLB(x2*x2*x2*x2)) - absLB(y2*y2*x2*x2));
  
  printf("delta : %f, %f\n", delta_LB, delta_UB);
  printf("zeta : %f, %f\n", zeta_LB, zeta_UB);
  printf("omega : %f, %f\n", omega_LB, omega_UB);
  printf("e_minues_2_perp : %f, %f\n", e_minus_2_re_perp_LB_normed, e_minus_2_re_perp_UB_normed); 

  // TODO: For now, box must only contain interior points
  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
  printf("u + w normed 12 LB : %f, u + w normed 21 LB : %f and coshh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, ch_2_re_perp_LB_normed);
  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || ch_2_re_perp_LB_normed <= 1) {
    printf("Box contains non-interior points\n");
    return -1; // Box contains non-interior points
  }
  if (zeta_LB > 0) {
    if (upper) {
      // return upper bound
      double top = (1+EPS)*(delta_UB + (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(delta_UB*delta_UB) + omega_UB),0)));
      return (1+EPS)*(top/zeta_LB);
    } else {
      // return lower bound
      double top = (1-EPS)*(delta_LB + (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(delta_LB*delta_LB) + omega_LB),0)));
      return (1-EPS)*(top/zeta_UB);
    }     
  } else if (zeta_UB < 0) {
    T e2pn = cpnsq + spnsq + sqrt(spnsq*cpnsq)*2; 
    if (upper) {
      return (1+EPS)*(absUB(e2pn)/exp_2_t(w2,w1,false));
    } else {
      return (1-EPS)*(absLB(e2pn)/exp_2_t(w2,w1,true));
    }
  } else { // u = v somehwere in the box... which is bad
    return -3;
  }
}; 

#endif // __IsomH3_h
