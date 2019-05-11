#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Params.hh"
#include "roundoff.h"
#include "types.hh"

template<typename T>
const T apply_mobius(const SL2<T>& w, const T& z) {
  // TODO: check if we care about division by zero
  return (w.a*z + w.b)/(w.c*z + w.d);
}




template<typename T>
const double four_cosh_re_length_LB(const SL2<T>& w) {
  T tr = w.a + w.d;
  return (1-EPS)*(absLB(tr*tr) + absLB(tr*tr - 4));
};

template<typename T>
const double four_cosh_re_length_UB(const SL2<T>& w) {
  T tr = w.a + w.d;
  return (1+EPS)*(absUB(tr*tr) + absUB(tr*tr - 4));
};

template<typename T>
const T cosh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T sh1 = sqrt((tr1*tr1 - 4));
  T sh2 = sqrt((tr2*tr2 - 4));
  if (absLB(tr1 + sh1) < 4) { sh1 = -sh1; }
  if (absLB(tr2 + sh2) < 4) { sh2 = -sh2; }
  return (td1 * td2 + (w1.b * w2.c + w1.c * w2.b)*2)/(sh1 * sh2);
};

template<typename T>
const T sinh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp(w1,w2);
  T sh = sqrt(ch * ch - 1);
  if (absLB(ch + sh) < 1) { sh = -sh; }
  return sh;
};

template<typename T>
const double e_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  return absUB(cosh_perp(w1,w2) + sinh_perp(w1,w2));
};

template<typename T>
const double e_re_perp_x_UB(const SL2<T>& w) {
  // Note : 2 a d - 1 = +/- cosh(dist((0,infty),(w . 0, w . infrty)))
  T z = w.a * w.d;
  T sh = sqrt(z*(z-1));
  if (absLB((z + sh)*2 - 1) < 1) { sh = -sh; }
  return absUB((z + sh) * 2 - 1);
};

template<typename T>
const double e_re_perp_y_UB(const SL2<T>& w, const T& coshP, const T& sinhP) {
  // Note: this seems too complecated
  T z = w.a * w.d;
  T f = z*4 - w.a*w.a - w.b*w.b - w.c*w.c - w.d*w.d - 2;
  T ch = z*2 + f*sinhP*sinhP/2 + (w.d - w.a)*(w.b - w.c)*sinhP*coshP - 1;
  fprintf(stderr, "Cosh LB %f and UB %f\n", absLB(ch), absUB(ch));
  T sh = sqrt(ch*ch-1);
  fprintf(stderr, "Sinh LB %f and UB %f\n", absLB(sh), absUB(sh));
  if (absLB(ch + sh) < 1) { sh = -sh; }
  return absUB(ch + sh);
};

template<typename T>
const double e_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  return absLB(cosh_perp(w1,w2) + sinh_perp(w1,w2));
};

template<typename T>
const double e_re_perp_x_LB(const SL2<T>& w) {
  // Note : 2 a d - 1 = +/- cosh(dist((0,infty),(w . 0, w . infrty)))
  T z = w.a * w.d;
  T sh = sqrt(z*(z-1));
  if (absLB((z + sh)*2 - 1) < 1) { sh = -sh; }
  return absLB((z + sh) * 2 - 1);
};

template<typename T>
const double e_re_perp_y_LB(const SL2<T>& w, T& coshP, T& sinhP) {
  // Note: this seems too complecated
  T z = w.a * w.d;
  T f = z*4 - w.a*w.a - w.b*w.b - w.c*w.c - w.d*w.d - 2;
  T ch = z*2 + f*sinhP*sinhP/2 + (w.d - w.a)*(w.b - w.c)*sinhP*coshP - 1;
  fprintf(stderr, "Cosh LB %f and UB %f\n", absLB(ch), absUB(ch));
  T sh = sqrt(ch*ch-1);
  fprintf(stderr, "Sinh LB %f and UB %f\n", absLB(sh), absUB(sh));
  if (absLB(ch + sh) < 1) { sh = -sh; }
  return absLB(ch + sh);
};

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + (w1.b * w2.c + w2.b * w1.c)*2;
};

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T sh1 = sqrt((tr1*tr1 - 4));
  T sh2 = sqrt((tr2*tr2 - 4));
  if (absLB(tr1 + sh1) < 4) { sh1 = -sh1; }
  if (absLB(tr2 + sh2) < 4) { sh2 = -sh2; }
  T norm = (tr1*tr1 - 4)*(tr2*tr2 - 4);
  T ch = td1 * td2 + (w1.b * w2.c + w2.b * w1.c)*2; 
  T sh = sqrt(ch * ch - norm);
  if (absLB((ch + sh)/(sh1*sh2)) < 1) { sh = -sh; }
  return sh;
};

template<typename T>
const T cosh_perp_sq(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T ch = td1 * td2 + (w1.b * w2.c + w2.b * w1.c)*2; 
  return (ch*ch)/((tr1*tr1 - 4)*(tr2*tr2 - 4));
};

template<typename T>
const T cosh_perp_sq_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T ch = td1 * td2 + (w1.b * w2.c + w2.b * w1.c)*2; 
  return ch*ch; 
};

template<typename T>
const T sinh_perp_sq_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z = td1 * td2 + (w1.b * w2.c + w2.b * w1.c)*2; 
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
  T sh = sqrt(tr*tr - 4);
  if (absLB(tr + sh) < 4) { sh = -sh; }
  return td/sh; 
};

template<typename T>
const T cosh_perp_x_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  return td*params.sinhL2*2; 
};

template<typename T>
const T sinh_perp_x_normed(const SL2<T>& w, const Params<T>& params) {
  // TODO: check if correct branch of sqrt
  return sqrt(-(w.b * w.c))*params.sinhL2*4; 
};

template<typename T>
const T cosh_perp_x_sq(const SL2<T>& w) {
  // return cosh(R)^2 where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return (td*td)/(tr*tr - 4); 
};


template<typename T>
const T cosh_perp_x_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  return (td*td)*(params.coshL2*params.coshL2-1)*4; 
};

template<typename T>
const T sinh_perp_x_sq_normed(const SL2<T>& w, const Params<T>& params) {
  return (1 - w.a * w.d)*(params.coshL2*params.coshL2-1)*16; 
};

template<typename T>
const T cosh_perp_y(const SL2<T>& w, const Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr*tr - 4); 
};

template<typename T>
const T cosh_perp_y_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)*params.sinhD2*2;
};

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  // TODO check if correct sqrt
  return sqrt(ch*ch-(tr*tr - 4))*params.sinhD2*2;
};

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch = td * params.coshP + dd * params.sinhP;
  return (ch*ch)/(tr*tr - 4); 
};

template<typename T>
const T cosh_perp_y_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  return (ch*ch)*(params.coshD2-1)*4;
};

template<typename T>
const T sinh_perp_y_sq_normed(const SL2<T>& w, const Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T ch  = td * params.coshP + dd * params.sinhP;
  return ((ch*ch)-(tr*tr - 4))*(params.coshD2*params.coshD2-1)*4;
};

template<typename T>
const double cosh_2_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double sinh_2_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  double c2rp_sq_LB = (1-EPS)*((1-EPS)*(absLB(cp_sq*cp_sq) + absLB((cp_sq-1)*(cp_sq-1))) + absLB(cp_sq*(cp_sq - 1))*2);  
  // Lemma 7.0 in GMT and sqrt > 0 so we should be OK
  return (1-EPS)*sqrt(max((1-EPS)*(c2rp_sq_LB - 1),0));
};

template<typename T>
const double cosh_2_re_perp_LB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_LB_normed_f(const SL2<T>& w1, const SL2<T>& w2, const T& factor) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2)*factor;
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2)*factor;
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double sinh_2_re_perp_LB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cpsqn = cosh_perp_sq_normed(w1,w2);
  T spsqn = sinh_perp_sq_normed(w1,w2);
  double ch_sq_LB = (1-EPS)*((1-EPS)*(absLB(cpsqn*cpsqn) + absLB(spsqn*spsqn)) + absLB(spsqn*cpsqn)*2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4);
  double norm_UB = absUB(x2y2*x2y2);
  fprintf(stderr, "norm_UB %f\n", norm_UB);
  fprintf(stderr, "cosh_sq_LB %f\n", ch_sq_LB);
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*(ch_sq_LB - norm_UB),0));
};

template<typename T>
const double sinh_2_re_perp_LB_normed_f(const SL2<T>& w1, const SL2<T>& w2, const T& factor) {
  T cpsqn = cosh_perp_sq_normed(w1,w2)*factor;
  T spsqn = sinh_perp_sq_normed(w1,w2)*factor;
  double ch_sq_LB = (1-EPS)*((1-EPS)*(absLB(cpsqn*cpsqn) + absLB(spsqn*spsqn)) + absLB(spsqn*cpsqn)*2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4)*factor*factor;
  double norm_UB = absUB(x2y2*x2y2);
  fprintf(stderr, "norm_UB %f\n", norm_UB);
  fprintf(stderr, "cosh_sq_LB %f\n", ch_sq_LB);
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*(ch_sq_LB - norm_UB),0));
};

template<typename T>
const double cosh_2_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
};

template<typename T>
const double sinh_2_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  double c2rp_sq_UB = (1+EPS)*((1+EPS)*(absUB(cp_sq*cp_sq) + absUB((cp_sq-1)*(cp_sq-1))) + absUB(cp_sq*(cp_sq - 1))*2);  
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*(c2rp_sq_UB-1),0));
};

template<typename T>
const double cosh_2_re_perp_UB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_UB_normed_f(const SL2<T>& w1, const SL2<T>& w2, const T& factor) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2)*factor;
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2)*factor;
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double sinh_2_re_perp_UB_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cpsqn = cosh_perp_sq_normed(w1,w2);
  T spsqn = sinh_perp_sq_normed(w1,w2);
  double ch_sq_UB = (1+EPS)*((1+EPS)*(absUB(cpsqn*cpsqn) + absUB(spsqn*spsqn)) + absUB(spsqn*cpsqn)*2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4);
  double norm_UB = absUB(x2y2*x2y2);
  fprintf(stderr, "norm_UB %f\n", norm_UB);
  fprintf(stderr, "cosh_sq_UB %f\n", ch_sq_UB);
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*(ch_sq_UB - norm_UB),0));
};

template<typename T>
const double sinh_2_re_perp_UB_normed_f(const SL2<T>& w1, const SL2<T>& w2, const T& factor) {
  T cpsqn = cosh_perp_sq_normed(w1,w2)*factor;
  T spsqn = sinh_perp_sq_normed(w1,w2)*factor;
  double ch_sq_UB = (1+EPS)*((1+EPS)*(absUB(cpsqn*cpsqn) + absUB(spsqn*spsqn)) + absUB(spsqn*cpsqn)*2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  T x2y2 = (tr1*tr1 - 4)*(tr2*tr2 - 4)*factor*factor;
  double norm_UB = absUB(x2y2*x2y2);
  fprintf(stderr, "norm_UB %f\n", norm_UB);
  fprintf(stderr, "cosh_sq_UB %f\n", ch_sq_UB);
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*(ch_sq_UB - norm_UB),0));
};

template<typename T>
const double cosh_2_re_perp_x_LB(const SL2<T>& w) {
  T cp_sq = cosh_perp_x_sq(w);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_x_LB_normed(const SL2<T>& w, const Params<T>& params) {
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
const double cosh_2_re_perp_x_UB_normed(const SL2<T>& w, const Params<T>& params) {
  T cp_sq_normed = cosh_perp_x_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_x_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_y_LB(const SL2<T>& w, const Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_y_LB_normed(const SL2<T>& w, const Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
};

template<typename T>
const double cosh_2_re_perp_y_UB(const SL2<T>& w, const Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
};

template<typename T>
const double cosh_2_re_perp_y_UB_normed(const SL2<T>& w, const Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
};

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
};

template<typename T>
const double jorgensen_xw_UB(const SL2<T>& w, const Params<T>& params) {
  T shL2 = params.sinhL2;
  T sh2bc = shL2 * shL2 * w.b * w.c; 
  return 4*((1+EPS)*(absUB(shL2*shL2) + absUB(sh2bc)));
};

template<typename T>
const double jorgensen_wx_UB(const SL2<T>& w, const Params<T>& params) {
  T shL2 = params.sinhL2;
  T sh2bc = shL2 * shL2 * w.b * w.c; 
  T tr = w.a + w.d;
  return (1+EPS)*(absUB(tr*tr - 4) + 4*absUB(sh2bc));
};

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
};

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
};

template<typename T>
const double four_cosh_margulis_simple(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // retuns 2 cosh( margulis ) for w1,w2
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // Normed cosh and sihn values
  T cpn = cosh_perp_normed(w1,w2);
  T spn = sinh_perp_normed(w1,w2);
  T cpnsq = cosh_perp_sq_normed(w1,w2);
  T spnsq = sinh_perp_sq_normed(w1,w2);
  T e2pn = cpnsq + spnsq + (spn * cpn) * 2;
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
  
  double beta_LB = (1-EPS)*(absLB(e2pn*y1) - absUB(y2*y2*x1));
  double beta_UB = (1+EPS)*(absUB(e2pn*y1) - absLB(y2*y2*x1));

  double kappa_LB = (1-EPS)*(e_2_re_perp_LB_normed - absUB(y2*y2));
  double kappa_UB = (1+EPS)*(e_2_re_perp_UB_normed - absLB(y2*y2));
 
  double eta_LB = (1-EPS)*((1-EPS)*(2*ch_2_re_perp_LB_normed - absUB(x2*x2)) - absUB(y2*y2));
  double eta_UB = (1+EPS)*((1+EPS)*(2*ch_2_re_perp_UB_normed - absLB(x2*x2)) - absLB(y2*y2));
  
  fprintf(stderr, "al : %f, %f\n", al_LB, al_UB);
  fprintf(stderr, "beta : %f, %f\n", beta_LB, beta_UB);
  fprintf(stderr, "kappa : %f, %f\n", kappa_LB, kappa_UB);
  fprintf(stderr, "eta : %f, %f\n", eta_LB, eta_UB);
  fprintf(stderr, "s : %f, %f\n", sh_2_re_perp_LB_normed, sh_2_re_perp_UB_normed); 
  fprintf(stderr, "cosh2rePnormed : %f, %f\n", ch_2_re_perp_LB_normed, ch_2_re_perp_UB_normed);

  // TODO: For now, box must only contain interior points
  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
  fprintf(stderr, "u + w normed 12 LB : %f, u + w normed 21 LB : %f and sinh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, sh_2_re_perp_LB_normed);
  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || sh_2_re_perp_LB_normed <= 0) {
    fprintf(stderr, "Box contains non-interior points\n");
    return -1; // Box contains non-interior points
  }

  if (al_LB >= 0 || kappa_LB > 0) {
    // TODO: verify that this increasing and decreasing behavior is correct within these bounds
    if (beta_LB < 0 || eta_LB < 0) {
      fprintf(stderr, "Error: non-implemented state: beta_LB  < 0 or eta < 0\n");
      return -3;
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
    return four_cosh_margulis_simple(w2, w1, upper);
  } else { 
    return -2;
  }
}; 

template<typename T>
const double four_cosh_re_length(const SL2<T>& w, bool upper) {
  if (upper) return four_cosh_re_length_UB(w);
  else return four_cosh_re_length_LB(w);
}

template<typename T>
const double four_cosh_margulis(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // TODO Optimize for x and y as w1 (or w2)
  double margulis = pow(2,100);
  int n = 1;
  int m = 1;
  SL2<T> A(w1);
  while (four_cosh_re_length(A,upper) < margulis) {
    SL2<T> B(w2);
    while (four_cosh_re_length(B,upper) < margulis) {
      double margulis_new = four_cosh_margulis_simple(A,B,upper);
      if (margulis_new >= 0) {
        margulis = fmin(margulis, margulis_new);
      } else {
        fprintf(stderr, "Failed Margulis computation with error %f\n", margulis_new); 
      }
      m += 1;
      B = pow(w2,m); // reducing number of powers needed, might be better to just accumuate
    }
    n += 1;
    A = pow(w2,n); // reducing the number of powers needed, might be better to just accumulate
  }
  return margulis; 
}

template<typename T>
const double exp_2_t(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // retuns exp(2t) bounds where t is the distance along orth(w1,w2) from axis(w1) to margulis point
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // Normed cosh and sihn values
  T cpn = cosh_perp_normed(w1,w2);
  T spn = sinh_perp_normed(w1,w2);
//  fprintf(stderr, "Cosh perp %f + %f I and Sinh perp %f + %f I\n", cpn.f.re, cpn.f.im, spn.f.re, spn.f.im);
  T cpnsq = cosh_perp_sq_normed(w1,w2);
  T spnsq = sinh_perp_sq_normed(w1,w2);
//  fprintf(stderr, "Cosh^2 perp %f + %f I and Sinh^2 perp %f + %f I\n", cpnsq.f.re, cpnsq.f.im, spnsq.f.re, spnsq.f.im);
  T em2pn = cpnsq + spnsq - (spn * cpn) * 2;
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
 
  double omega_LB = (1-EPS)*((1-EPS)*(2*cosh_2_re_perp_LB_normed_f(w1,w2,x2*x2) - absUB(x2*x2*x2*x2)) - absUB(y2*y2*x2*x2));
  double omega_UB = (1+EPS)*((1+EPS)*(2*cosh_2_re_perp_UB_normed_f(w1,w2,x2*x2) - absLB(x2*x2*x2*x2)) - absLB(y2*y2*x2*x2));
  
  fprintf(stderr, "delta : %f, %f\n", delta_LB, delta_UB);
  fprintf(stderr, "zeta : %f, %f\n", zeta_LB, zeta_UB);
  fprintf(stderr, "omega : %f, %f\n", omega_LB, omega_UB);
  fprintf(stderr, "e_minues_2_perp : %f, %f\n", e_minus_2_re_perp_LB_normed, e_minus_2_re_perp_UB_normed); 

  // TODO: For now, box must only contain interior points
  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
  fprintf(stderr, "u + w normed 12 LB : %f, u + w normed 21 LB : %f and coshh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, ch_2_re_perp_LB_normed);
  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || ch_2_re_perp_LB_normed <= 1) {
    fprintf(stderr, "Box contains non-interior points\n");
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
    // TODO : replace with a better formula for exp(2p-2t)
    T e2pn = cpnsq + spnsq + sqrt(spnsq*cpnsq)*2; 
    if (upper) {
      return (1+EPS)*(absUB(e2pn)/exp_2_t(w2,w1,false));
    } else {
      return (1-EPS)*(absLB(e2pn)/exp_2_t(w2,w1,true));
    }
  } else { // u = v somehwere in the box... which is bad
    return -2;
  }
}; 

template<typename T>
const bool margulis_is_inner(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // Check is the point realizing the margulis constant always lies in the interior fo the ortholine
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // cosh(2(re(P)) x2 y2
  double ch_2_re_perp_LB_normed = cosh_2_re_perp_LB_normed(w1,w2); 
  double ch_2_re_perp_UB_normed = cosh_2_re_perp_UB_normed(w1,w2); 

  // TODO: For now, box must only contain interior points
  double u_plus_w_12_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(y1*x2)) - absUB(x2*x2)) - absUB(x1*x2));
  double u_plus_w_21_LB = (1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed + absLB(x1*y2)) - absUB(y2*y2)) - absUB(y1*y2));
  fprintf(stderr, "u + w normed 12 LB : %f, u + w normed 21 LB : %f and coshh2RePx2y2 LB %f\n", u_plus_w_12_LB, u_plus_w_21_LB, ch_2_re_perp_LB_normed);
  if (u_plus_w_12_LB <= 0 || u_plus_w_21_LB <= 0 || ch_2_re_perp_LB_normed <= 1) {
    fprintf(stderr, "Box contains non-interior points\n");
    return false; // Box contains non-interior points
  } else {
    return true;
  }
}

#endif // __IsomH3_h
