#include "IsomH3.hh"
#include "roundoff.h"

double max(double a, double b) {
  if (a < b) return b;
  else return a;
}

template<typename T>
const T cosh_perp(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))/sqrt((tr1*tr1-4)*(tr2*tr2-4));
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b);
}

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z  = td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b); 
  return sqrt(z*z - (tr1*tr1 - 4)*(tr2*tr2 - 4));
}

template<typename T>
const T cosh_perp_sq(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z   = td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b);
  return (z*z)/((tr1*tr1 - 4)*(tr2*tr2 - 4));
}

template<typename T>
const T cosh_perp_sq_normed(const SL2<T>& w1, SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z   = td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b);
  return z*z; 
}

template<typename T>
const T sinh_perp_sq_normed(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  T z  = td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b); 
  return z*z - (tr1*tr1 - 4)*(tr2*tr2 - 4);
}

template<typename T>
const T cosh_2_perp(const SL2<T>& w1, SL2<T>& w2) {
  return 2 * cosh_perp_sq(w1,w2) - 1;
}

template<typename T>
const T cosh_perp_x(const SL2<T>& w) {
  // returns cosh(R) where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return td/sqrt(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return 2*td*params.sinhL2; 
}

template<typename T>
const T sinh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  return 4*sqrt(1 - w.a * w.d)*params.sinhL2; 
}

template<typename T>
const T cosh_perp_x_sq(const SL2<T>& w) {
  // return cosh(R)^2 where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return (td*td)/(tr*tr - 4); 
}


template<typename T>
const T cosh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return 4*(td*td)*(params.coshL2*params.coshL2-1); 
}

template<typename T>
const T sinh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  return 16*(1 - w.a * w.d)*(params.coshL2*params.coshL2-1); 
}

template<typename T>
const T cosh_perp_y(const SL2<T>& w, Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return 2*(td * params.coshP + dd * params.sinhP)*params.sinhD2;
}

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 2*sqrt(z*z-(tr*tr - 4))*params.sinhD2;
}

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return (z*z)/(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 4*(z*z)*(params.coshD2-1);
}

template<typename T>
const T sinh_perp_y_sq_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 4*((z*z)-(tr*tr - 4))*(params.coshD2*params.coshD2-1);
}

template<typename T>
const double cosh_2_re_perp_LB(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
}

template<typename T>
const double sinh_2_re_perp_LB(const SL2<T>& w1, SL2<T>& w2) {
  double ch_LB = cosh_2_re_perp_LB(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(ch_LB*ch_LB)-1),0));
}

template<typename T>
const double cosh_2_re_perp_LB_normed(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
}

template<typename T>
const double sinh_2_re_perp_LB_normed(const SL2<T>& w1, SL2<T>& w2) {
  double ch_LB = cosh_2_re_perp_LB(w1,w2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  double norm_UB = absUB((tr1*tr1 - 4)(tr2*tr2 - 4));
  // Lemma 7.0 in GMT
  return (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(ch_LB*ch_LB) - norm_UB),0));
}

template<typename T>
const double cosh_2_re_perp_UB(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
}

template<typename T>
const double sinh_2_re_perp_UB(const SL2<T>& w1, SL2<T>& w2) {
  double ch_UB = cosh_2_re_perp_UB(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(ch_UB*ch_UB)-1),0));
}

template<typename T>
const double cosh_2_re_perp_UB_normed(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
}

template<typename T>
const double sinh_2_re_perp_UB_normed(const SL2<T>& w1, SL2<T>& w2) {
  double ch_UB = cosh_2_re_perp_UB(w1,w2);
  T tr1 = w1.a + w1.d;
  T tr2 = w1.a + w1.d;
  double norm_LB = absLB((tr1*tr1 - 4)(tr2*tr2 - 4));
  // Lemma 7.0 in GMT
  return (1+EPS)*sqrt(max((1+EPS)*((1+EPS)*(ch_UB*ch_UB) - norm_LB),0));
}

template<typename T>
const double cosh_2_re_perp_x_LB(const SL2<T>& w) {
  T cp_sq = cosh_perp_x_sq(w);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_x_LB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_x_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_x_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
}

template<typename T>
const double cosh_2_re_perp_x_UB(const SL2<T>& w) {
  T cp_sq = cosh_perp_x_sq(w);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_x_UB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_x_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_x_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
}

template<typename T>
const double cosh_2_re_perp_y_LB(const SL2<T>& w, Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_y_LB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
}

template<typename T>
const double cosh_2_re_perp_y_UB(const SL2<T>& w, Params<T>& params) {
  T cp_sq = cosh_perp_y_sq(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_y_UB_normed(const SL2<T>& w, Params<T>& params) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
}

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
  T cpn = cosh_perp_normed(w1,w2);
  T spn = sinh_perp_normed(w1,w2);
  T e2pn = (cpn + spn) * (cpn + spn);
  // exp(2re(P)) x2 y2 
  double e_2_re_perp_LB_normed = absLB(e2pn);
  double e_2_re_perp_UB_normed = absUB(e2pn);
  // cosh(2(re(P)) x2 y2
  double ch_2_re_perp_LB_normed = cosh_2_re_perp_UB_normed(w1,w2); 
  double ch_2_re_perp_UB_normed = cosh_2_re_perp_UB_normed(w1,w2); 
  // sinh(2(re(P)) x2 y2
  double sh_2_re_perp_LB_normed = sinh_2_re_perp_LB_normed(w1,w2); 
  double sh_2_re_perp_UB_normed = sinh_2_re_perp_UB_normed(w1,w2);

  // TODO: For now, box must only contain interior points
  double upw_LB = (1-EPS)*((1-EPS)*((1-EPS)*((1-EPS)*(ch_2_re_perp_LB_normed*absLB(y2)) + absLB(y1)) - absUB(x2)) - absUB(x1));
  if (upw_LB <= 0 || sh_2_re_perp_LB_normed <= 0) {
    return -1; // Box contains non-interior points
  }

  // Variables used in formulas paper 
  double al_LB = (1-EPS)*(absLB(y1)-absUB(x1));
  double al_UB = (1+EPS)*(absUB(y1)-absLB(x1));
  
  double beta_LB = (1-EPS)*((1-EPS)*(e_2_re_perp_LB_normed*absLB(y1)) - absUB(y2*y2*x1));
  double beta_UB = (1+EPS)*((1+EPS)*(e_2_re_perp_UB_normed*absUB(y1)) - absLB(y2*y2*x1));

  double kappa_LB = (1-EPS)*(e_2_re_perp_LB_normed - absUB(y2*y2));
  double kappa_UB = (1+EPS)*(e_2_re_perp_UB_normed - absLB(y2*y2));
 
  double eta_LB = (1-EPS)*((1-EPS)*(2*ch_2_re_perp_LB_normed - absUB(x2*x2)) - absUB(y2*y2));
  double eta_UB = (1+EPS)*((1+EPS)*(2*ch_2_re_perp_UB_normed - absLB(x2*x2)) - absLB(y2*y2));

  if (al_LB >= 0 || kappa_LB > 0) {
    // TODO: verify that this increasing and decreasing behavior is correct within these bounds
    if (upper) {
      // return upper bound
      double denom = (1-EPS)*(al_LB + (1-EPS)*sqrt(max((1-EPS)*((1-EPS)*(al_LB*al_LB) + eta_LB),0)));
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
} 

//template<typename T>
//const double cosh_2_re_tube_UB(const SL2<T>& w1, const SL2<T>& w2) {
//} 

/*
const ACJpair fixed_points(const SL2ACJ& w) {
  // Returns the pair of fixed point of the Mobius transform w
  ACJ tr = w.a + w.d;
  ACJ td = w.a - w.d;
  return ACJpair((td + sqrt(tr*tr - 4))/(2*w.c),
                 (td - sqrt(tr*tr - 4))/(2*w.c));
}

const ACJ cosh_perp(const ACJpair& p1, ACJpair& p2) {
  return  ((p1.x + p1.y)(p2.x + p2.y) - 2 * p1.x * p2.y - 2 * p2.x * p2.y)/((p1.x - p1.y)(p2.x - p2.y));
}

const std::Cpair fixed_points(const SL2C& w); 
  // Returns the pair of fixed point of the Mobius transform w
  XComplex tr = (w.a + w.d).z;
  XComplex s = (w.a - w.d).z;
  return Cpair(((s + sqrt(((tr*tr)).z-4).z)/(2*w.c)).z,
               ((s - sqrt(((tr*tr)).z-4).z)/(2*w.c)).z);
}*/
