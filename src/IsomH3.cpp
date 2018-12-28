#include "IsomH3.h"

template<typename T>
const T cosh_perp(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))/sqrt((tr1^2-4)*(tr2^2-4));
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b));
}

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return sqrt((td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))^2 - (tr1^2-4)*(tr2^2-4));
}

template<typename T>
const T cosh_perp_sq(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))^2/((tr1^2-4)*(tr2^2-4));
}

template<typename T>
const T cosh_perp_sq_normed(const SL2<T>& w1, SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))^2;
}

template<typename T>
const T sinh_perp_sq_normed(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return (td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b))^2 - (tr1^2-4)*(tr2^2-4);
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
  return td/sqrt(tr^2-4); 
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
  return td^2/(tr^2-4); 
}


template<typename T>
const T cosh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return 4*td^2*(params.coshL2^2-1); 
}

template<typename T>
const T sinh_perp_x_sq_normed(const SL2<T>& w, Params<T>& params) {
  return 16*(1 - w.a * w.d)*(params.coshL2^2-1); 
}

template<typename T>
const T cosh_perp_y(const SL2<T>& w, Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr^2-4); 
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
  return 2*sqrt((td * params.coshP + dd * params.sinhP)^2-(tr^2-4))*params.sinhD2;
}

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)^2/(tr^2-4); 
}

template<typename T>
const T cosh_perp_y_sq_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return 4*(td * params.coshP + dd * params.sinhP)^2*(params.coshD2-1);
}

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return 4*((td * params.coshP + dd * params.sinhP)^2-(tr^2-4))*(params.coshD2^2-1);
}

template<typename T>
const double cosh_2_re_perp_LB(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq) + absLB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_LB_normed(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1-EPS)*(absLB(cp_sq_normed) + absLB(sp_sq_normed));
}

template<typename T>
const double cosh_2_re_perp_UB(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq = cosh_perp_sq(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq) + absUB(cp_sq-1));
}

template<typename T>
const double cosh_2_re_perp_UB_normed(const SL2<T>& w1, SL2<T>& w2) {
  T cp_sq_normed = cosh_perp_sq_normed(w1,w2);
  T sp_sq_normed = sinh_perp_sq_normed(w1,w2);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
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
const double cosh_2_re_perp_x_UB_normed(const SL2<T>& w) {
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
const double cosh_2_re_perp_y_UB_normed(const SL2<T>& w) {
  T cp_sq_normed = cosh_perp_y_sq_normed(w, params);
  T sp_sq_normed = sinh_perp_y_sq_normed(w, params);
  // Lemma 7.0 in GMT
  return (1+EPS)*(absUB(cp_sq_normed) + absUB(sp_sq_normed));
}

template<typename T>
const double cosh_2_re_tube_LB(const SL2<T>& w1, const SL2<T>& w2) {
  // retuns cosh(2t) where t is the distance to magulis point of w1,w2 from axis(w1)
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d; 
  T tr2 = w2.a + w2.d;
  double c2rpLB = cosh_2_re_perp_LB(w1,w2);
  double s2rpLB = sqrt(c2rpLB^2 - 1);
  double c2rpUB = cosh_2_re_perp_UB(w1,w2);
  double s2rpUB = sqrt(c2rpUB^2 - 1);
  double uLB = c2rpLB * absLB(tr2^2 - 4) - absUB(tr1^2 -4);
  double vLB = s2rpLB * absLB(tr2^2 - 4);
  double wLB = absLB(tr2^2) - absUB(tr1^2);
  double uUB = c2rpUB * absUB(tr2^2 - 4) - absLB(tr1^2 -4);
  double vUB = s2rpUB * absUB(tr2^2 - 4);
  double wUB = absUB(tr2^2) - absLB(tr1^2);
  if (wLB < 0) return -1; // either the box is too big or we need to swtich geods
  if (wLB^2 + vLB^2 - uUB^2 < 0) return -1; // sqrt in the final formula may not be positive - mox may contain critical pt
  if (uLB^2 < vUB^2 || vLB^2 < uUB^2) return -1; // denominator may be zero in the box -> may contain criical point
  double u,v,w; // pick for obtaining (estimating for now) LB
  if ((uUb <= 0 && 0 < uLB + wLB) || ( 0 <= uLB && uUB <= vLB)) { // inc function of w
    w = wLB;
  } else 
  if ((0 < uLB && uUB <= wLV && vUB < uLB) || (wUB^2 < uLB^2 && uUB^2 < vLB^2 + wLB^2) { // dec function of w
    w = wUB;
  } else {
    return -1; // box might contain critical point of this function
  }
  if (0 < uLB + wLB && uUB^2 < vLB^2 + wLB^2) { // inc function of w
    u = uLB;
  } else 
  if (uUB + wUB < 0 && uUB^2 < vLB^2 + wLB^2) { // dec function of w
    u = uUB;
  } else {
    return -1; // box might contain critical point of this function
  }
  // for v we only dec function or constant wrt to v
  v = vUB;
  return (u*w + v * sqrt(w^2 + v^2 - u^2))/(v^2 - u^2);
} 

template<typename T>
const double cosh_2_re_tube_UB(const SL2<T>& w1, const SL2<T>& w2) {
  // retuns cosh(2t) where t is the distance to magulis point of w1,w2 from axis(w1)
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d; 
  T tr2 = w2.a + w2.d;
  double c2rpLB = cosh_2_re_perp_LB(w1,w2);
  double s2rpLB = sqrt(c2rpLB^2 - 1);
  double c2rpUB = cosh_2_re_perp_UB(w1,w2);
  double s2rpUB = sqrt(c2rpUB^2 - 1);
  double uLB = c2rpLB * absLB(tr2^2 - 4) - absUB(tr1^2 -4);
  double vLB = s2rpLB * absLB(tr2^2 - 4);
  double wLB = absLB(tr2^2) - absUB(tr1^2);
  double uUB = c2rpUB * absUB(tr2^2 - 4) - absLB(tr1^2 -4);
  double vUB = s2rpUB * absUB(tr2^2 - 4);
  double wUB = absUB(tr2^2) - absLB(tr1^2);
  if (wLB < 0) return -1; // either the box is too big or we need to swtich geods
  if (wLB^2 + vLB^2 - uUB^2 < 0) return -1; // sqrt in the final formula may not be positive - mox may contain critical pt
  if (uLB^2 < vUB^2 || vLB^2 < uUB^2) return -1; // denominator may be zero in the box -> may contain criical point
  double u,v,w; // pick for obtaining (estimating for now) LB
  if ((uUb <= 0 && 0 < uLB + wLB) || ( 0 <= uLB && uUB <= vLB)) { // inc function of w
    w = wUB;
  } else 
  if ((0 < uLB && uUB <= wLV && vUB < uLB) || (wUB^2 < uLB^2 && uUB^2 < vLB^2 + wLB^2) { // dec function of w
    w = wLB;
  } else {
    return -1; // box might contain critical point of this function
  }
  if (0 < uLB + wLB && uUB^2 < vLB^2 + wLB^2) { // inc function of w
    u = uUB;
  } else 
  if (uUB + wUB < 0 && uUB^2 < vLB^2 + wLB^2) { // dec function of w
    u = uLB;
  } else {
    return -1; // box might contain critical point of this function
  }
  // for v we only dec function or constant wrt to v
  v = vLB;
  return (u*w + v * sqrt(w^2 + v^2 - u^2))/(v^2 - u^2);
} 

/*
const ACJpair fixed_points(const SL2ACJ& w) {
  // Returns the pair of fixed point of the Mobius transform w
  ACJ tr = w.a + w.d;
  ACJ td = w.a - w.d;
  return ACJpair((td + sqrt(tr^2-4))/(2*w.c),
                 (td - sqrt(tr^2-4))/(2*w.c));
}

const ACJ cosh_perp(const ACJpair& p1, ACJpair& p2) {
  return  ((p1.x + p1.y)(p2.x + p2.y) - 2 * p1.x * p2.y - 2 * p2.x * p2.y)/((p1.x - p1.y)(p2.x - p2.y));
}

const std::Cpair fixed_points(const SL2C& w); 
  // Returns the pair of fixed point of the Mobius transform w
  XComplex tr = (w.a + w.d).z;
  XComplex s = (w.a - w.d).z;
  return Cpair(((s + sqrt((tr^2).z-4).z)/(2*w.c)).z,
               ((s - sqrt((tr^2).z-4).z)/(2*w.c)).z);
}*/

