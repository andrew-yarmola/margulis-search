#ifndef __Params_h
#define __Params_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "AJ.h"
#include "Generators.hh"
#include "roundoff.h"
#include "types.hh"
#include "assert.h"

// Exp of distance from ortho endpts to the special point on axis
template<typename T>
double exp_dist_to_ortho_x(const T& zm, const T& zp, const Params<T>& p) {
  T one(1.0); 
  T z = ((p.expdx * zm + one) * (p.expdx * zp + one))/
    ((p.expdx * zm - one) * (p.expdx * zp - one));
  return absUB(sqrt(z));
}

// Exp of distance from ortho endpts to the special point on axis
template<typename T>
double exp_dist_to_ortho_y(const T& zm, const T& zp, const Params<T>& p) {
  T z = ((p.expdyf + zm) * (p.expdyf + zp))/
    ((p.expdyf - zm) * (p.expdyf - zp));
  return absUB(sqrt(z));
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
//  print_type(shLx2);
//  print_type(td);
//  print_type(z);
//  T s = abs_sqrd(shLx2);
//  print_type(s);
//  T m = abs(td * td - z * z) + 4;
//  print_type(m);
//  T ans = (abs(td * td - z * z) + 4) * abs_sqrd(shLx2);
//  print_type(ans);
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

// Complex distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_sinh_perp2_sq_ax_wax(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expmdx - w.b * p.expdx;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// Distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ax_wax(const SL2<T>& w, const Params<T>& p) {
  T four_sinh_sq_perp2 = four_sinh_perp2_sq_ax_wax(w, p);
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Complex distance between axis(y) and w(axis(y)) 
template<typename T>
const T four_sinh_perp2_sq_ay_way(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expdyf - w.b * p.expmdyf;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// Distance between axis(y) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ay_way(const SL2<T>& w, const Params<T>& p) {
  T four_sinh_sq_perp2 = four_sinh_perp2_sq_ay_way(w, p); 
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Complex distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_sinh_perp2_sq_ax_way(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
    ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
  return z - 2;
}

// Distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ax_way(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
    ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
  return  abs(z - 2) + abs(z + 2);
}

// Complex distance between axis(y) and w(axis(x)) 
template<typename T>
const T four_sinh_perp2_sq_ay_wax(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expmdx - (w.b * w.b) * p.expdx ) * p.expmdyf +
    ((w.d * w.d) * p.expdx  - (w.c * w.c) * p.expmdx) * p.expdyf;
  return  z - 2;
}

// Distance between axis(y) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ay_wax(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expmdx - (w.b * w.b) * p.expdx ) * p.expmdyf +
    ((w.d * w.d) * p.expdx  - (w.c * w.c) * p.expmdx) * p.expdyf;
  return  abs(z - 2) + abs(z + 2);
}

#endif // __Params_h
