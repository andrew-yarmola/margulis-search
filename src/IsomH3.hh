#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Params.h"

template<typename T>
const T cosh_perp(const SL2<T>& w1, SL2<T>& w2);

template<typename T>
const T cosh_perp_sq(const SL2<T>& w1, SL2<T>& w2);

template<typename T>
const T cosh_2_perp(const SL2<T>& w1, SL2<T>& w2);

template<typename T>
const T cosh_perp_x(const SL2<T>& w);

template<typename T>
const T cosh_perp_x_sq(const SL2<T>& w);

template<typename T>
const T cosh_perp_y(const SL2<T>& w, Params<T>& params);

template<typename T>
const T cosh_perp_y_sq(const SL2<T>& w, Params<T>& params);

template<typename T>
const double cosh_2_re_perp_LB(const SL2<T>& w1, SL2<T>& w2);

template<typename T>
const double cosh_2_re_perp_UB(const SL2<T>& w1, SL2<T>& w2);

template<typename T>
const double cosh_2_re_perp_x_LB(const SL2<T>& w);

template<typename T>
const double cosh_2_re_perp_x_UB(const SL2<T>& w);

template<typename T>
const double cosh_2_re_perp_y_LB(const SL2<T>& w, Params<T>& params);

template<typename T>
const double cosh_2_re_perp_y_UB(const SL2<T>& w, Params<T>& params);

template<typename T>
const double cosh_2_re_tube_LB(const SL2<T>& w1, const SL2<T>& w2);

template<typename T>
const double cosh_2_re_tube_UB(const SL2<T>& w1, const SL2<T>& w2);

#endif // __IsomH3_h
