#ifndef __Params_h
#define __Params_h
#include <math.h>
#include <string>
#include "SL2.hh"

/*
* We will work in sinh space. The parameters of interest are
* Let x and y realize the margulis constant such that re(length(x)) =< re(length(y))
* cosh(L/2) where L = length(x)
* cosh(D/2) where D = length(y)
* cosh(P) where P is the length of the ortholine from axis(x) to axis(y)
* ASSUMPTION : L,D have imaginary part between -pi and pi, and P has imaginary part 
* between -pi/2 and pi/2. This is to guarnatee that sqrt(cosh^2(t)-1) = sinh(t) for
* t = L/2,D/2, and P, where the branch of sqrt is take to be the negative real axis.
*/

template<typename T> struct Params {
	T sinhP;
	T coshP;
	T sinhD2;
	T coshD2; // may not always be filled
	T sinhL2; // may not always be filled
	T coshL2; // may not always be filled
};

/*
* Let x and y realize the margulis constant such that re(length(x)) =< re(length(y))
* We will conjugate x to have axis(x) = {0,\infty}
*/

template<typename T>
SL2<T> construct_x(const Params<T>& params) {
  // We want the matrix {{ exp(L/2), 0 }, { 0, exp(-L/2) }}
  T zero = T(0);
  T sl2 = params.sinhL2;
  T cl2 = params.coshL2;
	return SL2<T>(cl2 + sl2, zero,
                zero, cl2 - sl2);
};

template<typename T>
SL2<T> construct_y(const Params<T>& params) {
  // We want the matrix {{ cosh(D/2) + cosh(P)sinh(D/2), -sinh(D/2)sinh(P) }, 
  //                     { sinh(D/2)sinh(P), cosh(D/2) - cosh(P)sinh(D/2) }}
  T sd2 = params.sinhD2; 
  T cd2 = params.coshD2;
  T sp = params.sinhP;
  T cp = params.coshP;
	return SL2<T>(cd2 + cp * sd2, -(sd2 * sp),
                sd2 * sp, cd2 - cp * sd2);
};

//int x_power(std::string w) {
//  int count = 0;
//  for (std::string::size_type p = 0; p < w.size(); ++p) {
//      if (w[p] == 'x' || w[p] == 'X') ++count;
//  }
//  return count;
//}; 
//
//int y_power(std::string w) {
//  int count = 0;
//  for (std::string::size_type p = 0; p < w.size(); ++p) {
//      if (w[p] == 'y' || w[p] == 'Y') ++count;
//  }
//  return count;
//}; 
//
//bool x_power_sort(std::string a, std::string b) { return x_power(a) < x_power(b); };
//
//bool y_power_sort(std::string a, std::string b) { return y_power(a) < y_power(b); };

#endif // __Params_h
