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
	T coshD2;
	T sinhL2;
	T coshL2;
};

template<typename T>
const SL2<T> construct_x(const Params<T>& params);

template<typename T>
const SL2<T> construct_y(const Params<T>& params);

int x_power(std::string w);
bool x_power_sort(std::string a, std::string b);

int y_power(std::string w);
bool y_power_sort(std::string a, std::string b);

#endif // __Params_h
