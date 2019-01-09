#ifndef __types_h
#define __types_h
#include <complex>

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


typedef std::complex<double> Complex;

double absUB(Complex& x);

double absLB(Complex& x);

double max(double a, double b); 

#endif // __types_h
