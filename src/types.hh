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
	T sinhD2;
	T sinhL2;
	T coshP; // may not always be filled
	T coshD2; // may not always be filled
	T coshL2; // may not always be filled
};


typedef std::complex<double> Complex;

inline const Complex operator+(const Complex&x,double y) { return x + Complex(y,0); };
inline const Complex operator-(const Complex&x,double y) { return x - Complex(y,0); };
inline const Complex operator*(const Complex&x,double y) { return x * Complex(y,0); };
inline const Complex operator/(const Complex&x,double y) { return x / Complex(y,0); };

double absUB(const Complex& x);

double absLB(const Complex& x);

double max(double a, double b); 

std::string repeat(std::string s, int n);

int x_power(std::string w);

int y_power(std::string w);

bool x_power_sort(std::string a, std::string b);

bool y_power_sort(std::string a, std::string b);

#endif // __types_h
