#include "types.hh"
#include "roundoff.h"

double absUB(Complex& x) {
  return (1+EPS)*((1+EPS)*(x.real()*x.real()) + (1+EPS)*(x.imag()*x.imag()));
}

double absLB(Complex& x) {
  return (1-EPS)*((1-EPS)*(x.real()*x.real()) + (1-EPS)*(x.imag()*x.imag()));
}

double max(double a, double b) {
  if (a < b) return b;
  else return a;
}
