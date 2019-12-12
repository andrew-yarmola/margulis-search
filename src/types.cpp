#include "types.hh"
#include "roundoff.h"

double absUB(const Complex& x) {
  return (1+2*EPS)*hypot(x.real(), x.imag());
}

double absLB(const Complex& x) {
  return (1-2*EPS)*hypot(x.real(), x.imag());
}

double max(double a, double b) {
  if (a < b) return b;
  else return a;
}

std::string repeat(std::string s, int n) {
  std::string t = "";
  if (n <= 0) { return t; }
  std::string r = s;
  while (n > 1) {
    if (n & 1) { // n odd
      t += r;
    }
    r += r;
    n /= 2; // int division
  } 
  return t+r;
}

int x_power(std::string w) {
  int count = 0;
  for (std::string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'x' || w[p] == 'X') ++count;
  }
  return count;
} 

int y_power(std::string w) {
  int count = 0;
  for (std::string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'y' || w[p] == 'Y') ++count;
  }
  return count;
} 

bool x_power_sort(std::string a, std::string b) { return x_power(a) < x_power(b); }

bool y_power_sort(std::string a, std::string b) { return y_power(a) < y_power(b); }

template<>
void print_type<const Complex>(const Complex& x) {
  printf("%f + %f I\n", x.real(), x.imag());
}

template<>
void print_type<Complex>(Complex& x) {
  print_type((const Complex) x);
}

template<>
bool comp_type<const Complex>(const Complex& a, const Complex& b) {
  return std::abs(a) < std::abs(b);
}

template<>
bool comp_type<Complex>(Complex& a, Complex& b) {
  return comp_type((const Complex) a, (const Complex) b);
}
