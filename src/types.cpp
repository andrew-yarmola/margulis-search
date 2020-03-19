#include "types.hh"
#include "roundoff.h"

using namespace std;

double pi = atan(1.0) * 4;
double twopi = atan(1.0) * 8;

double absUB(const Complex& x) {
  return (1+2*EPS)*hypot(x.real(), x.imag());
}

double absLB(const Complex& x) {
  return (1-2*EPS)*hypot(x.real(), x.imag());
}

void split_string(const string &str, const string &delims, vector<string> &out) {
  size_t start;
  size_t end = 0;
  while ((start = str.find_first_not_of(delims, end)) != string::npos) {
    end = str.find_first_of(delims, start);
    out.push_back(str.substr(start, end - start));
  }
}

Complex shift_imag_around_zero(const Complex &a) {
  double im = a.imag();
  while (im < -pi) { im += twopi; }
  while (im >  pi) { im -= twopi; }
  return Complex(a.real(), im);
}

Complex shift_imag_around_pi(const Complex &a) {
  double im = a.imag();
  while (im <     0) { im += twopi; }
  while (im > twopi) { im -= twopi; }
  return Complex(a.real(), im);
}

Complex parse_complex(const string &complex_str) {
  vector<string> parts;
  // printf("complex: %s\n", complex_str.c_str());
  split_string(complex_str, "(-+j)", parts);
  if (complex_str.find("-") != string::npos) {
    return Complex(stod(parts[0]), -stod(parts[1]));
  } else {
    return Complex(stod(parts[0]), stod(parts[1]));
  }
}

string repeat(string s, int n) {
  string t = "";
  if (n <= 0) { return t; }
  string r = s;
  while (n > 1) {
    if (n & 1) { // n odd
      t += r;
    }
    r += r;
    n /= 2; // int division
  } 
  return t+r;
}

int x_power(string w) {
  int count = 0;
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'x' || w[p] == 'X') ++count;
  }
  return count;
} 

int y_power(string w) {
  int count = 0;
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'y' || w[p] == 'Y') ++count;
  }
  return count;
} 

bool x_power_sort(string a, string b) { return x_power(a) < x_power(b); }

bool y_power_sort(string a, string b) { return y_power(a) < y_power(b); }

template<>
void print_type<const Complex>(const Complex& x) {
  printf("%f + %f I\n", x.real(), x.imag());
}

template<>
void print_type<Complex>(Complex& x) {
  print_type((const Complex) x);
}

template<>
void print_center<Complex>(const Complex& x) {
  print_type((const Complex) x);
}

template<>
bool comp_type<const Complex>(const Complex& a, const Complex& b) {
  return absUB(a) < absUB(b);
}

template<>
bool comp_type<Complex>(Complex& a, Complex& b) {
  return comp_type((const Complex) a, (const Complex) b);
}
