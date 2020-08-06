#ifndef __types_h
#define __types_h
#include <complex>
#include <vector>
#include <string.h>
#include <algorithm>

typedef std::complex<double> Complex;
typedef std::pair<double, double> float_pair;
typedef std::pair<std::string, std::string> word_pair;

inline const Complex eye(const Complex&x) { return Complex(0,1); };
inline const Complex operator+(const Complex&x,double y) { return x + Complex(y,0); };
inline const Complex operator-(const Complex&x,double y) { return x - Complex(y,0); };
inline const Complex operator*(const Complex&x,double y) { return x * Complex(y,0); };
inline const Complex operator/(const Complex&x,double y) { return x / Complex(y,0); };
inline const Complex abs(const Complex& x) { return Complex(hypot(x.real(), x.imag()), 0); };
inline const Complex abs_sqrd(const Complex& x) { return Complex(x.real()*x.real()+x.imag()*x.imag(), 0); };

Complex parse_complex(const std::string &complex_str);
Complex shift_imag_around_zero(const Complex &a);
Complex shift_imag_around_pi(const Complex &a);

double absUB(const Complex& x);
double absLB(const Complex& x);

std::string repeat(std::string s, int n);
void split_string(const std::string &str, const std::string &delims, std::vector<std::string> &out);

int x_power(std::string w);
int y_power(std::string w);

bool x_power_sort(std::string a, std::string b);
bool y_power_sort(std::string a, std::string b);

template<typename T>
inline T max_local(const T& x, const T& y) {
  return ((x+y)+abs(x-y))/2;
} 

template<typename T>
inline T min_local(const T& x, const T& y) {
  return ((x+y)-abs(x-y))/2;
} 

/*
 * We will work with 6 real parameters and several derived parameters.
 * coshmu = cosh of margulis constat
 * sinhdx = sinh(d_x) where d_x is the hyperbolic distance from the margulis point to axis(x)
 * sinhdy = sinh(d_y) where d_y is the hyperbolic distance from the margulis point to axis(y)
 * cosf = cos(phi) where phi is the twist angle from the oriented axis(x) to axis(y)
note: we can assume 0 <= phi <= pi by using a reflection symmtery
 * sintx2 = sin(theta_x/2) where -pi <= theta_x <= pi is the loxodromic rotation angle of x
 * sinty2 = sin(theta_y/2) where -pi <= theta_y <= pi is the loxodromic rotation angle of y
 * 
 * The derived paramters are
 * coshdx = sqrt(1+sinh(d_x)^2)
 * coshdx = sqrt(1+sinh(d_y)^2)
 * sinf = sqrt(1-cos(phi)^2)
 * costx2 = sqrt(1-sin(t_x/2)^2)
 * costy2 = sqrt(1-sin(t_y/2)^2)
 *
 * costx = 1 - 2*sin(t_x/2)^2
 * costy = 1 - 2*sin(t_y/2)^2
 *
 * coshlx = (cosh(mu) + cos(tx)*sinh(dx)^2)/cosh(d2)^2
 * coshly = (cosh(mu) + cos(ty)*sinh(dy)^2)/cosh(d2)^2
 *
 * coshlx2 = sqrt((coshlx + 1)/2)
 * sinhlx2 = sqrt((coshlx - 1)/2)
 * coshly2 = sqrt((coshly + 1)/2)
 * sinhly2 = sqrt((coshly - 1)/2)
 *
 * With L = l_x + i theta_x and D = l_y + i theta_y
 * coshLx2 = cosh(L_x/2) = cosh(l_x/2)cos(theta_x/2) + i sinh(l_x/2)sin(theta_x/2)
 * sinhLx2 = sinh(L_x/2) = sinh(l_x/2)cos(theta_x/2) + i cosh(l_x/2)sin(theta_x/2)
 * coshLy2 = cosh(L_y/2) = cosh(l_y/2)cos(theta_y/2) + i sinh(l_y/2)sin(theta_y/2)
 * sinhLy2 = sinh(L_y/2) = sinh(l_y/2)cos(theta_y/2) + i cosh(l_y/2)sin(theta_y/2)
 * 
 * expdx = coshdx + sinhdx
 * expmdx = coshdx - sinhdx
 * expdy = coshdy + sinhdy
 * expmdy = coshdy - sinhdy
 * 
 * expif = cosf + i sinf
 * expmif = cosf - i sinf
 */

template<typename T> struct Params {
  T coshmu;
  T sinhdx;
  T sinhdy;
  T cosf;
  T sintx2;
  T sinty2;
  T sinhsdx; // derived parameter
  T sinhsdy; // derived parameter
  T coshsdx; // derived parameter
  T coshsdy; // derived parameter
  T coshdx; // derived parameter
  T coshdy; // derived parameter
  T sinf; // derived parameter
  T costx2; // derived parameter
  T costy2; // derived parameter
  T costx; // derived parameter
  T costy; // derived parameter
  T coshlx; // derived parameter
  T coshly; // derived parameter
  T cosh2dx; // derived parameter
  T cosh2dy; // derived parameter
  T coshdxdy; // derived parameter
  T sinhdxdy; // derived parameter
  T coshlx2; // derived parameter
  T sinhlx2; // derived parameter
  T coshly2; // derived parameter
  T sinhly2; // derived parameter
  T coshLx2; // derived parameter
  T sinhLx2; // derived parameter
  T coshLy2; // derived parameter
  T sinhLy2; // derived parameter
  T sinhperp; //derived parameter
  T expdx; // derived parameter
  T expmdx; // derived parameter
  T expdy; // derived parameter
  T expmdy; // derived parameter
  T expif; // derived parameter
  T expmif; // derived parameter
  T expdyf; // derived parameter
  T expmdyf; // derived parameter
  T explx; // derived parameter
  T expmlx; // derived parameter
  T explx2; // derived parameter
  T expmlx2; // derived parameter
  T exply; // derived parameter
  T expmly; // derived parameter
  T exply2; // derived parameter
  T expmly2; // derived parameter
};

template<typename T>
void fill_derived(Params<T>& p) {
  T one = T(1);
  T i = eye(one); // HACK

  p.sinhsdx = p.sinhdx * p.sinhdx;
  p.sinhsdy = p.sinhdy * p.sinhdy;
  p.coshsdx = p.sinhdx * p.sinhdx + one;
  p.coshsdy = p.sinhdy * p.sinhdy + one;
  p.coshdx = sqrt(p.coshsdx);  
  p.coshdy = sqrt(p.coshsdy);

  p.cosh2dx = p.sinhsdx * 2 + one;
  p.cosh2dy = p.sinhsdy * 2 + one;

  p.coshdxdy = p.coshdx * p.coshdy + p.sinhdx * p.sinhdy;
  p.sinhdxdy = p.sinhdx * p.coshdy + p.coshdx * p.sinhdy;

  p.sinf = sqrt(one - p.cosf * p.cosf);

  p.costx2 = sqrt(one - p.sintx2 * p.sintx2); 
  p.costy2 = sqrt(one - p.sinty2 * p.sinty2); 

  p.costx = one - (p.sintx2 * p.sintx2) * 2; 
  p.costy = one - (p.sinty2 * p.sinty2) * 2;

  // main formula relating margulis, real length, twist, and distance to axis
  p.coshlx = (p.coshmu + p.costx * p.sinhsdx) / p.coshsdx;
  p.coshly = (p.coshmu + p.costy * p.sinhsdy) / p.coshsdy;

  p.coshlx2 = sqrt((p.coshlx + one) / 2); 
  p.sinhlx2 = sqrt((p.coshlx - one) / 2);
  p.coshly2 = sqrt((p.coshly + one) / 2); 
  p.sinhly2 = sqrt((p.coshly - one) / 2);

  p.coshLx2 = p.coshlx2 * p.costx2 + i * (p.sinhlx2 * p.sintx2); 
  p.sinhLx2 = p.sinhlx2 * p.costx2 + i * (p.coshlx2 * p.sintx2); 
  p.coshLy2 = p.coshly2 * p.costy2 + i * (p.sinhly2 * p.sinty2); 
  p.sinhLy2 = p.sinhly2 * p.costy2 + i * (p.coshly2 * p.sinty2); 

  p.sinhperp = p.sinhdxdy * p.cosf + (p.coshdxdy * p.sinf) * i;

  p.expdx = p.coshdx + p.sinhdx; 
  p.expmdx = p.coshdx - p.sinhdx;
  p.expdy = p.coshdy + p.sinhdy; 
  p.expmdy = p.coshdy - p.sinhdy;

  p.expif = p.cosf + i * p.sinf;
  p.expmif = p.cosf - i * p.sinf;

  p.expdyf = p.expdy * p.expif;
  p.expmdyf = p.expmdy * p.expmif;

  p.explx = p.coshlx + sqrt(p.coshlx * p.coshlx - one);
  p.expmlx = p.coshlx - sqrt(p.coshlx * p.coshlx - one);
  p.explx2 = p.coshlx2 + p.sinhlx2;
  p.expmlx2 = p.coshlx2 - p.sinhlx2;

  p.exply = p.coshly + sqrt(p.coshly * p.coshly - one);
  p.expmly = p.coshly - sqrt(p.coshly * p.coshly - one);
  p.exply2 = p.coshly2 + p.sinhly2;
  p.expmly2 = p.coshly2 - p.sinhly2;
}

inline double re_center(const Complex& x) {
  return x.real();
}

template<typename T>
void print_type(T& x); 

template<typename T>
void print_center(const T& x); 

template<typename T>
bool sort_comp(const T& a, const T& b); 

template<typename T>
inline bool strictly_pos(const T& diff) {
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
void print_type(const char desc[], const T& x) {
  printf("%s\n", desc);
  print_type(x);
} 

template<typename T>
void print_center(const char desc[], const T& x) {
  printf("%s\n", desc);
  print_center(x);
} 

#endif // __types_h
