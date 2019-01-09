#ifndef _SL2_h_
#define _SL2_h_

template<typename T>
struct SL2 {
  SL2<T>() : a(1), b(0), c(0), d(1) {}
  SL2<T>(const T &aa, const T &bb, const T &cc, const T &dd) : a(aa), b(bb), c(cc), d(dd) {}
  T a, b, c, d;
};

template<typename T>
const SL2<T> operator*(const SL2<T> &x, const SL2<T> &y) {
  return SL2<T>(x.a*y.a+x.b*y.c, x.a*y.b+x.b*y.d,
                x.c*y.a+x.d*y.c, x.c*y.b+x.d*y.d);
};

template<typename T>
const SL2<T> inverse(const SL2<T> &x) {
  return SL2<T>(x.d,-x.b,-x.c,x.a);
};

#endif

/* Don't need and not in SL2 
template<typename T>
const SL2<T> operator+(const SL2<T> &x, const SL2<T> &y) {
  return SL2<T>(x.a+y.a, x.b+y.b, x.c+y.c, x.d+y.d);
}


template<typename T>
const SL2<T> operator-(const SL2<T> &x, const SL2<T> &y) {
  return SL2<T>(x.a-y.a, x.b-y.b, x.c-y.c, x.d-y.d);
} */
