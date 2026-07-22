#ifndef _SL2_h_
#define _SL2_h_

#include "types.hh"

template<typename T>
struct SL2 {
  SL2<T>() : a(1), b(0), c(0), d(1) {}
  SL2<T>(const T &aa, const T &bb, const T &cc, const T &dd) : a(aa), b(bb), c(cc), d(dd) {}
  T a, b, c, d;
};

template<typename T>
const SL2<T> operator-(const SL2<T> &M, const SL2<T> &N) {
  return SL2<T>(M.a - N.a, M.b - N.b,
                M.c - N.c, M.d - N.d);
};

template<typename T>
const SL2<T> operator*(const SL2<T> &M, const SL2<T> &N) {
  return SL2<T>(M.a*N.a+M.b*N.c, M.a*N.b+M.b*N.d,
                M.c*N.a+M.d*N.c, M.c*N.b+M.d*N.d);
};

template<typename T>
const SL2<T> operator*(const SL2<T> &M, const T &s) {
  return SL2<T>(M.a*s, M.b*s,
                M.c*s, M.d*s);
};

template<typename T>
const SL2<T> inverse(const SL2<T> &M) {
  return SL2<T>(M.d,-M.b,-M.c,M.a);
};

template<typename T>
const T dist(const SL2<T> &M1, const SL2<T> &M2) {
  return abs(M1.a - M2.a) + abs(M1.b - M2.b) + abs(M1.c - M2.c) + abs(M1.d - M2.d); 
};

// v0 = 0
// v1 = 1
// v_(n+1) = x v_n - v_(n_1)
// v_2 = x 1 - 0 = x
// v_3 = x^2 - 1
// v_4 = x(x^2-1) - x = x^3 - 2 x
// v_5 = x (x^3 - 2 x) - (x^2 - 1) = x^4 - 3 x^2 + 1 

// C^n = V_n(tr(C)) * C - V_(n-1)(tr(C)) * I_2
template<typename T>
SL2<T> pow_sl2(const T& t, const SL2<T>& w, int n_pow)
{
  SL2<T> I; // identity is default
  if (n_pow == 0) {
    return I;
  }
  SL2<T> B; // identity is default
  if (n_pow < 0) { 
    B = inverse(w);
    n_pow = -n_pow;
  } else {
    B = w;
  }
  if (n_pow == 1) {
    return B;
  }

  // C^n = V_n(t) C - V_{n-1}(t) I  (Chebyshev acceleration), where V_k are
  // the Chebyshev-like polynomials V_0 = 0, V_1 = 1, V_{k+1} = t V_k - V_{k-1}
  // (V_k(t) = U_{k-1}(t/2), Chebyshev U of the second kind).
  //
  // For 2 <= n <= 10 the pair (V_n, V_{n-1}) is precomputed in closed form and
  // evaluated by Horner in s = t^2. Each V_k has definite parity, so factoring
  // out t and expressing the rest as a polynomial in s uses fewer arithmetic
  // operations than the recurrence and has shallower dependency depth, giving a
  // tighter rigorous enclosure under the AJCC (1+EPS) error model. Higher powers
  // fall back to the O(n) recurrence.
  T v_curr; // V_n     (default-constructed to 0)
  T v_prev; // V_{n-1}
  if (n_pow <= 10) {
    T s = t * t;
    switch (n_pow) {
      case 2: // V_2 = t,                         V_1 = 1
        v_curr = t;
        v_prev = T(1);
        break;
      case 3: // V_3 = s - 1,                     V_2 = t
        v_curr = s - T(1);
        v_prev = t;
        break;
      case 4: // V_4 = t(s - 2),                  V_3 = s - 1
        v_curr = t * (s - T(2));
        v_prev = s - T(1);
        break;
      case 5: // V_5 = (s - 3)s + 1,              V_4 = t(s - 2)
        v_curr = (s - T(3)) * s + T(1);
        v_prev = t * (s - T(2));
        break;
      case 6: // V_6 = t((s - 4)s + 3),           V_5 = (s - 3)s + 1
        v_curr = t * ((s - T(4)) * s + T(3));
        v_prev = (s - T(3)) * s + T(1);
        break;
      case 7: // V_7 = ((s - 5)s + 6)s - 1,       V_6 = t((s - 4)s + 3)
        v_curr = ((s - T(5)) * s + T(6)) * s - T(1);
        v_prev = t * ((s - T(4)) * s + T(3));
        break;
      case 8: // V_8 = t(((s - 6)s + 10)s - 4),   V_7 = ((s - 5)s + 6)s - 1
        v_curr = t * (((s - T(6)) * s + T(10)) * s - T(4));
        v_prev = ((s - T(5)) * s + T(6)) * s - T(1);
        break;
      case 9: // V_9 = (((s - 7)s + 15)s - 10)s + 1
        v_curr = (((s - T(7)) * s + T(15)) * s - T(10)) * s + T(1);
        v_prev = t * (((s - T(6)) * s + T(10)) * s - T(4)); // V_8
        break;
      case 10: // V_10 = t((((s - 8)s + 21)s - 20)s + 5)
        v_curr = t * ((((s - T(8)) * s + T(21)) * s - T(20)) * s + T(5));
        v_prev = (((s - T(7)) * s + T(15)) * s - T(10)) * s + T(1); // V_9
        break;
    }
  } else {
    v_prev = T(0); // V_0
    v_curr = T(1); // V_1
    T v_temp;
    for (int i = 1; i < n_pow; ++i) {
      v_temp = v_curr;
      v_curr = t * v_curr - v_prev;
      v_prev = v_temp;
    }
  }

  return B * v_curr - I * v_prev;
};

template<typename T>
inline const SL2<T> pow_sl2(const SL2<T>& M, int n) {
  return pow_sl2(M.a + M.d, M, n);
};

template<typename T>
inline const SL2<T> pow_old(const SL2<T> &M, int n) {
  SL2<T> A; // identity
  if (n == 0) { return A; }
  SL2<T> B;
  if (n < 0) { 
    B = inverse(M);
    n = -n;
  } else {
    B = M;
  }
  while (n > 1) {
    if (n & 1) { // n odd
      A = A*B;
    }
    B = B*B;
    n /= 2; // int division
  } 
  return A*B;
};

template<typename T>
void print_SL2(const SL2<T>& x) {
  print_type("a =", x.a);
  print_type("b =", x.b);
  print_type("c =", x.c);
  print_type("d =", x.d);
}

template<typename T>
const SL2<T> operator-(const SL2<T> &x) {
  return SL2<T>(-x.a, -x.b, -x.c, -x.d);
}

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
