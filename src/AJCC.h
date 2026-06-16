#ifndef __AJCC_h_
#define __AJCC_h_
#include "Complex.h"
#include <assert.h>
#include <stdio.h>
#include "roundoff.h"

// ℓ¹-norm of a complex number, used to bound the gradient magnitude. The true
// absolute value costs sqrt, so we use |re|+|im| ≥ |z| as a cheap upper bound.
// The caller must account for the factor of (1+EPS/2) from rounding.
inline double abs_taxi(const XComplex& x) {
  return fabs(x.re) + fabs(x.im);
}

// AJCC — Affine 1-Jet Complex Complex: a rigorously verified set of complex-valued
// functions over the 2-complex-parameter box { (ζ₀,ζ₁) ∈ ℂ² : |ζᵢ| ≤ 1 }.
//
// The box coordinates ζ₀ and ζ₁ are normalized so that the two primary parameters
// satisfy:  sinh(L/2) = f_sinhL2 + z0·ζ₀   and   sinh(D/2) = f_sinhD2 + z1·ζ₁.
// (The antiholomorphic fields w0,w1 arise from operations like conj() and abs().)
//
// An AJCC {f, z0, z1, w0, w1, e} represents the set of functions g: ℂ² → ℂ such that
//
//   |g(ζ₀,ζ₁) − f − z0·ζ₀ − z1·ζ₁ − w0·ζ̄₀ − w1·ζ̄₁| ≤ e   for all |ζᵢ| ≤ 1.
//
// Fields:
//   f    — value at the box center
//   z0   — holomorphic partial w.r.t. the sinh(L/2) coordinate
//   z1   — holomorphic partial w.r.t. the sinh(D/2) coordinate
//   w0   — antiholomorphic partial w.r.t. sinh(L/2) (nonzero after conj/abs)
//   w1   — antiholomorphic partial w.r.t. sinh(D/2)
//   e    — total rounding + approximation error across the whole box
//   size — precomputed ℓ¹ upper bound on |z0|+|z1|+|w0|+|w1| (gradient magnitude)
//
// absUB(x) = (1+2ε)(|f| + size + e) is a rigorous upper bound on |g| over the box.
// absLB(x) = max(0, (1-ε)(|f| - (1+ε)(size + e))) is a rigorous lower bound.
//
// See GMT Annals of Math. 157 (2003) §7-8 for the arithmetic derivations.
struct AJCC {
	AJCC(const XComplex& f  = 0,
	   const XComplex& z0 = 0,
	   const XComplex& z1 = 0,
	   const XComplex& w0 = 0,
	   const XComplex& w1 = 0,
	   double err = 0) : f{f},z0{z0},z1{z1},w0{w0},w1{w1},e{err},
		size((1+3*EPS) * ((abs_taxi(z0) + abs_taxi(z1)) +
                          (abs_taxi(w0) + abs_taxi(w1)))) {}
	XComplex f;   // center value
	XComplex z0;  // holomorphic partial w.r.t. sinh(L/2)
	XComplex z1;  // holomorphic partial w.r.t. sinh(D/2)
	XComplex w0;  // antiholomorphic partial w.r.t. sinh(L/2)
	XComplex w1;  // antiholomorphic partial w.r.t. sinh(D/2)
	double e;     // error bound: |g − affine approximation| ≤ e over the box
	double size;  // ℓ¹ upper bound on gradient magnitude, precomputed for speed
};

inline const AJCC eye(const AJCC& x) { return AJCC(XComplex(0,1)); };
inline const AJCC operator-(const AJCC& x);
inline const AJCC conj(const AJCC& x);
inline const AJCC re(const AJCC& x);
inline const AJCC im(const AJCC& x);
inline const AJCC operator+(const AJCC& x,const AJCC& y);
inline const AJCC operator-(const AJCC& x,const AJCC& y);
inline const AJCC operator+(const AJCC& x,double y);
inline const AJCC operator-(const AJCC& x,double y);
inline const AJCC operator*(const AJCC& x,double y);
inline const AJCC operator/(const AJCC& x,double y);
inline const double absUB(const AJCC& x);
inline const double absLB(const AJCC& x);
inline const AJCC abs_sqrd(const AJCC& x);
inline const AJCC abs(const AJCC& x);
inline const double size(const AJCC& x);
const AJCC operator*(const AJCC& x,const AJCC& y);
const AJCC operator/(const AJCC& x,const AJCC& y);
const AJCC operator/(double x,const AJCC& y);
const AJCC sqrt(const AJCC& x);

#include "AJCC.inline"
#endif // __AJCC_h
