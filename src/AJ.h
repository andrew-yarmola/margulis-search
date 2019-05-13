#ifndef __AJ_h_
#define __AJ_h_
#include "Complex.h"
#include <assert.h>
#include <stdio.h>
#include "roundoff.h"

struct AJ {
	AJ(const XComplex& f  = 0,
	   const XComplex& z0 = 0,
	   const XComplex& z1 = 0,
	   const XComplex& z2 = 0,
	   const XComplex& w0 = 0,
	   const XComplex& w1 = 0,
	   const XComplex& w2 = 0,
	   double err = 0) : f{f},z0{z0},z1{z1},z2{z2},w0{w0},w1{w1},w2{w2},e{err},
		size((1+3*EPS)*((absUB(z0)+absUB(w0))+((absUB(z1)+absUB(w1))+(absUB(z2)+absUB(w2))))) {}

	XComplex f;
	XComplex z0;
	XComplex z1;
	XComplex z2;
	XComplex w0;
	XComplex w1;
	XComplex w2;
	double e;

	double size;
};

inline const AJ operator-(const AJ& x);
inline const AJ conj(const AJ& x);
inline const AJ operator+(const AJ& x,const AJ& y);
inline const AJ operator-(const AJ& x,const AJ& y);
inline const AJ operator+(const AJ& x,double y);
inline const AJ operator-(const AJ& x,const AJ& y);
inline const AJ operator*(const AJ& x,double y);
inline const AJ operator/(const AJ& x,double y);
inline const double absUB(const AJ& x);
inline const double absLB(const AJ& x);
inline const AJ abs_sqrd(const AJ& x);
inline const AJ abs(const AJ& x);
inline const double size(const AJ& x);
const AJ operator*(const AJ& x,const AJ& y);
const AJ operator/(const AJ& x,const AJ& y);
const AJ operator/(double x,const AJ& y);
const AJ sqrt(const AJ& x);

#include "AJ.inline"
#endif // __AJ_h
