#include "AJ.h"

const AJ operator*(const AJ&x,const AJ&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.a), ay = absUB(y.a);
	AComplex r_a = x.a*y.a;
	AComplex r_z0 = x.a*y.z0+x.z0*y.a;
	AComplex r_z1 = x.a*y.z1+x.z1*y.a;
	AComplex r_z2 = x.a*y.z2+x.z2*y.a;
	AComplex r_w0 = x.a*y.w0+x.w0*y.a;
	AComplex r_w1 = x.a*y.w1+x.w1*y.a;
	AComplex r_w2 = x.a*y.w2+x.w2*y.a;
	double A = (xdist+x.e)*(ydist+y.e);
	double B = ax*y.e+ay*x.e;
	double C = (r_a.e+(r_z0.e+r_w0.e))+((r_z1.e+r_w1.e)+(r_z2.e+r_w2.e));
	double r_error = (1+3*EPS)*((A+B)+C));
	return AJ(r_a.z,r_z0.z,r_z1.z,r_z2.z,r_w0.z,r_w1.z,r_w2.z,r_error);

}
const AJ operator/(const AJ&x,const AJ&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.a), ay = absLB(y.a);
	double D = ay-(1+EPS)*(y.e+ydist);
	if(!(D > 0))return AJ(0,0,0,0,0,0,0,infinity());
	AComplex den = (y.a*y.a);
	AComplex r_a = x.a/y.a;
	AComplex r_z0 = (x.z0*y.a-x.a*y.z0)/den;
	AComplex r_z1 = (x.z1*y.a-x.a*y.z1)/den;
	AComplex r_z2 = (x.z2*y.a-x.a*y.z2)/den;
	AComplex r_w0 = (x.w0*y.a-x.a*y.w0)/den;
	AComplex r_w1 = (x.w1*y.a-x.a*y.w1)/den;
	AComplex r_w2 = (x.w2*y.a-x.a*y.w2)/den;
	double A = (ax+(xdist+x.e))/D;
	double B = (ax/ay+xdist/ay)+(ydist*ax)/(ay*ay);
	double C = (r_a.e+(r_z0.e+r_w0.e))+((r_z1.e+r_w1.e)+(r_z2.e+r_w2.e));
	double r_error = (1+3*EPS)*(((1+3*EPS)*A-(1-3*EPS)*B)+C);
	return AJ(r_a.z,r_z0.z,r_z1.z,r_z2.z,r_w0.z,r_w1.z,r_w2.z,r_error);

}
const AJ operator/(double x,const AJ&y) {

	double ydist = size(y);
	double ax = fabs(x), ay = absLB(y.a);
	double D = ay-(1+EPS)*(y.e+ydist);
	if(!(D > 0))return AJ(0,0,0,0,0,0,0,infinity());
	AComplex den = (y.a*y.a);
	AComplex r_a = x/y.a;
	AComplex r_z0 = (-x*y.z0)/den;
	AComplex r_z1 = (-x*y.z1)/den;
	AComplex r_z2 = (-x*y.z2)/den;
	AComplex r_w0 = (-x*y.w0)/den;
	AComplex r_w1 = (-x*y.w1)/den;
	AComplex r_w2 = (-x*y.w2)/den;
	double B = ax/ay+(ydist*ax)/(ay*ay);
	double C = (r_a.e+(r_z0.e+r_w0.e))+((r_z1.e+r_w1.e)+(r_z2.e+r_w2.e));
	double r_error = (1+3*EPS)*(((1+2*EPS)*(ax/D)-(1-3*EPS)*B)+C);
	return AJ(r_a.z,r_z0.z,r_z1.z,r_z2.z,r_w0.z,r_w1.z,r_w2.z,r_error);

}
const AJ sqrt(const AJ&x) {

	double xdist = size(x);
	double ax = absUB(x.a);
	double D = ax-(1+EPS)*(xdist+x.e);
	if(!(D > 0)) {
		return AJ(0,0,0,0,0,0,0,(1+2*EPS)*sqrt(ax+(xdist+x.e)));
	} else {
		AComplex r_a = sqrt(x.a);
		AComplex t = r_a+r_a;
		AComplex r_z0 = AComplex(x.z0.re,x.z0.im,0)/t;
		AComplex r_z1 = AComplex(x.z1.re,x.z1.im,0)/t;
		AComplex r_z2 = AComplex(x.z2.re,x.z2.im,0)/t;
		AComplex r_w0 = AComplex(x.w0.re,x.w0.im,0)/t;
		AComplex r_w1 = AComplex(x.w1.re,x.w1.im,0)/t;
		AComplex r_w2 = AComplex(x.w2.re,x.w2.im,0)/t;
	  double C = (r_a.e+(r_z0.e+r_w0.e))+((r_z1.e+r_w1.e)+(r_z2.e+r_w2.e));
		double r_error = (1+3*EPS)*(
							((1+EPS)*sqrt(ax)-(1-3*EPS)*(xdist/(2*sqrt(ax))+sqrt(D)))+C);
	  return AJ(r_a.z,r_z0.z,r_z1.z,r_z2.z,r_w0.z,r_w1.z,r_w2.z,r_error);
	}

}
