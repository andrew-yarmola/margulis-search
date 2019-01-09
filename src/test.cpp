#include <stdio.h>
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "ACJ.h"

int main() {
//  SL2<Complex> M = SL2<Complex>(1,Complex(7,5),0,1);
//  printf("%f + i %f\n", (M*M).b.real(), (M*M).b.imag());
//  printf("%f + i %f\n", inverse(M).b.real(), inverse(M).b.imag());
//  ACJ one =  ACJ(1,0,0,0,0);
//  ACJ zero = ACJ(0,0,0,0,0);
//  SL2<ACJ> N = SL2<ACJ>(one,ACJ(XComplex(7,5),0,0,0,0),zero,one);
//  printf("%f + i %f with err %f\n", (N*N).b.f.re, (N*N).b.f.im, (N*N).b.e );
//  printf("%f + i %f with err %f\n", inverse(N).b.f.re, inverse(N).b.f.im, inverse(N).b.e );
//  Params<ACJ> p;
//  p.sinhP = ACJ(1,1,0,1);
//  printf("%f + i %f\n", p.sinhP.f.re, p.sinhP.e);

  // char code[] = "0110111001001001000111001001001001001010110111000110"; 
  char code[] = "011011100100100100011100100100100100101011011100011001"; 
  // char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000"; 
  //char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111110001110001000100010000001100100010";
	Box box;
	for (char* dir = (char *) &code; *dir; ++dir) {
//    printf("%c\n", *dir);
		if (*dir == '0') {
			box = box.child(0);
		} else if (*dir == '1') {
			box = box.child(1);
		}
	}
  const Params<ACJ> params = box.cover();
  ACJ sinhP = params.sinhP;
  ACJ sinhD2 = params.sinhD2;
  ACJ sinhL2 = params.sinhL2; 
  printf("sinhP %f + %f I with size %f, absLB %f, and absUB %f\n", sinhP.f.re, sinhP.f.im, sinhP.size, absLB(sinhP), absUB(sinhP));
  printf("sinhD2 %f + %f I with size %f, absLB %f, and absUB %f\n", sinhD2.f.re, sinhD2.f.im, sinhD2.size, absLB(sinhD2), absUB(sinhD2));
  printf("sinhL2 %f + %f I with size %f, absLB %f, and absUB %f\n", sinhL2.f.re, sinhL2.f.im, sinhL2.size, absLB(sinhL2), absUB(sinhL2));

  SL2<ACJ> x = construct_x(params);
  SL2<ACJ> y = construct_y(params);

  double m_xy_upper = four_cosh_margulis(x,y,true);
  double m_xy_lower = four_cosh_margulis(x,y,false);
  double t_xy_upper = exp_2_t(x,y,true);
  double t_xy_lower = exp_2_t(x,y,false);
//  double m_yx_upper = four_cosh_margulis(y,x,true);
//  double m_yx_lower = four_cosh_margulis(y,x,false);

  printf("Margulis xy between %f and %f\n", m_xy_lower, m_xy_upper);
  printf("Exp 2t xy between %f and %f\n", t_xy_lower, t_xy_upper);
//  printf("Margulis yx between %f and %f\n", m_yx_lower, m_yx_upper);
}
