#include <stdio.h>
#include "SL2.hh"
#include "Params.hh"
#include "ACJ.h"

int main() {
  SL2<Complex> M = SL2<Complex>(1,Complex(7,5),0,1);
  printf("%f + i %f\n", (M*M).b.real(), (M*M).b.imag());
  printf("%f + i %f\n", inverse(M).b.real(), inverse(M).b.imag());
  Params<ACJ> p;
  p.sinhP = ACJ(1,1,0,1);
  printf("%f + i %f\n", p.sinhP.f.re, p.sinhP.e);
}
