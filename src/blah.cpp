#include <stdio.h>
#include "SL2.hh"

int expo(int a, int n){
  int x = 1;
  int y = a;
  while (n > 1) {
    if (n & 1){
      x *= y;
    }
    y *= y; 
    n >>=1;
  }
  return x*y;
}

int main() {
  SL2<int> A(1,2,0,1);
  SL2<int> B = pow(A,6);
  printf("A = (%d,%d,%d,%d) and A^6 = (%d,%d,%d,%d)\n",A.a,A.b,A.c,A.d,B.a,B.b,B.c,B.d);
}
