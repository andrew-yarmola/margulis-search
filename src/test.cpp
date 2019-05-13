#include <stdio.h>
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJ.h"
#include "TubeSearch.hh"
#include "TestCollection.hh"

using namespace std;

int main() {
  SL2<Complex> M = SL2<Complex>(1,Complex(7,5),0,1);
  printf("%f + i %f\n", (M*M).b.real(), (M*M).b.imag());
  printf("%f + i %f\n", inverse(M).b.real(), inverse(M).b.imag());
  AJ one =  AJ(1,0,0,0,0,0,0,0);
  AJ zero = AJ(0,0,0,0,0,0,0,0);
  AJ t = abs(AJ(XComplex(3,4),0,0,0,0,0,0));
  printf("%f + i %f with err %f\n", t.f.re, t.f.im, t.e);
  SL2<AJ> N = SL2<AJ>(one,AJ(XComplex(7,5),0,0,0,0,0,0),zero,one);
  printf("%f + i %f with err %f\n", (N*N).b.f.re, (N*N).b.f.im, (N*N).b.e );
  printf("%f + i %f with err %f\n", inverse(N).b.f.re, inverse(N).b.f.im, inverse(N).b.e );
  Params<AJ> p;
  p.sinhP = AJ(1,1,0,1,0,0,0);
  printf("%f + i %f\n", p.sinhP.f.re, p.sinhP.e);

//  char code[] = "";
  // char code[] = "0110111001001001000111001001001001001010110111000110"; 
//  char code[] = "011011100100100100011100100100100100101011011100011001"; 
//  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000"; 
//  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111";
// char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111110001110001000100010000001100100010";
  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111";
// char code[] = "01101110010010010001110010010010010101100010101011001011010001010111100010100110101100101010011011011000101101010010001001011";

	Box box;
	for (char* dir = (char *) &code; *dir; ++dir) {
//    printf("%c\n", *dir);
		if (*dir == '0') {
			box = box.child(0);
		} else if (*dir == '1') {
			box = box.child(1);
		}
	}

  Params<AJ> params = box.cover();
  AJ sinhP = params.sinhP;
  AJ coshP = params.coshP;

  printf("Box: %s", box.desc().c_str());

  SL2<AJ> x = construct_x(params);
  SL2<AJ> y = construct_y(params);
  SL2<Complex> c_x = construct_x(box.center());
  SL2<Complex> c_y = construct_y(box.center());

  double m_xy_lower = four_cosh_margulis(x,y,false);
  double c_m_xy_lower = four_cosh_margulis(c_x,c_y,false);
  double m_xy_upper = four_cosh_margulis(x,y,true);
  double c_m_xy_upper = four_cosh_margulis(c_x,c_y,true);

  printf("x.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", x.a.f.re, x.a.f.im, x.a.size, absLB(x.a), absUB(x.a));
  printf("y.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", y.a.f.re, y.a.f.im, y.a.size, absLB(y.a), absUB(y.a));
  printf("re_len(x) LB %f\n", four_cosh_re_length_LB(x));
  printf("re_len(y) LB %f\n", four_cosh_re_length_LB(y));
  printf("re_len(x) UB %f\n", four_cosh_re_length_UB(x));
  printf("re_len(y) UB %f\n", four_cosh_re_length_UB(y));
  printf("re_len(p) LB %f\n", absLB(coshP + sinhP));
  printf("re_len(p) UB %f\n", absUB(coshP + sinhP));
  double t_xy_upper = exp_2_t(x,y,true);
  double t_xy_lower = exp_2_t(x,y,false);
  double m_yx_upper = four_cosh_margulis(y,x,true);
  double m_yx_lower = four_cosh_margulis(y,x,false);

  printf("Margulis xy between %f and %f\n", m_xy_lower, m_xy_upper);
  printf("Exp 2t xy between %f and %f\n", t_xy_lower, t_xy_upper);
//  printf("Margulis yx between %f and %f\n", m_yx_lower, m_yx_upper);

  return 1;

  // printf("%s\n",repeat("abc",5).c_str());
  SL2<Complex> f = construct_word("YxYYxxyy", box.center());
  SL2<AJ> g = construct_word("YxYYxxyy", box.cover());
  SL2<AJ> x1 = construct_word("x", box.cover());
  SL2<AJ> y1 = construct_word("y", box.cover());
  SL2<AJ> xy = construct_word("xy", box.cover());
  SL2<AJ> xxxxx = construct_word("xxxxx", box.cover());
  printf("f.a is  %f + %f I\n", f.a.real(), f.a.imag());
  printf("g.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", g.a.f.re, g.a.f.im, g.a.size, absLB(g.a), absUB(g.a));
  printf("x1.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", x1.a.f.re, x1.a.f.im, x1.a.size, absLB(x1.a), absUB(x1.a));
  printf("y1.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", y1.a.f.re, y1.a.f.im, y1.a.size, absLB(y1.a), absUB(y1.a));
  printf("xy.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", xy.a.f.re, xy.a.f.im, xy.a.size, absLB(xy.a), absUB(xy.a));
  printf("(x*y).a is  %f + %f I with size %f, absLB %f, and absUB %f\n", (x*y).a.f.re, (x*y).a.f.im, (x*y).a.size, absLB((x*y).a), absUB((x*y).a));
  printf("xxxxx.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", xxxxx.a.f.re, xxxxx.a.f.im, xxxxx.a.size, absLB(xxxxx.a), absUB(xxxxx.a));
  printf("(x*x*x*x*x).a is  %f + %f I with size %f, absLB %f, and absUB %f\n", (x*x*x*x*x).a.f.re, (x*x*x*x*x).a.f.im, (x*x*x*x*x).a.size, absLB((x*x*x*x*x).a), absUB((x*x*x*x*x).a));

  vector<string> known;
  known.push_back("XXYXXY");
  vector<word_pair> pairs = findPairs(box.center(), known, 50, 20, vector<string>());  
 
  vector<word_pair>::iterator it;
  for (it = pairs.begin(); it != pairs.end(); ++it) {
    printf("(%s, %s)\n", (*it).first.c_str(), (*it).second.c_str());
  }  
}
