#ifndef __Params_h
#define __Params_h
#include <math.h>
#include <unordered_map>
#include <string>
#include "SL2.hh"
#include "types.hh"

/*
* Let x and y realize the margulis constant such that re(length(x)) =< re(length(y))
* We will conjugate x to have axis(x) = {0,\infty}
*/

template<typename T>
SL2<T> construct_x(const Params<T>& params) {
  // We want the matrix {{ exp(L/2), 0 }, { 0, exp(-L/2) }}
  T zero = T(0);
  T sl2 = params.sinhL2;
  T cl2 = params.coshL2;
	return SL2<T>(cl2 + sl2, zero,
                zero, cl2 - sl2);
};

template<typename T>
SL2<T> construct_y(const Params<T>& params) {
  // We want the matrix {{ cosh(D/2) + cosh(P)sinh(D/2), -sinh(D/2)sinh(P) }, 
  //                     { sinh(D/2)sinh(P), cosh(D/2) - cosh(P)sinh(D/2) }}
  T sd2 = params.sinhD2; 
  T cd2 = params.coshD2;
  T sp = params.sinhP;
  T cp = params.coshP;
	return SL2<T>(cd2 + cp * sd2, -(sp * sd2),
                sp * sd2      , cd2 - cp * sd2);
};

template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params)
{
  // TODO : add caching of powers if we are too slow
  std::string recover = "";
  SL2<T> w; // identity
	SL2<T> x(construct_x(params));
	SL2<T> y(construct_y(params));

  char h;
  int x_pow = 0;
  int y_pow = 0;	
  std::string::reverse_iterator rit;
  for (rit = word.rbegin(); rit != word.rend(); ++rit) {
    h = *rit;
	  switch(h) {
      case 'x': ++x_pow; break;
      case 'X': --x_pow; break;
      case 'y': ++y_pow; break;
      case 'Y': --y_pow; break;
    }
    if (y_pow != 0 && x_pow != 0) {
      if (h == 'y' || h == 'Y' ) {
        w = pow(x, x_pow) * w;
        if (x_pow < 0) {
          recover = repeat("X",-x_pow) + recover;
        } else {
          recover = repeat("x",x_pow) + recover;
        }
        x_pow = 0;
      } else {
        w = pow(y, y_pow) * w;
        if (y_pow < 0) {
          recover = repeat("Y",-y_pow) + recover;
        } else {
          recover = repeat("y",y_pow) + recover;
        }
        y_pow = 0;
      }
    }
  }
  // Only one of these should be true
  if (x_pow != 0) { 
    w = pow(x, x_pow) * w;
    if (x_pow < 0) {
      recover = repeat("X",-x_pow) + recover;
    } else {
      recover = repeat("x",x_pow) + recover;
    }
  }
  if (y_pow != 0) { 
    w = pow(y, y_pow) * w;
    if (y_pow < 0) {
      recover = repeat("Y",-y_pow) + recover;
    } else {
      recover = repeat("y",y_pow) + recover;
    }
  }
  // fprintf(stderr, "Original %s vs %s\n", word.c_str(), recover.c_str());
	return w;
};

template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params, std::unordered_map< std::string, SL2<T> >& word_cache) {
  // TODO: fill and make use of cache
  return construct_word(word, params);  
}


#endif // __Params_h
