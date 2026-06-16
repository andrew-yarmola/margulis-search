#ifndef __Generators_h
#define __Generators_h
#include <math.h>
#include <unordered_map>
#include <string>
#include "SL2.hh"
#include "types.hh"

// cosh of the hyperbolic distance that w moves the basepoint j = (0,0,1) ∈ ℍ³.
// The standard identity is cosh(d(j, w·j)) = (|a|²+|b|²+|c|²+|d|²)/2.
// The code uses the algebraically-equivalent rearrangement (valid because det = 1)
//   (|ac̄+bd̄|² + (|c|²+|d|²-1)²) / (2(|c|²+|d|²)) + 1,
// which one can verify reduces to (|a|²+|b|²+|c|²+|d|²)/2 using |ad-bc|² = 1.
// This form keeps the dominant cancellation inside the (q-1)² term, where q=|c|²+|d|².
// This is the "move" quantity used for the Margulis condition: a word w with
// cosh_move_j(w) < cosh(μ) certifies the Margulis number is smaller than μ.
template<typename T>
const T cosh_move_j(const SL2<T>& w) {
  T q = abs_sqrd(w.c) + abs_sqrd(w.d);
  T z = w.a * conj(w.c) + w.b * conj(w.d);
  return (abs_sqrd(z) + (q - 1) * (q - 1))/(q * 2) + 1;
}

// Generator x: loxodromic with complex half-length sinh(L/2), axis endpoints ±e^{-D/2}.
// In the upper half-space model, axis(x) is the geodesic from -e^{-D/2} to e^{-D/2}.
// The matrix acts as: fixed points ±e^{-D/2}, complex translation length L.
template<typename T>
SL2<T> construct_x(const Params<T>& params) {
	return SL2<T>(params.coshL2, params.expmD2 * params.sinhL2,
                params.expD2 * params.sinhL2, params.coshL2);
};

// Generator y: loxodromic with the same complex half-length sinh(L/2), axis endpoints ±e^{D/2}.
// axis(y) is the geodesic from -e^{D/2} to e^{D/2}. The parameter D is the complex
// distance between axis(x) and axis(y); their common perpendicular runs along the
// {0, ∞} geodesic. (axis(x) and axis(y) are NOT orthogonal in general — they coincide
// when D = 0 and are only perpendicular for special complex D.)
template<typename T>
SL2<T> construct_y(const Params<T>& params) {
	return SL2<T>(params.coshL2, params.expD2 * params.sinhL2,
                params.expmD2 * params.sinhL2, params.coshL2);
};

// Builds the SL2 matrix for a word in the free group on {x,X,y,Y} (uppercase = inverse).
// Parses right-to-left, batching consecutive powers of the same generator before a
// generator-switch so each run x^a or y^b costs one matrix power instead of |a| or |b|
// multiplications. At most one of x_pow or y_pow is nonzero at any given time.
template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params)
{
  SL2<T> w; // identity
  SL2<T> x = construct_x(params);
  SL2<T> y = construct_y(params);

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
    // Flush accumulated power of the opposite generator when we switch.
    if (y_pow != 0 && x_pow != 0) {
      if (h == 'y' || h == 'Y') {
        w = pow(x, x_pow) * w;
        x_pow = 0;
      } else {
        w = pow(y, y_pow) * w;
        y_pow = 0;
      }
    }
  }
  // Flush any remaining leading run (only one of these can be nonzero).
  if (x_pow != 0) { w = pow(x, x_pow) * w; }
  if (y_pow != 0) { w = pow(y, y_pow) * w; }
  return w;
};

// Cache-aware version: returns the cached matrix if present, otherwise computes
// it, stores the result, and returns it. The cache is keyed by word string and
// is valid only for a single box (it is cleared in Box::child()).
template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params,
                      std::unordered_map<std::string, SL2<T>>& cache) {
  auto it = cache.find(word);
  if (it != cache.end()) {
    return it->second;
  }
  SL2<T> result = construct_word(word, params);
  cache.emplace(word, result);
  return result;
}


#endif // __Generators_h
