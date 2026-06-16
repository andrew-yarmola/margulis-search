#ifndef __Box_h
#define __Box_h
#include "types.hh"
#include "SL2.hh"
#include "AJCC.h"
#include "QuasiRelators.h"
#include <unordered_map>

// DIM: number of real parameters (4 for the symmetric Margulis search).
// The parameter space is ℂ² with coordinates (sinh(D/2), sinh(L/2)), each complex,
// giving 4 real dimensions: (re sinh D/2, re sinh L/2, im sinh D/2, im sinh L/2).
//
// SCL: the initial half-side length of the root box in "digit" units.
// The k-th dimension is scaled by 2^(-k/DIM), so the root box spans:
//   [-SCL·2^(-k/DIM), +SCL·2^(-k/DIM)] in each dimension k.
// This irrational ratio ensures that repeated halving keeps the box "round",
// reducing error accumulation compared to a cubic box. See GMT §7.
#define DIM 4
#define SCL 4

// A Box represents one node in the binary subdivision tree of the parameter space.
// Each box is identified by a binary string `name` (boxcode): '0' means take the
// left (negative) child and '1' the right (positive) child at each level.
//
// The box stores three views of its parameters:
//   center() — exact floating-point center; used for cheap approximate tests.
//   nearer() — for each coordinate, the value closest to 0 inside the box (or 0 if
//              the box straddles the origin); used for monotone-in-|param| bounds.
//   cover()  — AJCC 1-jet sets that rigorously enclose all parameter values in the box;
//              used for verified (proof-generating) elimination tests.
//
// short_words_cache stores computed SL2<AJCC> matrices for words already evaluated
// at this box; it is cleared when descending to a child (child() calls .clear()).
struct Box {
  Box();
  std::string name;
  std::string desc();
  QuasiRelators qr;
  Box child(int dir) const;
  Params<Complex> center() const { return _center; }
  // Returns, for each parameter, the value closest to 0 that lies in the box.
  // Returns 0 when the box straddles the origin (so bounds monotone in |param| work).
  Params<Complex> nearer() const { return _nearer; }
  // Returns AJCC jet-sets that rigorously cover all parameter values in the box.
  Params<AJCC> cover() const { return _cover; }
  SL2<Complex> x_center() const { return _x_center; }
  SL2<Complex> y_center() const { return _y_center; }
  SL2<AJCC> x_cover() const { return _x_cover; }
  SL2<AJCC> y_cover() const { return _y_cover; }
  // Per-box cache of SL2<AJCC> matrices for words; avoids recomputing the same word
  // multiple times during a single box evaluation. Cleared in child().
  std::unordered_map<std::string, SL2<AJCC>> short_words_cache;
  // Per-box memoization of the vol3/sym3 relator sweeps. These tests read only the box
  // cover (NOT the word under evaluation), so they give the same answer for every word
  // tested against this box. Computing them once per box instead of once per word is a
  // large saving on frontier boxes. Reset in child(). See evaluate_vol3/evaluate_sym3.
  bool vol3_done = false;
  bool sym3_done = false;
  TestResult vol3_cached{};
  TestResult sym3_cached{};
private:
  int pos;               // which dimension to split next (cycles 0..DIM-1)
  double center_digits[DIM]; // box center in "digit" units (before scale)
  double size_digits[DIM];   // box half-size in "digit" units
  double box_center[DIM];    // center in true parameter units (with IEEE rounding)
  double box_size[DIM];      // rigorous half-size in true units (includes rounding error)
  Params<Complex> _center;
  Params<Complex> _nearer;
  Params<AJCC> _cover;
  void compute_center_and_size();
  void compute_cover();
  void compute_nearer();
  SL2<Complex> _x_center;
  SL2<Complex> _y_center;
  SL2<AJCC> _x_cover;
  SL2<AJCC> _y_cover;
};

Box get_box(std::string code);

#endif // __Box_h
