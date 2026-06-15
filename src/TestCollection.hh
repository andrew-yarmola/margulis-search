#ifndef _TestCollection_ 
#define _TestCollection_
#include <unordered_map>
#include <string>
#include <set>
#include <vector>
#include "types.hh"
#include "Box.h"
#include "SL2.hh"
#include "IsomH3.hh"
#include "Params.hh"
#include "RelatorTest.hh"

extern bool g_debug;

struct RelatorTest;

struct TestCollection {
  int size();
  box_state evaluate_center(int index, Box& box);
  TestResult evaluate_box(int index, Box& box);
  TestResult evaluate_AJCC(word_pair& pair, Box& box);
  TestResult evaluate_qrs(Box& box);
  TestResult evaluate_vol3(Box& box);
  TestResult evaluate_sym3(Box& box);
  const std::string get_name(int index);
  word_pair get_pair(int index);
  int add(word_pair pair);
  int add(std::string pair);
  void load(const char* file_path);
  void load_relator_test(const char* impos_path,
                         const char* bad_rel_path);
  RelatorTest *relator_test;
  std::map<std::string, int> seen_words;
private:
  word_pair parse_word_pair(std::string buf);
  std::map<word_pair, int> pair_index;
  std::vector<word_pair> pair_vector;
  box_state evaluate_approx(word_pair pair, const Box& params);
  bool ready_for_elliptics_test(SL2<AJCC>& w);
  bool only_elliptics(SL2<AJCC>& w, Params<AJCC>& params);
};

// Jørgensen's inequality: |tr(w1)²-4| + |tr([w1,w2])-2| ≥ 1 for any discrete
// non-elementary group. If this sum is provably < 1 and at least one of w1,w2
// is not parabolic, then in a discrete group they must commute (variety neighbourhood).
template<typename T>
inline const bool inside_var_nbd(const SL2<T>& w1, const SL2<T>& w2) {
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_parabolic(w1) || not_parabolic(w2));
}

// Stricter form: both words are loxodromic, so inside the variety neighbourhood they
// commute and share an axis — proving a relator exists between them.
template<typename T>
inline const bool inside_var_nbd_ne(const SL2<T>& w1, const SL2<T>& w2) {
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_elliptic_or_parabolic(w1) && not_elliptic_or_parabolic(w2));
}

// A Möbius transformation is parabolic iff tr = ±2 (real). Certifies tr ≠ ±2.
template<typename T>
inline const bool not_parabolic(const SL2<T>& w) {
  T tr = w.a + w.d;
  return absLB(im(tr)) > 0 || (absLB(re(tr) - 2) > 0 && absLB(re(tr) + 2) > 0);
}

// Certifies that w is loxodromic: tr ∉ [-2,2] ⊂ ℝ (i.e., not elliptic or parabolic).
template<typename T>
inline const bool not_elliptic_or_parabolic(const SL2<T>& w) {
  T tr = w.a + w.d;
  return absLB(im(tr)) > 0 || absLB(re(tr)) > 2;
}

// Certifies w ≠ ±I by showing at least one off-diagonal entry is nonzero, or that
// the diagonal entries are not ±1.
template<typename T>
inline const bool not_identity(const SL2<T>& w) {
  return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
    ((absLB(w.a-1) > 0 || absLB(w.d-1) > 0) && (absLB(w.a+1) > 0 || absLB(w.d+1) > 0));
}

template<typename T>
inline const bool tube_hits_axis_two(const T& two_sinh_p2sq, const T& two_cosh_re_tube) {
  T tcd = two_cosh_dist(two_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube - tcd) && absLB(two_sinh_p2sq + 1) > 0;
}

template<typename T>
inline const bool tube_hits_axis_four(const T& four_sinh_p2sq, const T& two_cosh_re_tube) {
  T fcd = four_cosh_dist(four_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube * 2 - fcd) && absLB(four_sinh_p2sq + 2) > 0;
}

template<typename T>
inline const bool wx_hits_sym_axis(const SL2<T>& w, const Params<T>& p) {
  T tsp2sq_inf = two_sinh_perp2_sq_wax_zero_inf(w, p);
  T fsp2sq_one = four_sinh_perp2_sq_wax_mp_one(w, p);
  T fsp2sq_iye = four_sinh_perp2_sq_wax_mp_iye(w, p);
  if (g_debug) {
    if (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2)) {
      fprintf(stderr, "Hit zero inf\n");
      print_type(tsp2sq_inf);
    }
    if (tube_hits_axis_four(fsp2sq_one, p.twocoshreD2)) {
      fprintf(stderr, "Hit +- one\n");
      print_type(fsp2sq_one);
    }
    if (tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2)){
      fprintf(stderr, "Hit +- iye\n");
      print_type(fsp2sq_iye);
    }
  }
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
          tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
          tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

template<typename T>
inline const bool wy_hits_sym_axis(const SL2<T>& w, const Params<T>& p) {
  T tsp2sq_inf = two_sinh_perp2_sq_way_zero_inf(w, p);
  T fsp2sq_one = four_sinh_perp2_sq_way_mp_one(w, p);
  T fsp2sq_iye = four_sinh_perp2_sq_way_mp_iye(w, p);
  if (g_debug) {
    if (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2)) {
      fprintf(stderr, "Hit zero inf\n");
      print_type(tsp2sq_inf);
    }
    if (tube_hits_axis_four(fsp2sq_one, p.twocoshreD2)) {
      fprintf(stderr, "Hit +- one\n");
      print_type(fsp2sq_one);
    }
    if (tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2)){
      fprintf(stderr, "Hit +- iye\n");
      print_type(fsp2sq_iye);
    }
  }
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
          tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
          tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

template<typename T>
inline const bool does_not_fix_sym_axis(const SL2<T>& w) {
  T one(1);
  T zero(0);
  T iye(0,1);
  return (absLB(w.b) > 0 && absLB(w.d) > 0) 
    || (absLB(w.a) > 0 && absLB(w.c) > 0) ||
    (absLB(mobius(w, one) - one) > 0 && absLB(mobius(w, one) + one) > 0) || 
    (absLB(mobius(w, iye) - iye) > 0 && absLB(mobius(w, iye) + iye) > 0) || 
    (absLB(mobius(w, -one) - one) > 0 && absLB(mobius(w, -one) + one) > 0) || 
    (absLB(mobius(w, -iye) - iye) > 0 && absLB(mobius(w, -iye) + iye) > 0); 
}

#define LERR 0.00000000001

template<typename T>
inline const bool really_cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  if (g_debug && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
    fprintf(stderr, "Realy can't fix x_axis");
    print_type(fsp2sq);
    fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  return absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR;
}

template<typename T>
inline const bool cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
    fprintf(stderr, "********** CANNOT FIX X AXIS ***********\n");
    T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
    print_type(fsp2sq);
    fprintf(stderr, "Can't fix x axis LB values %f and %f\n",
      absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr,"Can't fix x axis LB away from %d and %d\n",
      absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
    fprintf(stderr, "*******************************\n");
  }
  // return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  return absLB(fsp2sq) > 0; 
}

// Returns true if w maps axis(x) to an axis that is provably closer to axis(x) than
// the tube radius 2·Re(D/2). In a discrete group this forces w to fix axis(x) (as an
// oriented geodesic), because the tube about axis(x) would otherwise be violated.
// Only valid for AJCC (verified) arithmetic.
template<typename T>
inline const bool must_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshreD * 4 - four_cosh_dist_ax_wax(w, p);
  if (g_debug && 
      //std::is_same<T, AJCC>::value && strictly_pos(diff)) {
      std::is_same<T, AJCC>::value) {
    fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
    print_SL2(w);
    print_type("4 cosh 2 dx:", p.coshreD * 4);
    print_type("four_sinh_perp2_sq_ax_wax:", four_sinh_perp2_sq_ax_wax(w, p));
    print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w, p));
    print_type("diff:", diff);
    fprintf(stderr, "*******************************\n");
  }
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline const bool really_cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  if (g_debug && absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR) {
      fprintf(stderr, "Realy can't fix y_axis");
    print_type(fsp2sq);
    fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  return absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR; 
}

template<typename T>
inline const bool cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
   fprintf(stderr, "********** CANNOT FIX Y AXIS ***********\n");
   print_type(fsp2sq);
   fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
   fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  // return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  return absLB(fsp2sq) > 0; 
}

template<typename T>
inline const bool must_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJCC tests
  T diff = p.coshreD * 4 - four_cosh_dist_ay_way(w, p);
  if (g_debug && 
      //std::is_same<T, AJCC>::value && strictly_pos(diff)) {
      std::is_same<T, AJCC>::value) {
    fprintf(stderr, "********** MUST FIX Y AXIS ***********\n");
    print_SL2(w);
    print_type("4 cosh 2 dy:", p.coshreD * 4);
    print_type("four_sinh_perp2_sq_ay_way:", four_sinh_perp2_sq_ay_way(w, p));
    print_type("4 cosh dist ay way:", four_cosh_dist_ay_way(w, p));
    print_type("diff:", diff);
    fprintf(stderr, "*******************************\n");
  }
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

// Returns true if w and x are in the variety neighbourhood of each other:
// either Jørgensen(w,x) < 1 or Jørgensen(x,w) < 1 (forcing them to commute in a
// discrete group), or w provably maps axis(x) inside the tube of x (must_fix_x_axis).
// Note: must_fix_x_axis is only valid for AJCC arithmetic when x has trace near ±2.
template<typename T>
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  if (g_debug && (absUB(jorgensen_wx(w, params)) < 1 ||
      absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params))) {
      fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w, params)),
        absUB(jorgensen_xw(w, params)), must_fix_x_axis(w, params));
  }
  return absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params);
}

// Analogous to inside_var_nbd_x but for y and the y-tube.
template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  if (g_debug && std::is_same<T, AJCC>::value && (absUB(jorgensen_wy(w, params)) < 1 ||
    absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params))) {
    fprintf(stderr, "UB Jwy %f, UB Jyw %f, must_fix %d\n", absUB(jorgensen_wy(w, params)),
      absUB(jorgensen_yw(w, params)), must_fix_y_axis(w, params));
  }
  return absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params);
}


// Returns true if w maps axis(y) provably closer to axis(x) than cosh(Re D):
// i.e., dist(axis(x), w·axis(y)) < Re(D). In a discrete group with embedded tubes
// this cannot happen, so the box is eliminated.
template<typename T>
inline const bool moves_y_axis_too_close_to_x(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshreD * 4 - four_cosh_dist_ax_way(w, p);
  if (g_debug && std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVES Y TOO CLOSE TO X\n");
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshreD * 4);
    print_type("4coshd(dist(x-axis, w(y-axis))):", four_cosh_dist_ax_way(w, p)); 
    T z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
      ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff);
}

// Analogous to moves_y_axis_too_close_to_x: certifies dist(axis(y), w·axis(x)) < Re(D).
template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshreD * 4 - four_cosh_dist_ay_wax(w, p);
  if (g_debug && std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVES X TOO CLOSE TO Y\n");
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshreD * 4);
    print_type("4coshd(dist(y-axis, w(x-axis))):", four_cosh_dist_ay_wax(w, p)); 
    T z = ((w.a * w.a) * p.expmD2 - (w.b * w.b) * p.expD2 ) * p.expmD2 +
      ((w.d * w.d) * p.expD2  - (w.c * w.c) * p.expmD2) * p.expD2;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moved_y_axis_not_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_way(w, p);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool moved_x_axis_not_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_wax(w, p);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

// Returns true if the Margulis function for (w1,w2) is provably less than cosh(μ),
// i.e., 4·cosh(margulis(w1,w2)) < 4·cosh(μ). This certifies that the pair (w1,w2)
// witnesses a Margulis number smaller than μ at every point in the box.
template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const Params<T>& p) {
  T diff = p.coshmu * 4 - four_cosh_margulis_simple(w1, w2).first;
  return strictly_pos(diff);
}

// Returns true if w moves the basepoint j by less than the Margulis number μ.
// cosh(d(j, w·j)) < cosh(μ) certifies that w is in the Margulis region for j.
template<typename T>
inline bool move_less_than_marg(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshmu - cosh_move_j(w);
  if (g_debug && std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVE J\n");
    fprintf(stderr, "word");
    print_SL2(w);
    fprintf(stderr, "cosh_move_j(w)\n");
    print_type("cosh_move_j(w):", cosh_move_j(w));
    print_type("cosh(mu):", p.coshmu);
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff); 
}

// Given that w fixes the same oriented axis as x_or_y, in a discrete torsion-free group
// they must generate a cyclic subgroup, so [x_or_y, w] = id. Returns true if we can
// certify the commutator is not the identity, ruling out w being a power of x_or_y.
template<typename T>
inline bool non_cyclic_power(const SL2<T>& w, const SL2<T>& x_or_y) {
  SL2<T> commutator = x_or_y * w * inverse(w * x_or_y);
  if (g_debug && std::is_same<T, AJCC>::value && not_identity(commutator)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "NOT CYCLIC POWER\n");
    fprintf(stderr, "x or y\n");
    print_SL2(x_or_y);
    fprintf(stderr, "(x or y)^2\n");
    print_SL2(x_or_y * x_or_y);
    fprintf(stderr, "w\n");
    print_SL2(w);
    fprintf(stderr, "commutator\n");
    print_SL2(commutator);
    fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
    fprintf(stderr, "****************************************\n");
  }
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return not_identity(commutator); 
}

// Meyerhoff K-test: certifies that the tube about axis(x) or axis(y) is larger than
// the tube bound encoded in four_cosh_tube_diam_UB.
//
// ch_o = cosh(Re L), cs_o = cos(Im L) for the generator's complex length L.
// Uses the Chebyshev recurrence cosh(nθ) = 2cosh(θ)cosh((n-1)θ) - cosh((n-2)θ)
// (and the same for cos) to step through multiples kL until cosh(Re(kL)) ≥ √2.
// For the smallest such k, the Meyerhoff formula gives a lower bound on the tube
// diameter. Returns true if that lower bound exceeds four_cosh_tube_diam_UB,
// proving the box can be killed. See Meyerhoff, "A lower bound for the volume of
// hyperbolic 3-manifolds."
#define MAX_MEYER 8
template<typename T>
bool meyerhoff_k_test(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  T ch_prev = T(1);
  T cs_prev = T(1);
  T ch = ch_o;
  T cs = cs_o;
  T temp;
  T four_cosh_tube_diam_LB;
  T meyer_k = T(1024); // sentinel: larger than any real k value we will find
  int count = 0;
  while (absUB(ch * ch) < 2 && count < MAX_MEYER) {
    temp = ch - cs;
    if (strictly_pos(meyer_k - temp) && absUB((temp + 1) * (temp + 1)) < 2) {
      meyer_k = temp;
      four_cosh_tube_diam_LB = sqrt(-(meyer_k * 32) + 16) / meyer_k;
      if (strictly_pos(four_cosh_tube_diam_LB - four_cosh_tube_diam_UB)) {
        if (g_debug) {
          fprintf(stderr, "Meyer k %f with 4 cosh tube diam LB %f and UB %f\n",
              absLB(meyer_k), absUB(four_cosh_tube_diam_LB), absLB(four_cosh_tube_diam_UB));
        }
        return true;
      }
    }
    // Chebyshev recurrence: advance both ch and cs one step.
    temp = (ch_o * 2) * ch - ch_prev;
    ch_prev = ch;
    ch = temp;
    temp = (cs_o * 2) * cs - cs_prev;
    cs_prev = cs;
    cs = temp;
    count += 1;
  }
  return false; // inconclusive
}

// Finds the longest (in Re) primitive root of L that the Meyerhoff K-test cannot kill.
// Starting from cosh(Re L)/cos(Im L), iterates square-root halvings (L → L/2) up to
// MAX_ROOTS times. Returns the cosh(Re L) value at the deepest root level still
// surviving the K-test, or 0 if the first level already fails. This gives a lower
// bound on the real part of the complex length of any primitive element fixing the axis.
#define MAX_ROOTS 8
template<typename T>
T worst_primitive_cosh_re_len(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  T ch_prev = ch_o;
  T cs_prev = cs_o;
  for (int i = 0; i < MAX_ROOTS; ++i) {
    T ch = sqrt((ch_prev + 1) / 2);
    T cs = sqrt((cs_prev + 1) / 2); // note, - pi <= Im(L) <= pi, so sign is +
    if (meyerhoff_k_test(ch, cs, four_cosh_tube_diam_UB)) {
      return ch_prev;
    }
    ch_prev = ch;
    cs_prev = cs;
  }
  // no luck
  T zero(0);
  return zero; 
}

// Polynomial lower bound on cosh(μ) as a function of 2·sinh(Re D/2) (the tube radius).
// Coefficients are from the bilipschitz paper (Futer-Purcell-Schleimer). The polynomial
// was fitted to be a guaranteed lower bound over the relevant range; the unusual ordering
// of the summation reduces accumulated floating-point error.
template<typename T>
T cosh_marg_lower_bound(const T& two_sinh_r) {
  T s = two_sinh_r;
  T a8 = powT(s, 8) * (-0.0000014461700558); 
  T a7 = powT(s, 7) *   0.0000365880448817; 
  T a6 = powT(s, 6) * (-0.0003163830157272);
  T a5 = powT(s, 5) *   0.0005316504647188;
  T a4 = powT(s, 4) *   0.0086912125268823;
  T a3 = powT(s, 3) * (-0.061949675652791);
  T a2 = powT(s, 2) *   0.151649220047696;
  T a1 = s          * (-0.01513801009421);
  double a0 = 0.9999999;
  return ((a8 + (a1 + a0)) + (a4 + a5)) + ((a7 + a2) + (a6 + a3)); 
}

// Returns true only if EVERY word-specific kill condition in evaluate_AJCC's single-word
// path is provably negative over the entire box (re(diff) < 0, via strictly_pos(-diff)).
// Because each descendant box is contained in this box, a true result also holds for the
// whole subtree, so the word can be removed from the candidate set for all descendants.
//
// SOUNDNESS: a false positive can only cause a missed elimination (an extra hole) — never
// a false elimination — because the pruned word is merely skipped, not asserted to kill.
// So over-conservatism only loses pruning opportunities; it can never corrupt the proof.
//
// This mirrors exactly the four word-specific conditions in evaluate_AJCC (move,
// wy_hits_sym_axis, moves_y_axis_too_close_to_x, inside_var_nbd_y). It deliberately does
// NOT consider the box-level checks (qrs / vol3 / sym3): those do not depend on this word,
// and a receding word already has do_eval == false, so it never triggers them anyway.
template<typename T>
inline bool word_provably_receding(const SL2<T>& w, const std::string& word,
                                   const Params<T>& p) {
  // 1. move_less_than_marg fires on strictly_pos(coshmu - cosh_move_j).
  if (!strictly_pos(cosh_move_j(w) - p.coshmu)) return false;
  // 2. wy_hits_sym_axis: three tube tests, each fires on strictly_pos(two_cosh_re_tube - dist).
  //    Proving the first clause of each receding (dist > tube) is sufficient to bar the kill.
  T tsp2sq_inf = two_sinh_perp2_sq_way_zero_inf(w, p);
  if (!strictly_pos(two_cosh_dist(tsp2sq_inf) - p.twocoshreD2)) return false;
  T fsp2sq_one = four_sinh_perp2_sq_way_mp_one(w, p);
  if (!strictly_pos(four_cosh_dist(fsp2sq_one) - p.twocoshreD2 * 2)) return false;
  T fsp2sq_iye = four_sinh_perp2_sq_way_mp_iye(w, p);
  if (!strictly_pos(four_cosh_dist(fsp2sq_iye) - p.twocoshreD2 * 2)) return false;
  // Conditions 3 and 4 only fire when x_power(word) > 0 (a word-structural, box-invariant
  // property), exactly as guarded in evaluate_AJCC.
  if (x_power(word) > 0) {
    // 3. moves_y_axis_too_close_to_x fires on strictly_pos(coshreD*4 - four_cosh_dist_ax_way).
    if (!strictly_pos(four_cosh_dist_ax_way(w, p) - p.coshreD * 4)) return false;
    // 4. inside_var_nbd_y is a disjunction (jorgensen_wy<1 || jorgensen_yw<1 || must_fix_y);
    //    all three must be provably barred. jorgensen_* and four_cosh_dist_* are real-valued.
    if (!strictly_pos(jorgensen_wy(w, p) - 1.0)) return false;
    if (!strictly_pos(jorgensen_yw(w, p) - 1.0)) return false;
    if (!strictly_pos(four_cosh_dist_ay_way(w, p) - p.coshreD * 4)) return false;
  }
  return true;
}

// Attempts to find a conjugate of `word` (by prefixing up to MAX_ID_SHIFT copies of
// y or Y) that lies inside the variety neighbourhood of y and has real length provably
// less than the shortest primitive element fixing axis(y). Such a conjugate must be the
// identity in any torsion-free discrete group, giving a proven relator.
// Returns the conjugated word string on success, or "" if inconclusive.
// Only applies when x_power(word) > 0 (word involves x).
#define MAX_ID_SHIFT 5
template<typename T>
std::string proven_identity(const std::string& word, const Params<T>& p) {
  SL2<T> w = construct_word(word, p);
  std::string new_word;
  if (x_power(word) > 0 && inside_var_nbd_y(w, p)) {
    // construct_x and the primitive-length bound depend only on p; build them only when
    // a word actually reaches this guard (far-from-manifold words bail out above).
    SL2<T> x = construct_x(p);
    T four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, p);
    T cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshreL, p.cosimL, four_cosh_y_tube_UB);
    if (g_debug && std::is_same<T, AJCC>::value) {
      print_type("four_cosh_y_tube_UB", four_cosh_y_tube_UB);
      print_type("cosh_prim_re_len", cosh_prim_re_len);
    }
    for (auto s : {"y", "Y"}) {
      // new_word = y_strip(word);
      new_word = word;
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        if (g_debug && std::is_same<T, AJCC>::value) {
          fprintf(stderr, "Testing new word: %s\n", new_word.c_str());
        }
        SL2<T> new_w = construct_word(new_word, p); // order matters
        if (inside_var_nbd_y(new_w, p)) {
          T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
          if (strictly_pos(diff)) {
            if (g_debug && std::is_same<T, AJCC>::value) {
              fprintf(stderr, "Found proven identity: %s .\n", new_word.c_str());
            }
            return new_word;
          }
        }      
        new_word = s + new_word;
      }
    }
  }
  return "";
}

#endif //_TestCollection_
