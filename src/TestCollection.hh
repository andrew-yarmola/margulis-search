#ifndef _TestCollection_ 
#define _TestCollection_
#include <unordered_map>
#include <string>
#include <vector>
#include "types.hh"
#include "Box.h"
#include "SL2.hh"
#include "IsomH3.hh"
#include "ImpossibleRelations.h"

extern bool g_debug;

struct ImpossibleRelations;

struct TestCollection {
  int size();
  box_state evaluate_center(int index, Box& box);
  box_state evaluate_box(int index, Box& box, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<AJ> >& words_cache);
  const std::string get_name(int index);
  word_pair get_pair(int index);
  int add(word_pair pair);
  int add(std::string pair);
  void load(const char* file_path);
  void load_impossible_relations(const char* file_path);
  ImpossibleRelations *impossible;
  private:
  word_pair parse_word_pair(std::string buf);
  std::map<word_pair, int> pair_index;
  std::vector<word_pair> pair_vector;
  box_state evaluate_approx(word_pair pair, const Box& params);
  box_state evaluate_AJ(word_pair pair, const Box& params, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<AJ> >& words_cache);
  bool ready_for_elliptics_test(SL2<AJ>& w);
  bool only_elliptics(SL2<AJ>& w, Params<AJ>& params);
};

template<typename T>
inline const bool inside_var_nbd(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting non-parabolics when discrete
  // So either we have a relator or they are elliptic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_parabolic(w1) || not_parabolic(w2));
}

template<typename T>
inline const bool inside_var_nbd_ne(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting loxodromics when discrete, so we have a relator
  // Note,we must test both as elliptic can commute with loxodromic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_elliptic_or_parabolic(w1) && not_elliptic_or_parabolic(w2));
}

template<typename T>
inline const bool not_parabolic(const SL2<T>& w) {
  return absLB(w.a + w.d - 2) > 0 && absLB(w.a + w.d + 2) > 0;
}

template<typename T>
inline const bool not_elliptic_or_parabolic(const SL2<T>& w) {
  T tr = w.a + w.d;
  return absLB(tr) > 2 || absLB(tr - conj(tr)) > 0;
}

template<typename T>
inline const bool not_identity(const SL2<T>& w) {
  return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
    ((absLB(w.a-1) > 0 || absLB(w.d-1) > 0) && (absLB(w.a+1) > 0 || absLB(w.d+1) > 0));
}

#define LERR 0.000000000001

template<typename T>
inline const bool really_cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  // print_type(fsp2sq);
  // printf("LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
  // printf("LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  return absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR; 
}

template<typename T>
inline const bool cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
//  if (std::is_same<T, AJ>::value && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
//    fprintf(stderr, "********** CANNOT FIX X AXIS ***********\n");
//    T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
//    print_type(fsp2sq);
//    fprintf(stderr, "Can't fix x axis LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
//    fprintf(stderr,"Can't fix x axis LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
//    fprintf(stderr, "*******************************\n");
//  }
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool must_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJ tests
  T diff = p.cosh2dx * 4 - four_cosh_dist_ax_wax(w, p);
//  if (std::is_same<T, AJ>::value && strictly_pos(diff)) {
//    fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
//    print_SL2(w);
//    print_type("4 cosh 2 dx:", p.cosh2dx * 4);
//    print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w, p));
//    print_type("diff:", diff);
//    fprintf(stderr, "*******************************\n");
//  }
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline const bool really_cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  // print_type(fsp2sq);
  // printf("LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
  // printf("LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  return absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR; 
}

template<typename T>
inline const bool cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  // print_type(fsp2sq);
  // printf("LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
  // printf("LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  // return absLB(four_cosh_dist_ay_way(w, p)) > 4;
}

template<typename T>
inline const bool must_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJ tests
  T diff = p.cosh2dy * 4 - four_cosh_dist_ay_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when x has trace close to +/- 2
//  if (absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params)) {
//    fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w, params)),
//            absUB(jorgensen_xw(w, params)), must_fix_x_axis(w, params));
//  }
  return absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params);
}

template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when y has trace close to +/- 2
//  if (absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params)) {
//    fprintf(stderr, "UB Jwy %f, UB Jyw %f, must_fix %d\n", absUB(jorgensen_wy(w, params)),
//            absUB(jorgensen_yw(w, params)), must_fix_y_axis(w, params));
//  }
  return absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params);
}


template<typename T>
inline const bool moves_y_axis_too_close_to_x(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_ax_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
//  if (std::is_same<T, AJ>::value && strictly_pos(diff)) {
//    fprintf(stderr, "****************************************\n");
//    fprintf(stderr, "MOVES Y TOO CLOSE TO X\n");
//    print_SL2(w);
//    print_type("4cosh(dx+dy):", p.coshdxdy * 4);
//    print_type("4coshd(dist(x-axis, w(y-axis))):", four_cosh_dist_ax_way(w, p)); 
//    T z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
//           ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
//    print_type("4 sinh^2(dist/2) + 2:", z);
//    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
//    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
//    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
//    print_type("diff:", diff);
//    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
//    fprintf(stderr, "****************************************\n");
//  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_ay_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
//  if (std::is_same<T, AJ>::value && strictly_pos(diff)) {
//    fprintf(stderr, "****************************************\n");
//    fprintf(stderr, "MOVES X TOO CLOSE TO Y\n");
//    print_SL2(w);
//    print_type("4cosh(dx+dy):", p.coshdxdy * 4);
//    print_type("4coshd(dist(y-axis, w(x-axis))):", four_cosh_dist_ay_wax(w, p)); 
//    T z = ((w.a * w.a) * p.expmdx - (w.b * w.b) * p.expdx ) * p.expmdyf +
//      ((w.d * w.d) * p.expdx  - (w.c * w.c) * p.expmdx) * p.expdyf;
//    print_type("4 sinh^2(dist/2) + 2:", z);
//    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
//    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
//    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
//    print_type("diff:", diff);
//    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
//    fprintf(stderr, "****************************************\n");
//  }
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

template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const Params<T>& p) {
  T diff = p.coshmu * 4 - four_cosh_margulis_simple(w1, w2).first;
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline bool move_less_than_marg(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshmu - cosh_move_j(w);
  return strictly_pos(diff); 
}

template<typename T>
inline bool non_cylic_power(const SL2<T>& w, const SL2<T>& x_or_y) {
  // Assume word fixes the same axis as x or y, so it must live in a cyclic group with x or y.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of x or y
  SL2<T> commutator = x_or_y * w * inverse(w * x_or_y);
//  if (std::is_same<T, AJ>::value && not_identity(commutator)) {
//    fprintf(stderr, "****************************************\n");
//    fprintf(stderr, "NOT CYCLIC POWER\n");
//    fprintf(stderr, "x or y\n");
//    print_SL2(x_or_y);
//    fprintf(stderr, "(x or y)^2\n");
//    print_SL2(x_or_y * x_or_y);
//    fprintf(stderr, "w\n");
//    print_SL2(w);
//    fprintf(stderr, "commutator\n");
//    print_SL2(commutator);
//    fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
//    fprintf(stderr, "****************************************\n");
//  }
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return not_identity(commutator); 
}

// Meyerhoff K Test
// We stop computing if we fail the test
#define MAX_MEYER 8
template<typename T>
bool meyerhoff_k_test(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  T ch_prev = T(1);
  T cs_prev = T(1);
  T ch = ch_o;
  T cs = cs_o;
  T temp, four_cosh_tube_diam_LB;
  T meyer_k = T(1024); // arbitray large enough number
  int count = 0;
  while (absUB(ch * ch) < 2 && count < MAX_MEYER) {
    temp = ch - cs; 
    if (strictly_pos(meyer_k - temp) && absUB((temp + 1) * (temp + 1)) < 2) {
      meyer_k = temp;
      // See Meyerhoff paper on volume lowerbounds for hyperbolic 3-manifolds
      four_cosh_tube_diam_LB = sqrt(-(meyer_k * 32) + 16) / meyer_k;
      if (strictly_pos(four_cosh_tube_diam_LB - four_cosh_tube_diam_UB)) {
        if (g_debug) {
          fprintf(stderr, "Meyer k %f with 4 cosh tube diam LB %f and UB %f\n",
              absLB(meyer_k), absUB(four_cosh_tube_diam_LB), absLB(four_cosh_tube_diam_UB));
        }
        return true; // box can be killed
      }
    } 
    // Use Chebyshev recurrance relation
    temp = (ch_o * 2) * ch - ch_prev;
    ch_prev = ch;
    ch = temp;  
    T temp = (cs_o * 2) * cs - cs_prev;
    cs_prev = cs;
    cs = temp;
    count +=1;
  }
  return false; // inconclusive
}

#define MAX_ROOTS 8
template<typename T>
T worst_primitive_cosh_re_len(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
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

template<typename T>
T cosh_marg_lower_bound(const T& sinh_r) {
  T s = sinh_r;
  T a8 = powT(s, 8) * (-0.002012744207511); 
  T a7 = powT(s, 7) *   0.050422869707363; 
  T a6 = powT(s, 6) * (-0.2800449482233);
  T a5 = powT(s, 5) *   0.6738467122499;
  T a4 = powT(s, 4) * (-0.730897277659114);
  T a3 = powT(s, 3) *   0.1178833280583;
  T a2 = powT(s, 2) *   0.390674936173773;
  T a1 = s          *   0.001212870129678;
  double a0 = 0.999972595620724;
  return ((a8 + (a1 + a0)) + (a4 + a5)) + ((a7 + a2) + (a6 + a3)); 
}

#define MAX_ID_SHIFT 5
template<typename T>
std::string proven_identity(std::string word, const Params<T>& p) {
  SL2<T> w = construct_word(word, p);
  SL2<T> x = construct_x(p);
  SL2<T> y = construct_y(p);
  std::string new_word;
  if (inside_var_nbd_x(w, p)) {
    T four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, p);
    T cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshlx, p.costx, four_cosh_x_tube_UB); 
    for (auto s : {"x", "X"}) {
      new_word = x_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        new_word = s + new_word;
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          return new_word;
        }      
      }
    }
  }
  if (inside_var_nbd_y(w, p)) {
    T four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, p);
    T cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshly, p.costy, four_cosh_y_tube_UB); 
    for (auto s : {"y", "Y"}) {
      new_word = y_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        new_word = s + new_word;
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          return new_word;
        }      
      }
    }
  }
  return "";
}

#endif //_TestCollection_
