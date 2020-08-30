#include <unordered_map>
#include <string>
#include <vector>
#include "types.hh"
#include "Box.h"
#include "SL2.hh"
#include "IsomH3.hh"

typedef enum _box_state
{
  killed_bounds = 1,
  // killed_powers = 2, // power of x ot y moves 1j less than mu (see killed_move)
  killed_only_elliptic = 2, // var_nbd of w^k but w not id
  killed_x_hits_y = 3, // w(axis(x)) closer to axis(y) than dx + dy
  killed_y_hits_x = 4, // w(axis(y)) closer to axis(x) than dx + dy
  killed_x_tube = 5, // w(axis(x)) closer to axis(x) than 2dx but further than 0 (i.e. w not x^k)
  killed_y_tube = 6, // w(axis(y)) closer to axis(y) than 2dy but further than 0 (i.e. w not y^k)
  killed_lox_not_x_power = 7, // w(axis(x)) closer to axis(x) than 2dx and provably w not x^k // FIXME need to be non-cyclic
  killed_lox_not_y_power = 8, // w(axis(y)) closer to axis(y) than 2dy and provably w not y^k // FIXME need to be non-cyclic
  killed_move = 9, // w(1j) moved less than marg, note power not produces by word search
  killed_marg = 10, // w1 and w2 have (simple) margulis less than mu TODO should we do powers?
  variety_nbd_x = 11, // w and x fail Jorgensen 
  variety_nbd_y = 12, // w and y fail Jorgensen
  variety_nbd = 13, // w1 and w2 fail Jorgensen and one is not parabolic
  killed_failed_qr = 27,
  open_with_qr = 14,
  out_of_bounds_center = 15,
  variety_center = 16,
  var_x_center = 17,
  var_y_center = 18,
  x_hits_y_center = 19,
  y_hits_x_center = 20,
  bad_x_tube_center = 21,
  bad_y_tube_center = 22,
  bad_lox_x_center = 23,
  bad_lox_y_center = 24,
  bad_move_center = 25,
  bad_marg_center = 26,
  open = -1
} 
box_state;

struct ImpossibleRelations;

struct TestCollection {
  int size();
  box_state evaluate_center(int index, Box& box);
  box_state evaluate_box(int index, Box& box, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<AJ> >& words_cache);
  const std::string get_name(int index);
  int add(word_pair pair);
  int add(std::string pair);
  void load(const char* fileName);
  void load_impossible_relations(const char* fileName);
  private:
  std::map<word_pair, int> pair_index;
  std::vector<word_pair> pair_vector;
  box_state evaluate_approx(word_pair pair, const Box& params);
  box_state evaluate_AJ(word_pair pair, const Box& params, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<AJ> >& words_cache);
  bool ready_for_elliptics_test(SL2<AJ>& w);
  bool only_elliptics(SL2<AJ>& w, Params<AJ>& params);
  ImpossibleRelations *impossible;
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
  // print_type(fsp2sq);
  // printf("LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
  // printf("LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  // return absLB(four_cosh_dist_ax_wax(w, p)) > 4;
}

template<typename T>
inline const bool must_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJ tests
  T diff = p.cosh2dx * 4 - four_cosh_dist_ax_wax(w, p);
  // print_SL2(w);
  //if (strictly_pos(diff)) {
  //  print_type("4 cosh 2 dx:", p.cosh2dx * 4);
  //  print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w, p));
  //  print_type("diff:", diff);
  //}
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
  // if (absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params)) {
  //   fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w, params)),
  //          absUB(jorgensen_xw(w, params)), must_fix_x_axis(w, params));
  // }
  return absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params);
}

template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when y has trace close to +/- 2
  // print_type("jorg_wy", jorgensen_wy(w, params));
  // print_type("jorg_yw", jorgensen_yw(w, params));
  return absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params);
}


template<typename T>
inline const bool moves_y_axis_too_close_to_x(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_ax_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  //  if (strictly_pos(diff)) {
  //    printf("SL2 of word:\n");
  //    print_SL2(w);
  //    print_type("4cosh(dx+dy):", p.coshdxdy * 4);
  //    print_type("4coshd(dist(x-axis, w(x-axis))):", four_cosh_dist_ax_way(w, p)); 
  //    T z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
  //          ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
  //    print_type("4 sinh^2(dist/2) + 2:", z);
  //    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
  //    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
  //    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
  //  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_ay_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
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
  if (not_identity(commutator)) {
    return true;
  }
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return false; 
}

// Meyerhoff K Test
// We stop computing if we fail the test
#define MAX_MEYER 8
template<typename T>
bool meyerhoff_k_test(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB, bool debug) {
  // Assumed ch and cs are real valued jets
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
        if (debug) {
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


