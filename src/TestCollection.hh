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
  killed_only_elliptic = 3, // var_nbd of w^k but w not id
  killed_x_hits_y = 4, // w(axis(x)) closer to axis(y) than dx + dy
  killed_y_hits_x = 5, // w(axis(y)) closer to axis(x) than dx + dy
  killed_x_tube_dist = 6, // w(axis(x)) closer to axis(x) than 2dx but further than 0 (i.e. w not x^k)
  killed_y_tube_dist = 7, // w(axis(y)) closer to axis(y) than 2dy but further than 0 (i.e. w not y^k)
  killed_x_tube_not_power = 8, // w(axis(x)) closer to axis(x) than 2dx and provably w not x^k // FIXME need to be non-cyclic
  killed_y_tube_not_power = 9, // w(axis(y)) closer to axis(y) than 2dy and provably w not y^k // FIXME need to be non-cyclic
  killed_move = 10, // w(1j) moved less than marg
  killed_marg = 11, // w1 and w2 have (simple) margulis less than mu TODO should we do powers?
  variety_nbd_x = 12, // w and x fail Jorgensen 
  variety_nbd_y = 13, // w and y fail Jorgensen
  variety_nbd = 14, // w1 and w2 fail Jorgensen and one is not parabolic
  open_with_qr = 15,
  out_of_bounds_center = 16,
  variety_center = 17,
  x_hits_y_center = 18,
  y_hits_x_center = 19,
  bad_x_tube_center = 20,
  bad_y_tube_center = 21,
  bad_move_center = 22,
  bad_marg_center = 23,
  open = -1
} 
box_state;

struct ImpossibleRelations;

struct TestCollection {
	int size();
	box_state evaluateCenter(int index, Box& box);
	box_state evaluateBox(int index, Box& box, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<AJ> >& words_cache);
	const char* getName(int index);
	int add(const word_pair& pair);
	int add(std::string pair);
	void load(const char* fileName);
	void loadImpossibleRelations(const char* fileName);
private:
	std::map<word_pair, int> pairIndex;
	std::vector<word_pair> pairVector;
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

template<typename T>
inline const bool cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  return absLB(four_cosh_dist_ax_wax(w, p)) > 4;
}

template<typename T>
inline const bool must_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJ tests
  T diff = p.cosh2dx * 4 - four_cosh_dist_ax_wax(w, p);
  // print_SL2(w);
  // print_type("4 cosh 2 dx:", p.cosh2dx * 4);
  // print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w, p));
  // print_type("diff:", diff);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
inline const bool cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  return absLB(four_cosh_dist_ay_way(w, p)) > 4;
}

template<typename T>
inline const bool must_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJ tests
  T diff = p.cosh2dy * 4 - four_cosh_dist_ay_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when x has trace close to +/- 2
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
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_ay_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const Params<T>& p) {
  T diff = p.coshmu * 4 - four_cosh_margulis_simple(w1, w2).first;
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return absLB(diff) > 0 && re_center(diff) > 0;
}

template<typename T>
inline bool non_cylic_power_x(const SL2<T>& w, const Params<T>& p) {
  // Assume word fixes the same axis as x, so it must live in a cyclic group with x.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of x
  // TODO
  return false; 
}

template<typename T>
inline bool non_cylic_power_y(const SL2<T>& w, const Params<T>& p) {
  // Assume word fixes the same axis as y, so it must live in a cyclic group with y.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of y
  // TODO
  return false; 
}

