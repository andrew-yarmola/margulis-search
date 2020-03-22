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
  killed_powers = 2, // power of x ot y moves 1j less
  killed_only_elliptic = 3, // var_nbd of w^k but w not id
  killed_x_hits_y = 4, // w(axis(x)) closer to axis(y) than dx + dy
  killed_y_hits_x = 5, // w(axis(y)) closer to axis(x) than dx + dy
  killed_x_tube = 6, // w(axis(x)) closer to axis(x) than 2dx but w not x^k
  killed_y_tube = 7, // w(axis(y)) closer to axis(y) than 2dy but w not y^k
  killed_move = 8, // w(1j) moved less than marg
  killed_marg = 9, // w1 and w2 have (simple) margulis less than mu TODO should we do powers?
  variety_nbd_x = 10, // w and x fail Jorgensen 
  variety_nbd_y = 11, // w and y fail Jorgensen
  variety_nbd = 12, // w1 and w2 fail Jorgensen and one is not parabolic
  open_with_qr = 13,
  out_of_bounds_center = 14,
  variety_center = 15,
  x_hits_y_center = 16,
  y_hits_x_center = 17,
  bad_x_tube_center = 18,
  bad_y_tube_center = 19,
  bad_move_center = 20,
  bad_margulis_center = 21,
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
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when x has trace close to +/- 2
  return jorgensen_wx_UB(w, params) < 1 || jorgensen_xw_UB(w, params) < 1;
}

template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when y has trace close to +/- 2
  return jorgensen_wy_UB(w, params) < 1 || jorgensen_yw_UB(w, params) < 1;
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
inline bool tubes_intersect_x(const SL2<T>& w, const SL2<T>& x, const SL2<T>& y, const Params<T>& params) {
    double exp_2_t_x_LB = exp_2_t(x,y,false);
    fprintf(stderr, "Tube around x is %f\n", exp_2_t_x_LB); 
    if (exp_2_t_x_LB < 0) {
      fprintf(stderr, "Failed tube computation with error %f\n", exp_2_t_x_LB); 
      return false;
    } else {
      double exp_re_perp_x_UB = e_re_perp_x_UB(w);
      fprintf(stderr, "Dist to x conj %f\n", exp_re_perp_x_UB); 
      return exp_re_perp_x_UB < exp_2_t_x_LB;
    }
}

template<typename T>
inline bool tubes_intersect_y(const SL2<T>& w, const SL2<T>& x, const SL2<T>& y, const Params<T>& params) {
    double exp_2_t_y_LB = exp_2_t(y,x,false);
    fprintf(stderr, "Tube around y is %f\n", exp_2_t_y_LB); 
    if (exp_2_t_y_LB < 0) {
      fprintf(stderr, "Failed tube computation with error %f\n", exp_2_t_y_LB); 
      return false;
    } else {
      double exp_re_perp_y_UB = e_re_perp_y_UB(w, params.coshP, params.sinhP);
      fprintf(stderr, "Dist to y conj  %f\n", exp_re_perp_y_UB); 
      return exp_re_perp_y_UB < exp_2_t_y_LB;
    }
}

template<typename T>
inline bool margulis_larger_than_cutoff(const SL2<T>& w1, const SL2<T>& w2, double cutoff) {
  // Cutoff is given as 4 cosh(margulis_cutoff)
  // TODO Optimize for x and y as w1 (or w2)
  double margulis_LB = four_cosh_margulis(w1, w2, false);
  if (margulis_LB >= cutoff) return true;
  else return false;
}

template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const SL2<T>& x, const SL2<T>& y) {
  // TODO Compute x y margulis once for each box
  double margulis_UB = four_cosh_margulis(w1, w2, true);
  double margulis_xy_LB = four_cosh_margulis(x, y, false);
  if (margulis_UB < margulis_xy_LB) return true;
  else return false;
}
