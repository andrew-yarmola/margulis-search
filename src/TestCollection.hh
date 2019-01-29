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
  killed_bad_tubes = 2,
  killed_margulis = 3,
  killed_fake_elliptic = 4,
  killed_only_elliptic = 5,
  killed_failed_qr = 6,
  variety_nbd = 7,
  open_with_qr = 8,
  out_of_bounds_center = 9,
  variety_center = 9,
  bad_tubes_center = 10,
  bad_margulis_center = 11,
  open = -1
} 
box_state;

struct ImpossibleRelations;

struct TestCollection {
	int size();
	box_state evaluateCenter(int index, Box& box);
	box_state evaluateBox(int index, Box& box, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<ACJ> >& words_cache);
	const char* getName(int index);
	int add(const word_pair& pair);
	int add(std::string pair);
	void load(const char* fileName);
	void loadImpossibleRelations(const char* fileName);
private:
	std::map<word_pair, int> pairIndex;
	std::vector<word_pair> pairVector;
	box_state evaluate_approx(word_pair pair, const Box& params);
  box_state evaluate_ACJ(word_pair pair, const Box& params, std::string& aux_word, std::vector<std::string>& new_qrs, std::unordered_map<std::string,SL2<ACJ> >& words_cache);
  bool ready_for_elliptics_test(SL2<ACJ>& w);
  bool only_elliptics(SL2<ACJ>& w, Params<ACJ>& params);
	ImpossibleRelations *impossible;
};

template<typename T>
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  return jorgensen_xw_UB(w, params) < 1 || jorgensen_wx_UB(w, params) < 1;
        
}

template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  return jorgensen_yw_UB(w, params) < 1 || jorgensen_wy_UB(w, params) < 1;
        
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
