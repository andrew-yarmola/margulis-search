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
  variety_center = 9,
  bad_tubes_center = 10,
  bad_margulis_center = 11,
  open = -1
} 
box_state;

typedef std::pair<std::string, std::string> word_pair;

struct ImpossibleRelations;

struct TestCollection {
	int size();
	box_state evaluateCenter(int index, Box& box);
	box_state evaluateBox(int index, Box& box, std::string& aux_word, std::vector<std::string>& new_qrs,
                        std::unordered_map<int,ACJ>& para_cache,std::unordered_map<std::string,SL2<ACJ>>& words_cache);
	const char* getName(int index);
	int add(std::string word);
	void load(const char* fileName);
	void loadImpossibleRelations(const char* fileName);
private:
	std::map<std::string, int> stringIndex;
	std::vector<std::string> indexString;
	box_state evaluate_approx(word_pair pair, const Box& params);
  box_state evaluate_ACJ(word_pair pair, const Box& params, std::string& aux_word, std::vector<std::string>& new_qrs,
                         std::unordered_map<int,ACJ>& para_cache, std::unordered_map<std::string,SL2<ACJ>>& words_cache);
  bool ready_for_elliptics_test(SL2<ACJ>& w);
  bool only_elliptics(SL2<ACJ>& w, Params<ACJ>& params);
	ImpossibleRelations *impossible;
};

template<typename T>
inline const bool inside_var_nbd(const SL2<T>& w, const Params<T>& params) {
  return jorgensen_xw_UB(w, params) < 1 || jorgensen_wx_UB(w, params) < 1 ||
         jorgensen_yw_UB(w, params) < 1 || jorgensen_wy_UB(w, params) < 1;
        
}

template<typename T>
inline const bool not_identity(const SL2<T>& w) {
    return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
         ((absLB(w.a-1) > 0 || absLB(w.d-1) > 0) && (absLB(w.a+1) > 0 || absLB(w.d+1) > 0));
}

template<typename T>
inline bool tubes_intersect(SL2<T>& w, SL2<T>& x, SL2<T>& y, Params<T>& params) {
    double exp_2_t_x_LB = exp_2_t(x,y,false);
    double exp_2_t_y_LB = exp_2_t(y,x,false);
    double exp_2_re_perp_x_UB = e_re_perp_x_UB(w);
    double exp_2_re_perp_y_UB = e_re_perp_y_UB(w, params.coshP, params.sinhP);
    return (exp_2_re_perp_x_UB < exp_2_t_x_LB) || (exp_2_re_perp_y_UB < exp_2_t_y_LB);
}

template<typename T>
inline bool margulis_larger_than_cutoff(SL2<T>& w1, SL2<T>& w2, double cutoff) {
  // Cutoff is given as 4 cosh(margulis_cutoff)
  // TODO Optimize for x and y as w1 (or w2)
  double margulis = 10000;
  int n = 1;
  int m = 1;
  SL2<T> A(w1);
  while (four_cosh_re_length_LB(A) < cutoff) {
    SL2<T> B(w2);
    while (four_cosh_re_length_LB(B) < cutoff) {
      margulis = min(margulis, four_cosh_margulis(A,B,false)); 
      m += 1;
      B = pow(w2,m); // reduing number of powers needed, might be better to just accumuate
    }
    n += 1;
    A = pow(w2,n); // reducing the number of powers needed, might be better to just accumulate
  }
  return margulis >= cutoff; 
}
