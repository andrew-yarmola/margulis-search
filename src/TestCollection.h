/*
 *  TestCollection.h
 *  mom
 *
 *  Created by Nathaniel Thurston on 27/09/2007.
 *  Copyright 2007 THingith ehf.. All rights reserved.
 *
 */

#include <unordered_map>
#include <string>
#include <vector>
#include "Params.h"
#include "Box.h"
#include "SL2.hh"

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
  open = -1
} 
box_state;

typedef enum _center_box_state
{
 variety_center = 1,
 out_of_bounds_center = 2,
 open = -1
} 
center_box_state;

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
	box_state evaluate_approx(std::string word, Params<Complex>& params);
  box_state evaluate_ACJ(std::string word, Params<ACJ>& params, std::string& aux_word, std::vector<std::string>& new_qrs,
                         std::unordered_map<int,ACJ>& para_cache, std::unordered_map<std::string,SL2<ACJ>>& words_cache);
  bool ready_for_elliptics_test(SL2<ACJ>& w);
  bool only_elliptics(SL2<ACJ>& w, Params<ACJ>& params);
	ImpossibleRelations *impossible;
};

inline const bool maybe_variety(const SL2C& w) {
    return (absUB(w.c) < 1) && (absUB(w.b) < 1);
}

inline const bool inside_var_nbd(const SL2<ACJ>& w)
{
    return (absUB(w.c) < 1) && (absUB(w.b) < 1);
}

inline const bool not_para_fix_inf(const SL2<ACJ>&x) {
    return absLB(x.c) > 0 
        || ((absLB(x.a-1) > 0 || absLB(x.d-1) > 0) && (absLB(x.a+1) > 0 || absLB(x.d+1) > 0));
}

inline const bool not_identity(const SL2<ACJ>&x) {
    return absLB(x.b) > 0 || not_para_fix_inf(x);
}

inline const bool maybe_large_horoball(const SL2C& w, const Params<Complex>& params) {
    return absUB((w.c / params.loxodromic_sqrt).z) < 1;
}

inline bool large_horoball(SL2<ACJ>& w, const Params<ACJ>& params) {
    return absUB(w.c / params.loxodromic_sqrt) < 1;
}
