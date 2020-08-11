#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TestCollection.hh"
#include "ImpossibleRelations.h"
#include <algorithm>
using namespace std;
// using namespace __gnu_cxx;

extern double g_cosh_marg_upper_bound;
extern double g_cosh_marg_lower_bound;
extern double g_sinh_d_bound; 

int num_bound_tests = 5;

int TestCollection::size()
{
  return num_bound_tests + pair_vector.size();
}

box_state TestCollection::evaluate_approx(word_pair pair, const Box& box)
{
  //  fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n     Word Pair: %s and %s\n +++++++++++++++++++++++++++++++++++++++++++\n", pair.first.c_str(), pair.second.c_str());
  Params<Complex> p = box.center();
  if (pair.second.length() == 0) {
    SL2<Complex> w = construct_word(pair.first, p);
    if (not_identity(w)) {
      if (move_less_than_marg(w, p)) {
        return bad_move_center;
      }
      if (moves_y_axis_too_close_to_x(w,p)) {
        return y_hits_x_center;
      }
      if (moves_x_axis_too_close_to_y(w,p)) {
        return x_hits_y_center;
      }
    }
    if (inside_var_nbd_x(w, p)) {
      if (cant_fix_x_axis(w,p)) {
        return bad_x_tube_center;
      } else if (non_cylic_power(w, box.x_center())) {
        return bad_lox_x_center;
      } else {
        return var_x_center;
      }
    }
    if (inside_var_nbd_y(w, p)) {
      if (cant_fix_y_axis(w,p)) {
        return bad_y_tube_center;
      } else if (non_cylic_power(w, box.y_center())) {
        return bad_lox_y_center;
      } else {
        return var_y_center;
      }
    }
    //    if (must_fix_x_axis(w,p)) {
    //      if (cant_fix_x_axis(w,p)) {
    //        return bad_x_tube_center;
    //      } else if (non_cylic_power(w, box.x_center())) {
    //        return bad_lox_x_center;
    //      }
    //    }
    //    if (must_fix_y_axis(w,p)) {
    //      if (cant_fix_y_axis(w,p)) {
    //        return bad_y_tube_center;
    //      } else if (non_cylic_power(w, box.y_center())) {
    //        return bad_lox_y_center;
    //      }
    //    }
  } else {
    SL2<Complex> w1 = construct_word(pair.first, p);
    SL2<Complex> w2 = construct_word(pair.second,p);
    if (margulis_smaller_than_xy(w1, w2, p)) {
      return bad_marg_center;
    }
    if (inside_var_nbd(w1, w2)) {
      return variety_center;
    }
  }
  return open;
}

box_state TestCollection::evaluate_AJ(word_pair pair, const Box& box, string& aux_word,
    vector<string>& new_qrs, unordered_map< string,SL2<AJ> >& words_cache)
{
  //    fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n     Word Pair: %s and %s\n +++++++++++++++++++++++++++++++++++++++++++\n", pair.first.c_str(), pair.second.c_str());
  Params<AJ> p = box.cover();
  if (pair.second.length() == 0) {
    SL2<AJ> w = construct_word(pair.first, p);
    if (not_identity(w)) {
      if (move_less_than_marg(w, p)) {
        return killed_move;
      }
      if (moves_y_axis_too_close_to_x(w,p)) {
        return killed_y_hits_x;
      }
      if (moves_x_axis_too_close_to_y(w,p)) {
        return killed_x_hits_y;
      }
    }
    if (inside_var_nbd_x(w, p)) {
      if (cant_fix_x_axis(w,p)) {
        return killed_x_tube;
      } else if (non_cylic_power(w, box.x_cover())) {
        return killed_lox_not_x_power;
      } else {
        //        printf("x:\n");
        //        print_SL2(box.x_cover());
        //        printf("y:\n");
        //        print_SL2(box.y_cover());
        //        printf("w:\n");
        //        print_SL2(w);
        return variety_nbd_x;
      }
    }
    if (inside_var_nbd_y(w, p)) {
      if (cant_fix_y_axis(w,p)) {
        return killed_y_tube;
      } else if (non_cylic_power(w, box.y_cover())) {
        return killed_lox_not_y_power;
      } else {
        //        printf("y:\n");
        //        print_SL2(box.y_cover());
        //        printf("x:\n");
        //        print_SL2(box.x_cover());
        //        printf("w:\n");
        //        print_SL2(w);
        return variety_nbd_y;
      }
    }
    //    if (must_fix_x_axis(w,p)) {
    //      if (cant_fix_x_axis(w,p)) {
    //        return killed_x_tube;
    //      } else if (non_cylic_power(w, box.x_cover())) {
    //        return killed_lox_not_x_power;
    //      } else {
    //        new_qrs.push_back(pair.first);
    //        return open_with_qr;
    //      }
    //    }
    //    if (must_fix_y_axis(w,p)) {
    //      if (cant_fix_y_axis(w,p)) {
    //        return killed_y_tube;
    //      } else if (non_cylic_power(w, box.y_cover())) {
    //        return killed_lox_not_y_power;
    //      } else {
    //        new_qrs.push_back(pair.first);
    //        return open_with_qr;
    //      }
    //    }
  } else {
    SL2<AJ> w1 = construct_word(pair.first, p);
    SL2<AJ> w2 = construct_word(pair.second,p);
    if (margulis_smaller_than_xy(w1, w2, p)) {
      return killed_marg;
    }
    if (inside_var_nbd(w1, w2)) {
      return variety_nbd;
    }
  }
  return open;
}

box_state check_bounds_center(bool result) {
  if (result) return out_of_bounds_center;
  else return open;
}

box_state check_bounds(bool result) {
  if (result) return killed_bounds;
  else return open;
}


box_state TestCollection::evaluate_center(int index, Box& box)
{
  //  fprintf(stderr, "Evaluating center test index %d\n", index);
  Params<Complex> center = box.center();
  switch(index) {
    case 0:	{ // 1.0052 < cosh(0.104) <= cosh(mu) <= 0.
              return check_bounds_center(absUB(center.coshmu) < g_cosh_marg_lower_bound ||
                  absLB(center.coshmu) > g_cosh_marg_upper_bound ||
                  strictly_pos(-center.coshmu));
            } 
    case 1: { //
              return check_bounds_center(absLB(center.sinhdx) > g_sinh_d_bound || strictly_pos(-center.sinhdx) ||
                  absLB(center.sinhdy) > g_sinh_d_bound || strictly_pos(-center.sinhdy));
            }
    case 2: { // sin/cos bounds between -1 and 1
              return check_bounds_center(absLB(center.cosf) > 1 || absLB(center.sintx2) > 1 || absLB(center.sinty2) > 1);
            }
    case 3: { // check lengths are not negative 
              return check_bounds_center(absUB(center.coshlx) < 1 || absLB(center.coshly) < 1);
            }
    case 4: { // Meyerhoff tube bound. Check if embeded tube radius is more than rad + marg/2 
              SL2<Complex> x = construct_x(center);
              SL2<Complex> y = construct_y(center);
              Complex four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, center);
              Complex four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, center);
              return check_bounds_center(meyerhoff_k_test(center.coshlx, center.costx, four_cosh_x_tube_UB, false) || 
                  meyerhoff_k_test(center.coshly, center.costy, four_cosh_y_tube_UB, false));
            }
    default:
            return evaluate_approx(pair_vector[index - num_bound_tests], box);
  }
}

box_state TestCollection::evaluate_box(int index, Box& box, string& aux_word, vector<string>& new_qrs, unordered_map< string,SL2<AJ> >& words_cache)
{
  //  fprintf(stderr, "Evaluating box test index %d\n", index);
  Params<AJ> cover = box.cover();
  switch(index) {
    case 0:	{ // 1.0052 < cosh(0.104) <= cosh(mu) <= 0.
              return check_bounds(absUB(cover.coshmu) < g_cosh_marg_lower_bound ||
                  absLB(cover.coshmu) > g_cosh_marg_upper_bound ||
                  strictly_pos(-cover.coshmu));
            } 
    case 1: { //
              return check_bounds(absLB(cover.sinhdx) > g_sinh_d_bound || strictly_pos(-cover.sinhdx) ||
                  absLB(cover.sinhdy) > g_sinh_d_bound || strictly_pos(-cover.sinhdy));
            }
    case 2: { // sin/cos bounds between -1 and 1
              return check_bounds(absLB(cover.cosf) > 1 || absLB(cover.sintx2) > 1 || absLB(cover.sinty2) > 1);
            }
    case 3: { // check lengths are not negative 
              return check_bounds(absUB(cover.coshlx) < 1 || absLB(cover.coshly) < 1);
            }
    case 4: { // Meyerhoff tube bound. Check if embeded tube radius is more than rad + marg/2 
              fprintf(stderr, "%s", box.desc().c_str());
              SL2<AJ> x = construct_x(cover);
              SL2<AJ> y = construct_y(cover);
              AJ four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, cover);
              AJ four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, cover);
              return check_bounds(meyerhoff_k_test(cover.coshlx, cover.costx, four_cosh_x_tube_UB, true) || 
                  meyerhoff_k_test(cover.coshly, cover.costy, four_cosh_y_tube_UB, true));
            }
    default:
            return evaluate_AJ(pair_vector[index - num_bound_tests], box, aux_word, new_qrs, words_cache);
  }
}

// Returns the index number for the first basic 2 tests
// or the quasi-relator if the index is 2 or above
const string TestCollection::get_name(int index)
{
  // static char buf[500];
  if (index < num_bound_tests) {
    return to_string(index);
    // sprintf(buf, "%d", index);
  } else {
    word_pair p = pair_vector[index - num_bound_tests];
    return "(" + p.first + "," + p.second + ")";
    //		sprintf(buf, "(%s,%s)", p.first.c_str(), p.second.c_str());
  }
}

int TestCollection::add(string buf)
{
  size_t start = buf.find('(');   
  size_t comma = buf.find(',');   
  string first;
  string second;
  // fprintf(stderr, "Adding test: %s\n", buf.c_str());
  if (start != string::npos) {
    size_t end = comma;
    first = buf.substr(start + 1, end - start - 1);
  } else {
    return -1;
  } 
  if (comma != string::npos) {
    size_t end = buf.find(')');
    second = buf.substr(comma + 1, end - comma - 1);
  } else {
    return -1;
  }
  word_pair p(first, second);
  return add(p);
}

int TestCollection::add(word_pair p) { 
  map< word_pair,int >::iterator it = pair_index.find(p);
  if (it == pair_index.end()) {
    //    fprintf(stderr, "Adding test: (%s,%s)\n", p.first.c_str(), p.second.c_str());
    pair_index[p] = pair_vector.size();
    pair_vector.push_back(p);
    return pair_vector.size() + num_bound_tests - 1;
  } else {
    return it->second + num_bound_tests;
  }
}

void TestCollection::load(const char* fileName)
{
  FILE *fp = fopen(fileName, "r");
  char buf[1024];
  while (fp && fgets(buf, sizeof(buf), fp)) {
    int n = strlen(buf);
    if (!isalpha(buf[n-1]))
      --n;
    add(string(buf, n));
  }
}

void TestCollection::load_impossible_relations(const char* file_name)
{
  impossible = ImpossibleRelations::create(file_name);
}

