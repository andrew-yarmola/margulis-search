#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TestCollection.hh"
#include <algorithm>
using namespace std;

set<pair<string, string> > vol3_rels = {
    {"XYXYxYXYXy", "yXYxYxyxYxYY"},
    {"XyxyXyXYXy", "XyxyxYxyxy"},
    {"XyxyXyXYXy", "XyXYXYxYXY"},
    {"XyxyXyXYXy", "XYXYxYXYXy"},
    {"XyxyxYxyxy", "XyXyxyXyXY"},
    {"XyxyxYxyxy", "XyXYXyXyxy"},
    {"XYXYxYXYXy", "yXYxYxyxYxYY"},
    {"XYXYxYXYXy", "yxYXYxYxyxYY"},
    {"XYXYxYXYXy", "yxYxyxYxYXYY"},
    {"XYXYxYXYXy", "yxyxYxYXYxYY"},
    {"XYXYxYXYXy", "yyXyXYXyXyxY"},
    {"XYXyXYXYxY", "yXYxYxyxYxYY"},
    {"XYXyXYXYxY", "yxYXYxYxyxYY"},
    {"XYXyXYXYxY", "yxYxyxYxYXYY"},
    {"XYXyXYXYxY", "yxyxYxYXYxYY"},
    {"XYXyXYXYxY", "yyXyXYXyXyxY"},
    {"XYXyXyxyXy", "XYXYxYXYXy"},
    {"XYXyXyxyXy", "XYXyXYXYxY"},
    {"XYXyXyxyXy", "XyxyxYxyxy"},
    {"XYXyXyxyXy", "yXYXyXYXYxYY"},
    {"XYXyXyxyXy", "yXYxYXYXyXYY"},
    {"XYXyXyxyXy", "yXyXYXYxYXYY"},
    {"XYXyXyxyXy", "yxYXYXyXYXYY"},
    {"XYXyXyxyXy", "yyXyxyxYxyxY"},
    {"XYXyXyxyXy", "yyxYxyxyXyxY"},
    {"XYXyXyxyXy", "yyxyXyxyxYxY"},
    {"XYxYXYXyXY", "XYXyXyxyXy"},
    {"XYxYXYXyXY", "yXYxYxyxYxYY"},
    {"XYxYXYXyXY", "yxYXYxYxyxYY"},
    {"XYxYXYXyXY", "yxYxyxYxYXYY"},
    {"XYxYXYXyXY", "yxyxYxYXYxYY"},
    {"XYxYXYXyXY", "yyXyXYXyXyxY"},
    {"XYxYxyxYxY", "XYXYxYXYXy"},
    {"XYxYxyxYxY", "XYXyXYXYxY"},
    {"XYxYxyxYxY", "XYxYXYXyXY"},
    {"XYxYxyxYxY", "yXYXyXYXYxYY"},
    {"XYxYxyxYxY", "yXYxYXYXyXYY"},
    {"XYxYxyxYxY", "yXyXYXYxYXYY"},
    {"XYxYxyxYxY", "yxYXYXyXYXYY"},
    {"XYxYxyxYxY", "yyXyxyxYxyxY"},
    {"XYxYxyxYxY", "yyxYxyxyXyxY"},
    {"XYxYxyxYxY", "yyxyXyxyxYxY"},
    {"XyXYXYxYXY", "XYXyXyxyXy"},
    {"XyXYXYxYXY", "XYxYxyxYxY"},
    {"XyXYXYxYXY", "yXYxYxyxYxYY"},
    {"XyXYXYxYXY", "yxYXYxYxyxYY"},
    {"XyXYXYxYXY", "yxYxyxYxYXYY"},
    {"XyXYXYxYXY", "yxyxYxYXYxYY"},
    {"XyXYXYxYXY", "yyXyXYXyXyxY"},
    {"XyXYXyXyxy", "XYXYxYXYXy"},
    {"XyXYXyXyxy", "XYXyXYXYxY"},
    {"XyXYXyXyxy", "XYxYXYXyXY"},
    {"XyXYXyXyxy", "XyXYXYxYXY"},
    {"XyXYXyXyxy", "yXYXyXYXYxYY"},
    {"XyXYXyXyxy", "yXYxYXYXyXYY"},
    {"XyXYXyXyxy", "yXyXYXYxYXYY"},
    {"XyXYXyXyxy", "yxYXYXyXYXYY"},
    {"XyXYXyXyxy", "yyXyxyxYxyxY"},
    {"XyXYXyXyxy", "yyxYxyxyXyxY"},
    {"XyXYXyXyxy", "yyxyXyxyxYxY"},
    {"XyXyxyXyXY", "XYXYxYXYXy"},
    {"XyXyxyXyXY", "XYXyXYXYxY"},
    {"XyXyxyXyXY", "XYxYXYXyXY"},
    {"XyXyxyXyXY", "XyXYXYxYXY"},
    {"XyXyxyXyXY", "yXYXyXYXYxYY"},
    {"XyXyxyXyXY", "yXYxYXYXyXYY"},
    {"XyXyxyXyXY", "yXyXYXYxYXYY"},
    {"XyXyxyXyXY", "yxYXYXyXYXYY"},
    {"XyXyxyXyXY", "yyXyxyxYxyxY"},
    {"XyXyxyXyXY", "yyxYxyxyXyxY"},
    {"XyXyxyXyXY", "yyxyXyxyxYxY"},
    {"XyxyXyXYXy", "XYXYxYXYXy"},
    {"XyxyXyXYXy", "XYxYXYXyXY"},
    {"XyxyXyXYXy", "yXYXyXYXYxYY"},
    {"XyxyXyXYXy", "yXYxYXYXyXYY"},
    {"XyxyXyXYXy", "yXyXYXYxYXYY"},
    {"XyxyXyXYXy", "yxYXYXyXYXYY"},
    {"XyxyXyXYXy", "yyXyxyxYxyxY"},
    {"XyxyXyXYXy", "yyxYxyxyXyxY"},
    {"XyxyXyXYXy", "yyxyXyxyxYxY"},
    {"XyxyxYxyxy", "XYXyXyxyXy"},
    {"XyxyxYxyxy", "XyXYXyXyxy"},
    {"XyxyxYxyxy", "XyXyxyXyXY"},
    {"XyxyxYxyxy", "XyxyXyXYXy"},
    {"XyxyxYxyxy", "yXYxYxyxYxYY"},
    {"XyxyxYxyxy", "yxYXYxYxyxYY"},
    {"XyxyxYxyxy", "yxYxyxYxYXYY"},
    {"XyxyxYxyxy", "yxyxYxYXYxYY"},
    {"XyxyxYxyxy", "yyXyXYXyXyxY"},
    {"xYxYXYxYxy", "XYXYxYXYXy"},
    {"xYxYXYxYxy", "XYXyXYXYxY"},
    {"xYxYXYxYxy", "XYxYXYXyXY"},
    {"xYxYXYxYxy", "XyXYXYxYXY"},
    {"xYxYXYxYxy", "XyxyxYxyxy"},
    {"xYxYXYxYxy", "yXYXyXYXYxYY"},
    {"xYxYXYxYxy", "yXYxYXYXyXYY"},
    {"xYxYXYxYxy", "yXyXYXYxYXYY"},
    {"xYxYXYxYxy", "yxYXYXyXYXYY"},
    {"xYxYXYxYxy", "yyXyxyxYxyxY"},
    {"xYxYXYxYxy", "yyxYxyxyXyxY"},
    {"xYxYXYxYxy", "yyxyXyxyxYxY"},
    {"xYxyxyXyxy", "XYXyXyxyXy"},
    {"xYxyxyXyxy", "XYxYxyxYxY"},
    {"xYxyxyXyxy", "XyXYXyXyxy"},
    {"xYxyxyXyxy", "XyXyxyXyXY"},
    {"xYxyxyXyxy", "XyxyXyXYXy"},
    {"xYxyxyXyxy", "xYxYXYxYxy"},
    {"xYxyxyXyxy", "yXYxYxyxYxYY"},
    {"xYxyxyXyxy", "yxYXYxYxyxYY"},
    {"xYxyxyXyxy", "yxYxyxYxYXYY"},
    {"xYxyxyXyxy", "yxyxYxYXYxYY"},
    {"xYxyxyXyxy", "yyXyXYXyXyxY"},
    {"xxyxYxxYxy", "XyxyxYxyxy"},
    {"xyXyxyxYxy", "XYXyXyxyXy"},
    {"xyXyxyxYxy", "XYxYxyxYxY"},
    {"xyXyxyxYxy", "XyXYXyXyxy"},
    {"xyXyxyxYxy", "XyXyxyXyXY"},
    {"xyXyxyxYxy", "XyxyXyXYXy"},
    {"xyXyxyxYxy", "xYxYXYxYxy"},
    {"xyXyxyxYxy", "yXYxYxyxYxYY"},
    {"xyXyxyxYxy", "yxYXYxYxyxYY"},
    {"xyXyxyxYxy", "yxYxyxYxYXYY"},
    {"xyXyxyxYxy", "yxyxYxYXYxYY"},
    {"xyXyxyxYxy", "yyXyXYXyXyxY"},
    {"xyxYxyxyXy", "XYXyXyxyXy"},
    {"xyxYxyxyXy", "XYxYxyxYxY"},
    {"xyxYxyxyXy", "XyXYXyXyxy"},
    {"xyxYxyxyXy", "XyXyxyXyXY"},
    {"xyxYxyxyXy", "XyxyXyXYXy"},
    {"xyxYxyxyXy", "xYxYXYxYxy"},
    {"xyxYxyxyXy", "yXYxYxyxYxYY"},
    {"xyxYxyxyXy", "yxYXYxYxyxYY"},
    {"xyxYxyxyXy", "yxYxyxYxYXYY"},
    {"xyxYxyxyXy", "yxyxYxYXYxYY"},
    {"xyxYxyxyXy", "yyXyXYXyXyxY"}
};

extern double g_cosh_marg_upper_bound;
extern double g_cosh_marg_lower_bound;
extern double g_cosh_r_bound;
extern bool g_symmetric;

int num_bound_tests = 4;
int relator_depth= 12;

int TestCollection::size()
{
  return num_bound_tests + pair_vector.size();
}


box_state TestCollection::evaluate_approx(word_pair pair, const Box& box)
{
  Params<Complex> p = box.center();
  if (pair.second.length() == 0) {
    string word = pair.first;
    SL2<Complex> w = construct_word(word, p);
    if (not_identity(w)) {
      if (move_less_than_marg(w, p)) {
        return maybe_killed_center;
      }
      if (moves_x_axis_too_close_to_y(w,p)) {
        return maybe_killed_center;
      }
      if (moves_y_axis_too_close_to_x(w,p)) {
        return maybe_killed_center;
      }
      if (wx_hits_sym_axis(w,p)) {
        return maybe_killed_center;
      }
      if (wy_hits_sym_axis(w,p)) {
        return maybe_killed_center;
      }
    }
    if (y_power(word) > 0) {
      string word_x = x_strip(word);
      SL2<Complex> w_x;
      if (word_x != word) {
        w_x = construct_word(word_x, p);
      } else {
        w_x = w;
      }
      if (inside_var_nbd_x(w_x, p)) {
        if (syllables(word_x) < 4) {
          return maybe_killed_center;
        }
        if (cant_fix_x_axis(w_x,p)) {
          return maybe_killed_center;
        }
        if (non_cyclic_power(w_x, box.x_center())) {
          return maybe_killed_center;
        }
      }
    }
    if (x_power(word) > 0) {
      string word_y = y_strip(word);
      SL2<Complex> w_y;
      if (word_y != word) {
        w_y = construct_word(word_y, p);
      } else {
        w_y = w;
      }
      if (inside_var_nbd_y(w_y, p)) {
        if (syllables(word_y) < 4) {
          return maybe_killed_center;
        }
        if (cant_fix_y_axis(w_y,p)) {
          return maybe_killed_center;
        }
        if (non_cyclic_power(w_y, box.y_center())) {
          return maybe_killed_center;
        }
      }
    }
  } else {
    return maybe_killed_center;
    SL2<Complex> w1 = construct_word(pair.first, p);
    SL2<Complex> w2 = construct_word(pair.second,p);
    if (margulis_smaller_than_xy(w1, w2, p)) {
      return maybe_killed_center;
    }
    if (inside_var_nbd(w1, w2)) {
      return maybe_killed_center;
    }
  }
  return open;
}

TestResult TestCollection::evaluate_vol3(Box& box) {
  TestResult result = {-1, open, word_pair()};
  Params<AJCC> p = box.cover();
  for (auto wp : vol3_rels) {  
    string first = proven_identity(wp.first, p);
    if (first == wp.first) {
      string second = proven_identity(wp.second, p);
      if (second == wp.second) {
          result.state = proven_vol3; 
          result.words.first.assign(first);
          result.words.second.assign(second);
          return result;
      }
    }
  }
  return result;
}

TestResult TestCollection::evaluate_qrs(Box& box) {
  TestResult result = {-1, open, word_pair()};
  Params<AJCC> p = box.cover();
  for (auto word : box.qr.word_classes()) {  
    result.state = open_with_qr;
    string proven = proven_identity(word, p);
    if (relator_test->is_good(proven)) {
      if (box.name.length() > relator_depth) {
        if (relator_test->is_sym(proven)) {
          result.state = killed_via_sym;
        } else {
          SL2<AJCC> w = construct_word(proven, p);
          if (inside_var_nbd_y(w, p)) {
            if (syllables(proven) < 5) {
              result.state = killed_impossible_relator;
            } else if (cant_fix_y_axis(w, p)) {
              if (g_debug) {
                fprintf(stderr, "Word %s is\n", word.c_str());
                print_SL2(w);
              }
              result.state = killed_y_hits_y;
            } else if (non_cyclic_power(w, box.y_cover())) {
              result.state = killed_y_hits_y;
            }
            // result.state = proven_relator;
            // fprintf(stderr, "Proven relator %s\n", proven.c_str());
            // print_SL2(construct_word(proven, p));
          }
        }
      }
    }
    vector<string> required;
    if (relator_test->is_impossible(proven, required)) {
      fprintf(stderr, "impossible with %s and size %d\n", proven.c_str(), required.size());
      if (required.size() == 0) {
        result.state = killed_impossible_relator;
        fprintf(stderr, "Impossible with %s\n", proven.c_str());
        fprintf(stderr, "Impossible relator:\n");
        print_SL2(construct_word(proven, p));
      } else {
        for (auto req : required) {
          SL2<AJCC> w_req = construct_word(req, p);
          if (not_identity(w_req)) {
            result.state = killed_impossible_relator;
            break;
          }
        }
      }
    }
    if (result.state != open &&
        result.state != open_with_qr) {
      result.words.first.assign(proven);
      break;
    }
  }
  return result;
}

TestResult TestCollection::evaluate_AJCC(word_pair& pair, Box& box)
{
  TestResult result = {-1, open, pair};
  Params<AJCC> p = box.cover();
  if (pair.second.length() == 0) {
    string word = pair.first;
    SL2<AJCC> w = construct_word(word, p);
    if (g_debug) {
      fprintf(stderr, "Testing word %s\n", word.c_str());
      print_SL2(w);
    }
    if (not_identity(w) &&
        move_less_than_marg(w, p) &&
        ((x_power(word) == 0 || y_power(word) == 0) ||
          does_not_fix_sym_axis(w))) {
      result.state = killed_move;
      return result;
    }
    if (wy_hits_sym_axis(w,p)) {
      result.state = killed_w_ay_hits_sym_axis;
      return result;
    }
    if (x_power(word) > 0) {
      // string word_yr = y_rstrip(word);
      string word_yr = word;
      SL2<AJCC> w_yr;
      if (word_yr != word) {
        w_yr = construct_word(word_yr, p);
      } else {
        w_yr = w;
      }
      if (moves_y_axis_too_close_to_x(w_yr,p)) {
        box.qr.get_name(word_yr);
        if (moved_y_axis_not_x_axis(w_yr, p)) {
          result.words.first.assign(word_yr);
          result.state = killed_y_hits_x;
          return result;
        }
      }
    }
    if (x_power(word) > 0) {
      // string word_y = y_strip(word);
      string word_y = word;
      SL2<AJCC> w_y;
      if (word_y != word) {
        w_y = construct_word(word_y, p);
      } else {
        w_y = w;
      }
      if (inside_var_nbd_y(w_y, p)) {
        box.qr.get_name(word_y);
        if (syllables(word_y) < 4 || cant_fix_y_axis(w_y, p)) {
          result.words.first.assign(word_y);
          if (g_debug) {
            fprintf(stderr, "Word %s is\n", word_y.c_str());
            print_SL2(w_y);
          }
          result.state = killed_y_hits_y;
          return result;
        }
        if (non_cyclic_power(w_y, box.y_cover())) {
          result.words.first.assign(word_y);
          result.state = killed_y_hits_y;
          return result;
        }
      }
    }
    result = evaluate_qrs(box);
    if (result.state != open && 
        result.state != open_with_qr) {
      return result;
    }
    // test vol3
    result = evaluate_vol3(box);
    if (result.state != open && 
        result.state != open_with_qr) {
      return result;
    }
  } else {
    result = evaluate_vol3(box);
    if (result.state != open && 
        result.state != open_with_qr) {
      return result;
    }
//    SL2<AJCC> w1 = construct_word(pair.first, p);
//    SL2<AJCC> w2 = construct_word(pair.second,p);
//    if (margulis_smaller_than_xy(w1, w2, p)) {
//      result.words = pair;
//      return result;
//    }
  }
  if (box.qr.word_classes().size() > 0) {
    result.state = open_with_qr;
  } else {
    result.state = open;
  }
  return result;
}

box_state check_bounds_center(bool ans) {
  if (ans) return out_of_bounds_center;
  else return open;
}

TestResult check_bounds(bool ans, TestResult& result) {
  if (ans) {
    result.state = killed_bounds;
  }
  return result;
}

box_state TestCollection::evaluate_center(int index, Box& box)
{
  Params<Complex> center = box.center();
  Complex one(1);
  switch(index) {
    case 0:	{ // 1.0052 < cosh(0.104) <= cosh(mu) <= 0.
              return check_bounds_center(absUB(center.coshmu) < g_cosh_marg_lower_bound ||
                  absLB(center.coshmu) > g_cosh_marg_upper_bound);
            }
    case 1: {
              return check_bounds_center(absLB(center.twocoshreD2) > g_cosh_r_bound * 2 ||
                      strictly_pos(re(-center.sinhL2)) || strictly_pos(re(-center.sinhD2)) ||
                      strictly_pos(one - center.coshreD) || strictly_pos(one - center.coshreL));
            }
    case 2: { // Meyerhoff tube bound. Check if embeded tube radius is more than rad + marg/2
              SL2<Complex> x = construct_x(center);
              SL2<Complex> y = construct_y(center);
              Complex four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, center);
              Complex four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, center);
              return check_bounds_center(
                  meyerhoff_k_test(center.coshreL, center.cosimL, four_cosh_x_tube_UB) ||
                  meyerhoff_k_test(center.coshreL, center.cosimL, four_cosh_y_tube_UB));
            }
    case 3: { // 4.26 in bilipschitz paper
              Complex cosh_mu_LB = cosh_marg_lower_bound(center.twosinhreD2);
              return check_bounds_center(strictly_pos(cosh_mu_LB - center.coshmu));
            }
   default:
            return evaluate_approx(pair_vector[index - num_bound_tests], box);
  }
}

TestResult TestCollection::evaluate_box(int index, Box& box)
{
  Params<AJCC> cover = box.cover();
  TestResult result = {index, open, word_pair()};
  AJCC one(1);
  switch(index) {
    case 0:	{ // 1.0052 < cosh(0.104) <= cosh(mu) <= 0.
              if (g_debug) {
                print_type(cover.coshmu);
              }
              return check_bounds(absUB(cover.coshmu) < g_cosh_marg_lower_bound ||
                  absLB(cover.coshmu) > g_cosh_marg_upper_bound, result);
            }
    case 1: {
              Params<Complex> nearer = box.nearer();
              return check_bounds(absLB(cover.twocoshreD2) > g_cosh_r_bound * 2 ||
                      strictly_pos(re(-nearer.sinhL2)) ||
                      strictly_pos(re(-nearer.sinhD2)) ||
                      strictly_pos(one - cover.coshreD) ||
                      strictly_pos(one - cover.coshreL),
                      result);
            }
    case 2: { // Meyerhoff tube bound. Check if embeded tube radius is 
              // more than rad + marg/2
              SL2<AJCC> x = construct_x(cover);
              SL2<AJCC> y = construct_y(cover);
              AJCC four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, cover);
              return check_bounds(
                  meyerhoff_k_test(cover.coshreL, cover.cosimL, four_cosh_y_tube_UB),
                  result);
            }
    case 3: { // 4.26 in bilipschitz paper
              return check_bounds(false, result);
              // AJCC cosh_mu_LB = cosh_marg_lower_bound(cover.twosinhreD2);
              // return check_bounds(strictly_pos(cosh_mu_LB - cover.coshmu), result);
            }
    default:
            return evaluate_AJCC(pair_vector[index - num_bound_tests], box);
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

word_pair TestCollection::get_pair(int index)
{
  if (index < num_bound_tests) {
    return word_pair();
  } else {
    word_pair p = pair_vector[index - num_bound_tests];
    return p;
  }
}

word_pair TestCollection::parse_word_pair(string buf)
{
  size_t start = buf.find('(');
  size_t comma = buf.find(',');
  string first;
  string second;
  // fprintf(stderr, "Adding test: %s\n", buf.c_str());
  if (start != string::npos) {
    size_t end = comma;
    if (comma == string::npos) {
      end = buf.find(')');
    }
    if (end == string::npos) {
      return word_pair();
    }
    first = buf.substr(start + 1, end - start - 1);
  } else {
    return word_pair();
  }
  if (comma != string::npos) {
    size_t end = buf.find(')');
    second = buf.substr(comma + 1, end - comma - 1);
  } else {
    return word_pair();
  }
  word_pair p(first, second);
  return p;
}


int TestCollection::add(string buf)
{
  return add(parse_word_pair(buf));
}

int TestCollection::add(word_pair p) {
  map< word_pair,int >::iterator it = pair_index.find(p);
  if (it == pair_index.end()) {
    // fprintf(stderr, "Adding test: (%s,%s)\n", p.first.c_str(), p.second.c_str());
    pair_index[p] = pair_vector.size();
    pair_vector.push_back(p);
    seen_words[p.first] = 1;
    seen_words[p.second] = 1;
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

void TestCollection::load_relator_test(
    const char* impos_path, const char* bad_rel_path)
{
  relator_test = RelatorTest::create(impos_path, bad_rel_path);
}

/*
void print_debug(SL2<AJCC>& w, Params<AJCC>& p, box_state state) {
  if (state == killed_x_hits_y) { 
    fprintf(stderr, "******* MOVED X TOO CLOSE TO Y *********\n");
    AJ diff = p.coshdxdy * 4 - four_cosh_dist_ay_wax(w, p);
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshdxdy * 4);
    print_type("4coshd(dist(y-axis, w(x-axis))):", four_cosh_dist_ay_wax(w, p));
    AJ z = ((w.a * w.a) * p.expmdx - (w.b * w.b) * p.expdx ) * p.expmdyf +
      ((w.d * w.d) * p.expdx  - (w.c * w.c) * p.expmdx) * p.expdyf;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    AJ fsp2sq = four_sinh_perp2_sq_ay_wax(w, p);
    print_type("4 sihn^2(perp/2):", fsp2sq);
    fprintf(stderr, "****************************************\n");
  }
  if (state == killed_y_hits_x) {
    fprintf(stderr, "******* MOVED Y TOO CLOSE TO X: %s *********\n", word.c_str());
    AJCC diff = p.coshreD * 4 - four_cosh_dist_ax_way(w, p);
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshreD * 4);
    print_type("4coshd(dist(x-axis, w(y-axis))):", four_cosh_dist_ax_way(w, p));
    AJCC z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
           ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    AJCC fsp2sq = four_sinh_perp2_sq_ax_way(w, p);
    print_type("4 sihn^2(perp/2):", fsp2sq);
    fprintf(stderr, "****************************************\n");
  }
  if (state == killed_x_hits_x) {
    fprintf(stderr, "********** KILLED  ***********\n");
    fprintf(stderr, "Word %s must but doesn't fix x-axis\n", word_x.c_str());
    print_SL2(w_x);
    fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
    fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w_x, p)),
            absUB(jorgensen_xw(w_x, p)), must_fix_x_axis(w_x, p));
    AJCC diff = p.coshreD * 4 - four_cosh_dist_ax_wax(w_x, p);
    print_type("4 cosh 2 dx:", p.coshreD * 4);
    print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w_x, p));
    print_type("diff:", diff);
    fprintf(stderr, "********** CANNOT FIX X AXIS ***********\n");
    AJCC fsp2sq = four_sinh_perp2_sq_ax_wax(w_x, p);
    print_type("4sinh^2(perp/2)", fsp2sq);
    fprintf(stderr, "Can't fix x axis LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr,"Can't fix x axis LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
    fprintf(stderr, "*******************************\n");
  }
  if (state == killed_y_hits_y) {
    fprintf(stderr, "********** DOES NOT COMMUTE  ***********\n");
    fprintf(stderr, "Word %s must fix x-axis but doesn't commute\n", word_x.c_str());
    print_SL2(w_x);
    fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
    fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w_x, p)),
            absUB(jorgensen_xw(w_x, p)), must_fix_x_axis(w_x, p));
    AJCC diff = p.coshreD * 4 - four_cosh_dist_ax_wax(w_x, p);
    print_type("4 cosh 2 dx:", p.coshreD * 4);
    print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w_x, p));
    print_type("diff:", diff);
    fprintf(stderr, "********** DOES NOTE COMMUTE ***********\n");
    SL2<AJCC> commutator = box.x_cover() * w * inverse(w * box.x_cover());
    fprintf(stderr, "commutator\n");
    print_SL2(commutator);
    fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
    fprintf(stderr, "****************************************\n");
  } 
}
*/
