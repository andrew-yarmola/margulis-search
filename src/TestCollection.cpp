#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TestCollection.hh"
// #include "ImpossibleRelations.h"
#include <algorithm>
using namespace std;
// using namespace __gnu_cxx;

extern double g_cosh_marg_upper_bound;
extern double g_cosh_marg_lower_bound;
extern double g_sinh_d_bound; 
extern bool g_symmetric; 

int num_bound_tests = 8;

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
    if (strictly_pos(p.coshlx * 4 - four_cosh_re_length(w))) {
      return bad_length_center;
    }
    if (not_identity(w)) {
      if (move_less_than_marg(w, p)) {
        return bad_move_center;
      }
      if (moves_x_axis_too_close_to_y(w,p)) {
        return x_hits_y_center;
      }
      if (moves_y_axis_too_close_to_x(w,p)) {
        return y_hits_x_center; 
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
          return bad_x_tube_center;
        }
        if (cant_fix_x_axis(w_x,p)) {
          if (g_debug) { 
            fprintf(stderr, "Word doesn't fix axis\n");
            print_SL2(w_x);
            fprintf(stderr, "x\n");
            print_SL2(box.x_center());
            fprintf(stderr, "&&&&&&&&&&&&&&&&&&&&\n");
          }
          return bad_x_tube_center;
        }
        if (non_cylic_power(w_x, box.x_center())) {
          return bad_lox_x_center;
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
          return bad_y_tube_center;
        }
        if (cant_fix_y_axis(w_y,p)) {
          return bad_y_tube_center;
        }
        if (non_cylic_power(w_y, box.y_center())) {
          return bad_lox_y_center;
        }
      }
    }
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
  if (g_debug) {
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n     Word Pair: %s and %s\n +++++++++++++++++++++++++++++++++++++++++++\n", pair.first.c_str(), pair.second.c_str());
  }
  Params<AJ> p = box.cover();
  vector<string> required;
  if (pair.second.length() == 0) {
    string word = pair.first;
    SL2<AJ> w = construct_word(word, p);
//    if (strictly_pos(p.coshlx * 4 - four_cosh_re_length(w))) {
//      return bad_length;
//    }
    if (not_identity(w) && move_less_than_marg(w, p)) {
      return killed_move;
    }
    string word_xr = x_rstrip(word);
    SL2<AJ> w_xr;
    if (word_xr != word) {
      w_xr = construct_word(word_xr, p);
    } else {
      w_xr = w;
    }
    if (moves_x_axis_too_close_to_y(w_xr,p)) {
      if (moved_x_axis_not_y_axis(w_xr, p)) {
        if (g_debug) {
          fprintf(stderr, "******* MOVES X TOO CLOSE TO Y *********\n");
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
        return killed_x_hits_y;
      }
      string proven = proven_identity(word_xr, p);
      if (proven.length() > 0) {
        if(impossible->is_impossible(proven, required)) {
          aux_word.assign(proven);
          return killed_failed_qr;
        } else { //HACK
          return killed_failed_qr;
        }
      }
      new_qrs.push_back(word_xr);
      return var_x_hits_y; 
    }
    string word_yr = y_rstrip(word);
    SL2<AJ> w_yr;
    if (word_yr != word) {
      w_yr = construct_word(word_yr, p);
    } else {
      w_yr = w;
    }
    if (moves_y_axis_too_close_to_x(w_yr,p)) {
      if (moved_y_axis_not_x_axis(w_yr, p)) {
        if (g_debug) {
          fprintf(stderr, "******* MOVES Y TOO CLOSE TO X: %s *********\n", word.c_str());
          AJ diff = p.coshdxdy * 4 - four_cosh_dist_ax_way(w, p);
          print_SL2(w);
          print_type("4cosh(dx+dy):", p.coshdxdy * 4);
          print_type("4coshd(dist(x-axis, w(y-axis))):", four_cosh_dist_ax_way(w, p)); 
          AJ z = ((w.a * w.a) * p.expdyf  - (w.b * w.b) * p.expmdyf) * p.expdx +
                 ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf ) * p.expmdx;
          print_type("4 sinh^2(dist/2) + 2:", z);
          print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
          print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
          print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
          print_type("diff:", diff);
          fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
          AJ fsp2sq = four_sinh_perp2_sq_ax_way(w, p);
          print_type("4 sihn^2(perp/2):", fsp2sq);
          fprintf(stderr, "****************************************\n");
        }
        return killed_y_hits_x;
      }
      string proven = proven_identity(word_yr, p);
      if (proven.length() > 0) {
        if(impossible->is_impossible(proven, required)) {
          aux_word.assign(proven);
          return killed_failed_qr;
        } else { //HACK
          return killed_failed_qr;
        }
      }
      new_qrs.push_back(word_yr);
      return var_y_hits_x;
    }
    if (y_power(word) > 0) {
      string word_x = x_strip(word);
      SL2<AJ> w_x;
      if (word_x != word) {
        w_x = construct_word(word_x, p);
      } else {
        w_x = w;
      }
      if (inside_var_nbd_x(w_x, p)) {
        if (syllables(word_x) < 4) {
          return killed_x_tube;
        }
        if (cant_fix_x_axis(w_x, p)) {
          if (g_debug) {
            fprintf(stderr, "********** KILLED  ***********\n");
            fprintf(stderr, "Word %s must but doesn't fix x-axis\n", word_x.c_str());
            print_SL2(w_x);
            fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
            fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w_x, p)),
                    absUB(jorgensen_xw(w_x, p)), must_fix_x_axis(w_x, p));
            AJ diff = p.cosh2dx * 4 - four_cosh_dist_ax_wax(w_x, p);
            print_type("4 cosh 2 dx:", p.cosh2dx * 4);
            print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w_x, p));
            print_type("diff:", diff);
            fprintf(stderr, "********** CANNOT FIX X AXIS ***********\n");
            AJ fsp2sq = four_sinh_perp2_sq_ax_wax(w_x, p);
            print_type("4sinh^2(perp/2)", fsp2sq);
            fprintf(stderr, "Can't fix x axis LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
            fprintf(stderr,"Can't fix x axis LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
            fprintf(stderr, "*******************************\n");
           }
          return killed_x_tube;
        }
        if (non_cylic_power(w_x, box.x_cover())) {
          if (g_debug) {
            fprintf(stderr, "********** DOES NOT COMMUTE  ***********\n");
            fprintf(stderr, "Word %s must fix x-axis but doesn't commute\n", word_x.c_str());
            print_SL2(w_x);
            fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
            fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w_x, p)),
                    absUB(jorgensen_xw(w_x, p)), must_fix_x_axis(w_x, p));
            AJ diff = p.cosh2dx * 4 - four_cosh_dist_ax_wax(w_x, p);
            print_type("4 cosh 2 dx:", p.cosh2dx * 4);
            print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w_x, p));
            print_type("diff:", diff);
            fprintf(stderr, "********** DOES NOTE COMMUTE ***********\n");
            SL2<AJ> commutator = box.x_cover() * w * inverse(w * box.x_cover());
            fprintf(stderr, "commutator\n");
            print_SL2(commutator);
            fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
            fprintf(stderr, "****************************************\n");
          }
          return killed_lox_not_x_power;
        }
        string proven = proven_identity(word_x, p);
        if (proven.length() > 0) {
          if(impossible->is_impossible(proven, required)) {
            aux_word.assign(proven);
            return killed_failed_qr;
          } else { //HACK
            return killed_failed_qr;
          }
        }
        new_qrs.push_back(word_x);
        return variety_nbd_x;
      }
    }
    if (x_power(word) > 0) {
      string word_y = y_strip(word);
      SL2<AJ> w_y;
      if (word_y != word) {
        w_y = construct_word(word_y, p);
      } else {
        w_y = w;
      }
      if (inside_var_nbd_y(w_y, p)) {
        if (syllables(word_y) < 4) {
          return killed_y_tube;
        }
        if (cant_fix_y_axis(w_y, p)) {
          return killed_y_tube;
        }
        if (non_cylic_power(w_y, box.y_cover())) {
          return killed_lox_not_y_power;
        }
        string proven = proven_identity(word_y, p);
        if (proven.length() > 0) {
          if(impossible->is_impossible(proven, required)) {
            aux_word.assign(proven);
            return killed_failed_qr;
          } else { //HACK
            return killed_failed_qr;
          }
        }
        new_qrs.push_back(word_y);
        return variety_nbd_y;
      }
    }
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
  if (new_qrs.size() > 0) {
    return open_with_qr;
  } else {
    return open;
  }
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
              return check_bounds_center(absUB(center.coshlx) < 1 || absUB(center.coshly) < 1);
            }
    case 4: { // Meyerhoff tube bound. Check if embeded tube radius is more than rad + marg/2 
              SL2<Complex> x = construct_x(center);
              SL2<Complex> y = construct_y(center);
              Complex four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, center);
              Complex four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, center);
              return check_bounds_center(
                  meyerhoff_k_test(center.coshlx, center.costx, four_cosh_x_tube_UB) || 
                  meyerhoff_k_test(center.coshly, center.costy, four_cosh_y_tube_UB));
            }
    case 5: { // diagonals
              return check_bounds_center(false);
            }
    case 6: { // 4.26 in bilipschitz paper
              Complex cosh_mu_LB_x = cosh_marg_lower_bound(center.sinhdx);
              Complex cosh_mu_LB_y = cosh_marg_lower_bound(center.sinhdy); // overkill in symmertic
              return check_bounds_center(
                  strictly_pos(cosh_mu_LB_x - center.coshmu) ||
                  strictly_pos(cosh_mu_LB_y - center.coshmu));
                    
            }
    case 7: { // x is shortest
              return check_bounds_center(strictly_pos(center.coshlx - center.coshly));
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
              return check_bounds(absUB(cover.coshlx) < 1 || absUB(cover.coshly) < 1);
            }
    case 4: { // Meyerhoff tube bound. Check if embeded tube radius is more than rad + marg/2 
              // fprintf(stderr, "%s", box.desc().c_str());
              SL2<AJ> x = construct_x(cover);
              SL2<AJ> y = construct_y(cover);
              AJ four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, cover);
              AJ four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, cover);
              return check_bounds(meyerhoff_k_test(cover.coshlx, cover.costx, four_cosh_x_tube_UB) || 
                  meyerhoff_k_test(cover.coshly, cover.costy, four_cosh_y_tube_UB));
            }
    case 5: { // diagonals
              return check_bounds(false);
            }
    case 6: { // 4.26 in bilipschitz paper
              AJ cosh_mu_LB_x = cosh_marg_lower_bound(cover.sinhdx);
              AJ cosh_mu_LB_y = cosh_marg_lower_bound(cover.sinhdy); // overkill in symmertic
              return check_bounds(
                  strictly_pos(cosh_mu_LB_x - cover.coshmu) ||
                  strictly_pos(cosh_mu_LB_y - cover.coshmu));
                    
            }
    case 7: { // x is shortest
              return check_bounds(strictly_pos(cover.coshlx - cover.coshly));
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

