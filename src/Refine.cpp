#include "Refine.hh"
#include "TestCollection.hh"
#include "TubeSearch.hh"
#include "QuasiRelators.h"

using namespace std;

typedef vector< vector< box_state > > TestHistory;

Options g_options;
TestCollection g_tests;
int g_boxes_visited = 0;

extern double g_cosh_marg_upper_bound;
extern double g_cosh_marg_lower_bound;
extern double g_sinh_d_bound; 

extern int num_bound_tests;

extern bool g_debug;

unordered_map<string, SL2<AJCC> > short_words_cache;

// TODO
//void parse_result(box_state result, Box& box, PartialTree& t) {
//}

bool refine_recursive(Box box, PartialTree& t, int depth, TestHistory& history, vector< Box >& place, int newDepth, int& searched_depth)
{
  place.push_back(box);
  int old_test_index = t.test_index;
  vector<string> new_qrs;
  short_words_cache.clear();

  string aux_word;
  if (t.test_index >= 0) {
    box_state result = g_tests.evaluate_box(t.test_index, box, aux_word, new_qrs, short_words_cache);
    if (result != open && result != open_with_qr) {
      t.aux_word.assign(aux_word);
      t.test_result = result;
      return true;
    } else if (result == open_with_qr) {
      for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
        box.qr.get_name(*it); // Also adds qr to the box's list
      }
    } else { 
      fprintf(stderr, "FAILED to eliminate %s with test %s with result %d\n", box.name.c_str(), g_tests.get_name(t.test_index).c_str(), result);
    }
  }

  if (t.test_index == -2 && !g_options.fill_holes) {
    return true;
  }

  // Check if the box is now small enough that some former qrs actually kill it
  Params<AJCC> p = box.cover();
  vector<string> quasi_relators = box.qr.word_classes();
  for (vector<string>::iterator it = quasi_relators.begin(); it != quasi_relators.end(); ++it) {
    SL2<AJCC> w = construct_word(*it, p, short_words_cache); 
    if (wx_hits_elliptic_axis(w,p)) {
        t.aux_word.assign(*it);
        t.aux_result = wx_hits_elliptic;
        t.test_result = killed_failed_qr;
        return true;
    }
    if (wy_hits_elliptic_axis(w,p)) {
        t.aux_word.assign(*it);
        t.aux_result = wy_hits_elliptic;
        t.test_result = killed_failed_qr;
        return true;
    }
    if (moves_x_axis_too_close_to_y(w,p) &&
        moved_x_axis_not_y_axis(w, p)) {
        t.aux_word.assign(*it);
        t.aux_result = killed_x_hits_y;
        t.test_result = killed_failed_qr;
        return true;
    }
    if (moves_y_axis_too_close_to_x(w,p) &&
        moved_y_axis_not_x_axis(w, p)) {
        t.aux_word.assign(*it);
        t.aux_result = killed_y_hits_x;
        t.test_result = killed_failed_qr;
        return true;
    }
    if (inside_var_nbd_x(w,p)) {
      if (cant_fix_x_axis(w,p)) {
        t.aux_word.assign(*it);
        t.aux_result = killed_x_tube;
        t.test_result = killed_failed_qr;
        return true;
      } 
      if (non_cylic_power(w, box.x_cover())) {
        t.aux_word.assign(*it);
        t.aux_result = killed_lox_not_x_power;
        t.test_result = killed_failed_qr;
        return true;
      }
    }
    if (inside_var_nbd_y(w,p)) {
      if (cant_fix_y_axis(w,p)) {
        t.aux_word.assign(*it);
        t.aux_result = killed_y_tube;
        t.test_result = killed_failed_qr;
        return true;
      } 
      if (non_cylic_power(w, box.y_cover())) {
        t.aux_word.assign(*it);
        t.aux_result = killed_lox_not_y_power;
        t.test_result = killed_failed_qr;
        return true;
      }
    }
    vector<string> required;
    string proven = proven_identity(*it, p);
    if (proven_is_good(proven, box)) {
      if (g_tests.impossible->is_impossible(proven, required)) {
        t.aux_word.assign(proven);
        t.aux_result = open;
        t.test_result = killed_failed_qr;
        return true;
      } else {
        t.aux_word.assign(proven);
        t.aux_result = proven_relator;
        t.test_result = killed_failed_qr;
        return true;
      }
    }
  }

  if (g_options.improve_tree || !t.l_child) {
    for (int i = 0; i < g_tests.size(); ++i) {
      if (i >= num_bound_tests && depth % 12 == 0) {
        break;
      } 
      vector<box_state>& th = history[i];
      while (th.size() <= depth) {
        box_state result = g_tests.evaluate_center(i, place[th.size()]);
        th.push_back(result);
      }
      bool do_eval = true;
      int s = th.size();
      for (int j = 1; j <= min(s, 7); j++) {
        do_eval = do_eval && th[th.size() - j] != open;
      }
      if (do_eval) {
        new_qrs.clear();
        box_state result = g_tests.evaluate_box(i, box, aux_word, new_qrs, short_words_cache);
        switch (result) {
          case killed_bounds :
          case killed_only_elliptic : 
          case killed_x_hits_y :
          case killed_y_hits_x :
          case killed_x_tube :
          case killed_y_tube :
          case killed_lox_not_x_power : 
          case killed_lox_not_y_power :
          case killed_move :
          case killed_marg :
          case variety_nbd_x :
          case variety_nbd_y :
          case variety_nbd : 
          case bad_length : 
          case wx_hits_elliptic :
          case wy_hits_elliptic :
          case var_x_hits_y :
          case var_y_hits_x : {
            t.test_index = i;
            t.test_result = result;
            return true;
          }
          case killed_failed_qr : {
              t.test_index = i;
              t.aux_word.assign(aux_word);
              t.aux_result = proven_relator;
              t.test_result = killed_failed_qr;
              return true;
          }
          case open_with_qr : {
            for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
              box.qr.get_name(*it); // Also adds qr to the box's list
            }
            break;
          }
          default : {
            continue;
          }
        }
      }
    }
  }

  if (g_options.word_search_depth > 0 && depth > 0 && (g_options.improve_tree || !t.l_child) && box.name.length() > 60 && depth % g_options.word_search_depth == 0) {
    Box& search_place = box;
    vector<word_pair> search_pairs_v1 = find_pairs(search_place.center(), vector<string>(), 1, g_options.max_word_length, box.qr.word_classes());
    vector<word_pair> search_pairs_v2 = find_words_v2(search_place.center(), 1, 7, box.qr.word_classes(), map<string, int>());
    vector<word_pair> search_pairs;
    search_pairs.insert(search_pairs.end(), search_pairs_v1.begin(), search_pairs_v1.end());
    search_pairs.insert(search_pairs.end(), search_pairs_v2.begin(), search_pairs_v2.end());
    while (search_pairs.size() > 0) {
      word_pair new_pair = search_pairs.back();

      int old_size = g_tests.size();
      int new_index = g_tests.add(new_pair);
      history.resize(g_tests.size());

      search_pairs.pop_back();

      if (old_size < g_tests.size()) {
        fprintf(stderr, "search (%s) found (%s,%s) at (%s)\n",
                search_place.qr.desc(box.cover()).c_str(), new_pair.first.c_str(), new_pair.second.c_str(), search_place.name.c_str());

        new_qrs.clear();
        box_state result = g_tests.evaluate_box(new_index, box, aux_word, new_qrs, short_words_cache);

        switch (result) {
          case killed_bounds :
          case killed_only_elliptic : 
          case killed_x_hits_y :
          case killed_y_hits_x :
          case killed_x_tube :
          case killed_y_tube :
          case killed_lox_not_x_power : 
          case killed_lox_not_y_power :
          case killed_move :
          case killed_marg :
          case variety_nbd_x :
          case variety_nbd_y :
          case variety_nbd : 
          case bad_length : 
          case wx_hits_elliptic :
          case wy_hits_elliptic :
          case var_x_hits_y :
          case var_y_hits_x : {
            t.test_index = new_index;
            t.test_result = result;
            return true;
          }
          case open_with_qr : {
            for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
              box.qr.get_name(*it); // Also adds qr to the box's list
            }
            break;
          }
          default : {
            continue;
          }
        }
      }
    }
  }

  t.test_index = -1;

  if (!t.l_child) {
    if (depth >= g_options.max_depth || ++g_boxes_visited >= g_options.max_size || ++newDepth > g_options.invent_depth) {
      fprintf(stderr, "HOLE %s (%s)\n", box.name.c_str(), box.qr.desc(box.cover()).c_str());
      return false;
    }
    t.l_child = new PartialTree();
    t.r_child = new PartialTree();
  }

  bool is_complete = true;

  is_complete = refine_recursive(box.child(0), *t.l_child, depth+1, history, place, newDepth, searched_depth) && is_complete;
  if (place.size() > depth+1)
    place.resize(depth+1);
  for (int i = 0; i < g_tests.size(); ++i) {
    if (history[i].size() > depth)
      history[i].resize(depth);
  }
  if (searched_depth > depth)
    searched_depth = depth;
  if (is_complete || depth < g_options.truncate_depth)
    is_complete = refine_recursive(box.child(1), *t.r_child, depth+1, history, place, newDepth, searched_depth) && is_complete;
  if (old_test_index >= 0 && t.test_index != old_test_index) {
    fprintf(stderr, "invalid box %s(%s) %d %s\n", g_tests.get_name(old_test_index).c_str(), box.name.c_str(),
      tree_size(t), is_complete ? "Patched" : "Unpatched");
  }
  if (!is_complete && depth >= g_options.truncate_depth) {
    truncate_tree(t);
  }
  return is_complete;
}

void refine_tree(Box box, PartialTree& t)
{
  TestHistory history(g_tests.size());
  vector<Box> place;
  int searched_depth = 0;
  refine_recursive(box, t, 0, history, place, 0, searched_depth);
}

void print_tree(PartialTree& t)
{
    char type = 'F';
    word_pair p = g_tests.get_pair(t.test_index);
    switch (t.test_result) {
      case open :
      case open_with_qr : {
        if (t.l_child && t.r_child) {
          printf("X\n");
          print_tree(*t.l_child);
          print_tree(*t.r_child);
        } else {
          printf("HOLE (%s)\n", t.qr_desc.c_str());
        }
        return;
      }
      case killed_bounds : {
        printf("%s\n", g_tests.get_name(t.test_index).c_str());
        return;
      }
      case killed_failed_qr : {
        p = word_pair(t.aux_word, "");
        switch (t.aux_result) {
          case wx_hits_elliptic : type = 'n'; break;
          case wy_hits_elliptic : type = 'N'; break;
          case killed_x_hits_y : type = 'a'; break;
          case killed_y_hits_x : type = 'A'; break;
          case killed_x_tube : type = 'x'; break;
          case killed_lox_not_x_power : type = 'p';  break;
          case killed_y_tube : type = 'y'; break;
          case killed_lox_not_y_power : type = 'P'; break;
          case proven_relator : type = 'R'; break;
          default: type = 'K'; break;
        }
        break;
      }
      case killed_only_elliptic : type = 'E'; break; 
      case killed_x_hits_y : {
                              type = 'a';
                              p = word_pair(x_rstrip(p.first), p.second);
                              break;
                             }
      case killed_y_hits_x : {
                              type = 'A';
                              p = word_pair(y_rstrip(p.first), p.second);
                              break;
                             }
      case killed_x_tube : {
                             type = 'x';
                             p = word_pair(x_strip(p.first), p.second);
                             break;
                           }
      case killed_lox_not_x_power : { 
                             type = 'p';
                             p = word_pair(x_strip(p.first), p.second);
                             break;
                           }
      case killed_y_tube : { 
                             type = 'y';
                             p = word_pair(y_strip(p.first), p.second);
                             break;
                           }
      case killed_lox_not_y_power : { 
                             type = 'P';
                             p = word_pair(y_strip(p.first), p.second);
                             break;
                           }
      case killed_move : type = 'm'; break;
      case killed_marg : type = 'M'; break;
      case bad_length : type = 'L'; break; 
      case variety_nbd_x : type = 'v'; break;
      case variety_nbd_y : type = 'V'; break;
      case variety_nbd : type = 'W'; break;
      case wx_hits_elliptic : type = 'n'; break;
      case wy_hits_elliptic : type = 'N'; break;
      case var_x_hits_y : type = 'c'; break;
      case var_y_hits_x : type = 'C'; break;
      default : return;
    }
    
    printf("%c(%s,%s)\n", type, p.first.c_str(), p.second.c_str()); 
}

