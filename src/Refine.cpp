#include "Refine.hh"
#include "TestCollection.hh"
#include "TubeSearch.hh"
#include "QuasiRelators.h"

using namespace std;

typedef vector< vector< box_state > > TestHistory;

Options g_options;
TestCollection g_tests;
int g_boxesVisited = 0;

double g_cosh_marg_upper_bound = 1.2947;
double g_cosh_marg_lower_bound = 1.0054;
double g_sinh_d_bound = 1.3426; 

unordered_map<string, SL2<AJ> > short_words_cache;

bool refine_recursive(Box box, PartialTree& t, int depth, TestHistory& history, vector< Box >& place, int newDepth, int& searched_depth)
{
  //fprintf(stderr, "rr: %s depth %d placeSize %lu\n", box.name.c_str(), depth, place.size());
  place.push_back(box);
  int old_test_index = t.test_index;
  vector<string> new_qrs;
  short_words_cache.clear();

  string aux_word;
  if (t.test_index >= 0) {
//    fprintf(stderr, "********************************* Validation *********************************\n");
//    fprintf(stderr, "%s", box.desc().c_str());
    box_state result = g_tests.evaluate_box(t.test_index, box, aux_word, new_qrs, short_words_cache);
    if (result != open && result != open_with_qr) {
      t.aux_word.assign(aux_word);
      t.test_result = result;
//      fprintf(stderr, "Eliminated %s with test %s with result %d\n", box.name.c_str(), g_tests.get_name(t.test_index), result);
//      fprintf(stderr, "********************************* End Validation *********************************\n");
//      fprintf(stderr, "Test %s kills\n %s", g_tests.get_name(t.test_index).c_str(), box.desc().c_str());
      return true;
    } else if (result == open_with_qr) {
//      fprintf(stderr,"Retested %d, new qrs len %lu\n", result, new_qrs.size());
      for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
        fprintf(stderr, "new quasirelator %s\n", (*it).c_str());
        box.qr.get_name(*it); // Also adds qr to the box's list
      }
      t.qr_desc = box.qr.min_pow_desc();
    } else { 
      fprintf(stderr, "FAILED to eliminate %s with test %s with result %d\n", box.name.c_str(), g_tests.get_name(t.test_index).c_str(), result);
    }
//    fprintf(stderr, "********************************* End Validation *********************************\n");
  }

  if (t.test_index == -2 && !g_options.fill_holes) {
    return true;
  }

  // Check if the box is now small enough that some former qrs actually kill it
  Params<AJ> p = box.cover();
  vector<string> quasi_relators = box.qr.word_classes();
  for (vector<string>::iterator it = quasi_relators.begin(); it != quasi_relators.end(); ++it) {
    // So not idenity and absUB(w.b) < 1
    SL2<AJ> w = construct_word(*it, p, short_words_cache); 
    if ((must_fix_x_axis(w,p) && (cant_fix_x_axis(w,p) || non_cylic_power(w, box.x_cover()))) ||
        (must_fix_y_axis(w,p) && (cant_fix_y_axis(w,p) || non_cylic_power(w, box.y_cover()))))
    {
      fprintf(stderr, "killed by failed quasirelator %s at %s\n", (*it).c_str(), box.name.c_str());
      t.aux_word.assign(*it);
      t.test_result = killed_failed_qr;
      return true;
    }
  }

  if (g_options.improve_tree || !t.l_child) {
    for (int i = 0; i < g_tests.size(); ++i) {
      vector<box_state>& th = history[i];
      while (th.size() <= depth && (th.size() < depth-6 || th.empty() || th.back() == open)) {
//        fprintf(stderr, "********************************* Center Test *********************************\n");
//        fprintf(stderr, "%s", box.desc().c_str());
        box_state result = g_tests.evaluate_center(i, place[th.size()]);
//        fprintf(stderr, "********************************* End Center Test *********************************\n");
        th.push_back(result);
      }
      if (th.back() != open) {
        new_qrs.clear();
//        fprintf(stderr, "********************************* Evaluate *********************************\n");
//        fprintf(stderr, "%s", box.desc().c_str());
        box_state result = g_tests.evaluate_box(i, box, aux_word, new_qrs, short_words_cache);
//        fprintf(stderr, "********************************* End Evaluate *********************************\n");

        switch (result) {
          case killed_bounds :
          case killed_only_elliptic : //TODO 
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
          case var_x_hits_y :
          case var_y_hits_x : {
            t.test_index = i;
            t.test_result = result;
            return true;
          }
          case open_with_qr : {
//          fprintf(stderr,"Result %d, new qrs len %d\n", result, new_qrs.size());
            for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
//            fprintf(stderr, "New QR is %s\n", (*it).c_str());
              box.qr.get_name(*it); // Also adds qr to the box's list
            }
            t.qr_desc = box.qr.min_pow_desc();
            break;
          }
          default : {
            continue;
          }
        }
      }
    }
  }

  if (g_options.word_search_depth > 0 && depth > 0 && (g_options.improve_tree || !t.l_child) && box.name.length() > 36 && depth % g_options.word_search_depth == 0) {
    // while (depth - searched_depth > g_options.word_search_depth) {
      //Box& search_place = place[++searched_depth];
      Box& search_place = box;
      // vector<word_pair> search_pairs = find_pairs(search_place.center(), vector<string>(), 1, g_options.max_word_length, box.qr.word_classes());
      vector<word_pair> search_pairs = find_words_v2(search_place.center(), 1, 8, box.qr.word_classes(), map<string, int>());
      //vector<word_pair> search_pairs;
      // fprintf(stderr, "Tube search ran at(%s\n", search_place.name.c_str());
      if (search_pairs.size() > 0) {
        word_pair new_pair = search_pairs.back();

        int old_size = g_tests.size();
        int new_index = g_tests.add(new_pair);
        history.resize(g_tests.size());

        if (old_size < g_tests.size()) {
          fprintf(stderr, "search (%s) found (%s,%s) at (%s)\n",
                  search_place.qr.desc().c_str(), new_pair.first.c_str(), new_pair.second.c_str(), search_place.name.c_str());

          new_qrs.clear();
          box_state result = g_tests.evaluate_box(new_index, box, aux_word, new_qrs, short_words_cache);

          switch (result) {
            case killed_bounds :
            case killed_only_elliptic : //TODO 
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
            case var_x_hits_y :
            case var_y_hits_x : {
              t.test_index = new_index;
              t.test_result = result;
              return true;
            }
            case open_with_qr : {
  //          fprintf(stderr,"Result %d, new qrs len %d\n", result, new_qrs.size());
              for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
  //            fprintf(stderr, "New QR is %s\n", (*it).c_str());
                box.qr.get_name(*it); // Also adds qr to the box's list
              }
              t.qr_desc = box.qr.min_pow_desc();
              break;
            }
            default : {
              // continue;
            }
          }
        }
    //  }
    }
  }

  t.test_index = -1;

  if (!t.l_child) {
    if (depth >= g_options.max_depth || ++g_boxesVisited >= g_options.max_size || ++newDepth > g_options.invent_depth) {
//    fprintf(stderr,"Deph %d, max depth %d, boxes_visited %d, max size %d, newDepth %d, invent depth %d\n", depth, g_options.max_depth, g_boxesVisited, g_options.max_size, newDepth, g_options.invent_depth);
      fprintf(stderr, "HOLE %s (%s)\n", box.name.c_str(), box.qr.desc().c_str());
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
//  printf("%d\n", t.test_result);
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
        printf("%c(%s,)\n", 'Q', t.aux_word.c_str());
        return;
      }
      case killed_only_elliptic : type = 'E'; break; 
      case killed_x_hits_y : type = 'a'; break;
      case killed_y_hits_x : type = 'A'; break;
      case killed_x_tube : type = 'x'; break;
      case killed_y_tube : type = 'y'; break;
      case killed_lox_not_x_power : type = 'p'; break; 
      case killed_lox_not_y_power : type = 'P'; break;
      case killed_move : type = 'm'; break;
      case killed_marg : type = 'M'; break;
      case variety_nbd_x : type = 'v'; break;
      case variety_nbd_y : type = 'V'; break;
      case variety_nbd : type = 'W'; break;
      case var_x_hits_y : type = 'c'; break;
      case var_y_hits_x : type = 'C'; break;
      default : return;
    }
    printf("%c%s\n", type, g_tests.get_name(t.test_index).c_str());
}

