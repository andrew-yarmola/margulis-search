#include <getopt.h>
#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Box.h"
#include "TestCollection.hh"
#include "TubeSearch.hh"
#include "QuasiRelators.h"

using namespace std;

struct Options {
  Options() :
  box_name(""), // Binary representation of box
  words_file("words"), // Previously generated words
  powers_file("null"), // Output from power parabolic.pl
  max_depth(24), // Maximum depth for a file
  truncate_depth(6), 
  invent_depth(12),
  max_size(1000000),
  improve_tree(false),
  word_search_depth(-1),
  fill_holes(false),
  max_word_length(40) {}
  const char* box_name;
  const char* words_file;
  const char* powers_file;
  int max_depth;
  int truncate_depth;
  int invent_depth;
  int max_size;
  bool improve_tree;
  int word_search_depth;
  bool fill_holes;
  int max_word_length;
};

Options g_options;
TestCollection g_tests;
typedef vector< vector< box_state > > TestHistory;
static int g_boxesVisited = 0;

struct PartialTree {
  PartialTree() : l_child(NULL), r_child(NULL), test_index(-1), test_result(open), aux_word(), qr_desc() {}
  PartialTree *l_child;
  PartialTree *r_child;
  int test_index;
  box_state test_result;
  string aux_word;
  string qr_desc;
};

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree read_tree()
{
  PartialTree t;
  char buf[1000];
  if (!fgets(buf, sizeof(buf), stdin)) {
    fprintf(stderr, "unexpected EOF\n");
    exit(1);
  }
  int n = strlen(buf);
  if (buf[n-1] == '\n')
    buf[n-1] = '\0';
  if (buf[0] == 'X') {
    t.l_child = new PartialTree(read_tree());
    t.r_child = new PartialTree(read_tree());
  } else if (strstr(buf, "HOLE") != NULL) {
    t.test_index = -2;
  } else {
    if (isdigit(buf[0])) {
      t.test_index = atoi(buf);
    } else {
      if (strchr("xXyY", buf[2]) != NULL) {
        t.test_index = g_tests.add(string(buf));
      }
    }
  }
  return t;
}

void truncate_tree(PartialTree& t)
{
  if (t.l_child) {
    truncate_tree(*t.l_child);
    delete t.l_child;
    t.l_child = 0;
  }
  if (t.r_child) {
    truncate_tree(*t.r_child);
    delete t.r_child;
    t.r_child = 0;
  }
}

int tree_size(PartialTree& t) {
  int size = 1;
  if (t.l_child)
    size += tree_size(*t.l_child);
  if (t.r_child)
    size += tree_size(*t.r_child);
  return size;
}

//double g_cosh_marg_upper_bound = 1.3175;
//double g_cosh_marg_upper_bound = 1.0454;
double g_cosh_marg_upper_bound = 1.0811;
//double g_cosh_marg_upper_bound = 1.0811;
//double g_cosh_marg_upper_bound = 1.186;
double g_cosh_marg_lower_bound = 1.0054;
// double g_sinh_d_bound = 1.474; 
// double g_sinh_d_bound = 5.0; 
// double g_sinh_d_bound = 0.434; 
double g_sinh_d_bound = 0.63; 

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
      vector<word_pair> search_pairs = find_words_v2(search_place.center(), 1, 6, box.qr.word_classes(), map<string, int>());
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

const char* g_programName;

static struct option long_options[] = {
  {"box",  required_argument, NULL, 'b' },
  {"words", required_argument, NULL, 'w' },
  {"powers", required_argument, NULL, 'p'},
  {"max_depth", required_argument, NULL, 'd' },
  {"invent_depth", required_argument, NULL, 'i' },
  {"improve_tree", no_argument, NULL, 'I'},
  {"truncate_depth", required_argument, NULL, 't' },
  {"max_size", required_argument, NULL, 's' },
  {"word_search_depth", required_argument, NULL, 'B'},
  {"fill_holes", no_argument, NULL, 'f'},
  {NULL, 0, NULL, 0}
};

static char opt_str[1000] = "";
void set_opt_str() {
  char* osp = opt_str;
  for (int i = 0; long_options[i].name; ++i) {
    *osp++ = long_options[i].val;
    if (long_options[i].has_arg != no_argument)
      *osp++ = ':';
  }
  *osp = '\0';
}

void usage()
{
  const char* longUsage = "\
       --box <box_id>\n\
           Box ID for the root of the input and output trees.\n\
    May include characters not in [01], which are ignored.\n\
\n\
Options controlling which relators to use:\n\
  [ --words  <words_file> ]\n\
    File containing starting list of words to try.\n\
  [ --word_search_depth <n> ]\n\
    Perform a search for relators when visiting a node at least n levels deep.\n\
\n\
Options controlling which relators eliminate boxes:\n\
  [ --powers <powers_file> ]\n\
    File containing impossible relator definitions.\n\
    See ImpossibleRelators::load(...)\n\
\n\
Options controlling tree manipulation:\n\
  [ --max_depth <n> ]\n\
    Don't descend more than n levels deeper than the root box.\n\
  [ --invent_depth <n> ]\n\
    Don't descend more than n levels deeper than the terminal node of the input tree.\n\
  [ --max_size <n> ]\n\
    Don't allow the output tree to have more than n nodes.\n\
  [ --truncate_depth <n> ]\n\
    Don't emit holes more than n levels deeper than the root node.\n\
    Instead, replace the subtree-with-holes with a single hole.\n\
  [ --improve_tree ]\n\
    If set, attempt to directly eliminate internal nodes of the input tree.\n\
  [ --fill_holes ]\n\
    If set, attempt to patch holes in the input tree.\n\
";
  fprintf(stderr, "Usage: %s %s\n\n%s", g_programName, opt_str, longUsage);
}

void load_words(set<string>& s, const char* fileName)
{
  FILE* fp = fopen(fileName, "r");
  char buf[10000];
  while (fp && fgets(buf, sizeof(buf), fp)) {
    int n = -1 + strlen(buf);
    if (buf[n] == '\n')
      buf[n] = '\0';
    s.insert(buf);
  }
}

int main(int argc, char** argv)
{
  set_opt_str();
  if (argc < 2) {
    usage();
    exit(1);
  }

  int ch;
  while ((ch = getopt_long(argc, argv, opt_str, long_options, NULL)) != -1) {
    fprintf(stderr,"Arg %c, %s\n", ch, optarg);
    switch(ch) {
    case 'b': g_options.box_name = optarg; break;
    case 'w': g_options.words_file = optarg; break;
    case 'p': g_options.powers_file = optarg; break;
    case 'd': g_options.max_depth = atoi(optarg); break;
    case 'i': g_options.invent_depth = atoi(optarg); break;
    case 'I': g_options.improve_tree = true; break;
    case 't': g_options.truncate_depth = atoi(optarg); break;
    case 's': g_options.max_size = atoi(optarg); break;
    case 'B': g_options.word_search_depth = atoi(optarg); break;
    case 'f': g_options.fill_holes = true; break;
    }
  }

  Box box;
  for (const char* boxP = g_options.box_name; *boxP; ++boxP) {
    if (*boxP == '0') {
      box = box.child(0);
    } else if (*boxP == '1') {
      box = box.child(1);
    }
  }
 
  g_tests.load(g_options.words_file);
  g_tests.load_impossible_relations(g_options.powers_file);

  fprintf(stderr, "%s", box.desc().c_str());
  PartialTree t = read_tree();
  refine_tree(box, t);
  print_tree(t);
  fprintf(stderr, "%d nodes added\n", g_boxesVisited);
}
