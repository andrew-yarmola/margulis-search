#ifndef __refine_h
#define __refine_h
#include <getopt.h>
#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Box.h"

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

struct PartialTree {
  PartialTree() : l_child(NULL), r_child(NULL), test_index(-1), test_result(open), aux_word(), qr_desc() {}
  PartialTree *l_child;
  PartialTree *r_child;
  int test_index;
  box_state test_result;
  std::string aux_word;
  std::string qr_desc;
};

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree read_tree();

void refine_tree(Box box, PartialTree& t);

void print_tree(PartialTree& t);

#endif // __refine_h
