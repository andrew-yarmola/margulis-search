#ifndef __partialtree_h
#define __partialtree_h

#include "types.hh"

struct PartialTree {
  PartialTree() : l_child(NULL), r_child(NULL), test_index(-1), test_result(open), aux_word(), qr_desc() {}
  PartialTree *l_child;
  PartialTree *r_child;
  int test_index;
  box_state test_result;
  box_state aux_result;
  std::string aux_word;
  std::string qr_desc;
};

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree read_tree();

void truncate_tree(PartialTree& t);

int tree_size(PartialTree& t);

#endif // __partialtree_h
