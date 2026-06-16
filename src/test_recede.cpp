// Standalone validation for word_provably_receding (subtree word-pruning predicate).
// Checks, on a given box:
//   (1) the predicate runs and is non-vacuous (true for some words, false for others);
//   (2) the soundness invariant: recede==true implies NO kill condition currently fires
//       (a kill needs diff>0 provably, recede needs diff<0 provably — they cannot coexist);
//   (3) the subtree property: if recede holds on a box, it holds on both children
//       (children are sub-boxes, so a provable bound on the parent must persist).
#include <stdio.h>
#include <string>
#include "Box.h"
#include "types.hh"
#include "SL2.hh"
#include "AJCC.h"
#include "IsomH3.hh"
#include "Params.hh"
#include "Generators.hh"
#include "TestCollection.hh"

bool g_debug = false;
double g_cosh_marg_upper_bound = 1.2947;
double g_cosh_marg_lower_bound = 1.0054;
double g_cosh_r_bound = 1e11;
bool g_symmetric = true;

using namespace std;

static int check_box(const string& pos, bool recurse, int& violations) {
  Box box;
  for (char d : pos) box = (d == '1') ? box.child(1) : box.child(0);
  Params<AJCC> p = box.cover();
  const char* words[] = {"x","y","X","Y","xy","xY","Xy","xxy","xyxy","xxxxy",
                         "xyXYxy","xyyyyX","XYxyXY","xYxYxY","xxxxxxyy","xyxyxyxy"};
  int n = sizeof(words) / sizeof(words[0]);
  int recede_true = 0;
  printf("box '%s'\n  %-10s %4s %4s %4s %4s | %s\n",
         pos.c_str(), "word", "move", "wyht", "mvY", "vnbY", "recede");
  for (int i = 0; i < n; ++i) {
    string wd = words[i];
    SL2<AJCC> w = construct_word(wd, p);
    bool mv = move_less_than_marg(w, p);
    bool wh = wy_hits_sym_axis(w, p);
    bool my = (x_power(wd) > 0) && moves_y_axis_too_close_to_x(w, p);
    bool vn = (x_power(wd) > 0) && inside_var_nbd_y(w, p);
    bool rec = word_provably_receding(w, wd, p);
    printf("  %-10s %4d %4d %4d %4d | %d\n", wd.c_str(), mv, wh, my, vn, rec);
    if (rec) ++recede_true;
    if (rec && (mv || wh || my || vn)) {
      printf("    *** VIOLATION: recede==true but a kill condition fires!\n");
      ++violations;
    }
    // For non-killing words, show which receding sub-conditions hold (1 = provably
    // receding for that condition). recede == AND of all applicable sub-conditions.
    if (!(mv || wh || my || vn)) {
      int c1 = strictly_pos(cosh_move_j(w) - p.coshmu);
      AJCC t0 = two_sinh_perp2_sq_way_zero_inf(w, p);
      int c2a = strictly_pos(two_cosh_dist(t0) - p.twocoshreD2);
      AJCC t1 = four_sinh_perp2_sq_way_mp_one(w, p);
      int c2b = strictly_pos(four_cosh_dist(t1) - p.twocoshreD2 * 2);
      AJCC t2 = four_sinh_perp2_sq_way_mp_iye(w, p);
      int c2c = strictly_pos(four_cosh_dist(t2) - p.twocoshreD2 * 2);
      int c3 = 1, c4a = 1, c4b = 1, c4c = 1;
      if (x_power(wd) > 0) {
        c3 = strictly_pos(four_cosh_dist_ax_way(w, p) - p.coshreD * 4);
        c4a = strictly_pos(jorgensen_wy(w, p) - 1.0);
        c4b = strictly_pos(jorgensen_yw(w, p) - 1.0);
        c4c = strictly_pos(four_cosh_dist_ay_way(w, p) - p.coshreD * 4);
      }
      printf("      sub[move=%d wy0=%d wy1=%d wyi=%d mvY=%d jwy=%d jyw=%d fixY=%d]\n",
             c1, c2a, c2b, c2c, c3, c4a, c4b, c4c);
    }
    // Subtree property: recede on parent must imply recede on both children.
    if (rec && recurse) {
      for (int dir = 0; dir < 2; ++dir) {
        Box c = box.child(dir);
        Params<AJCC> cp = c.cover();
        SL2<AJCC> cw = construct_word(wd, cp);
        if (!word_provably_receding(cw, wd, cp)) {
          printf("    *** VIOLATION: recede on parent but NOT on child %d for %s\n",
                 dir, wd.c_str());
          ++violations;
        }
      }
    }
  }
  printf("  recede=true for %d/%d words\n", recede_true, n);
  return recede_true;
}

int main(int argc, char** argv) {
  int violations = 0, any_recede = 0;
  if (argc > 1) {
    any_recede += check_box(argv[1], true, violations);
  } else {
    // Deep boxes (production runs at depth ~100-180): at this depth the AJCC cover is
    // tight enough that "far" words are provably receding, exercising the predicate.
    const char* seeds[] = {"1010110010", "1101001011", "1011001101", "1100101101",
                           "0110100110", "1001011010"};
    for (auto s : seeds) {
      string pos;
      for (int r = 0; r < 5; ++r) pos += s;  // depth 50
      any_recede += check_box(pos, true, violations);
    }
  }
  printf("\nTOTAL: invariant violations = %d ; words proven receding = %d\n",
         violations, any_recede);
  return violations > 0 ? 1 : 0;
}
