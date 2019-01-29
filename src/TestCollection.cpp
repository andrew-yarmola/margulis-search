#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TestCollection.hh"
#include "ImpossibleRelations.h"
#include <algorithm>
using namespace std;
// using namespace __gnu_cxx;

extern double g_4_cosh_margulis_bound;
extern double g_exp_half_margulis_bound; 

int num_bound_tests = 3;

int TestCollection::size()
{
	return num_bound_tests + pairVector.size();
}

box_state TestCollection::evaluate_approx(word_pair pair, const Box& box)
{
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n     Word Pair: %s and %s\n +++++++++++++++++++++++++++++++++++++++++++\n", pair.first.c_str(), pair.second.c_str());
    // we can assume that words pairs use canonical names
    SL2<Complex> w1 = construct_word(pair.first, box.center());
    SL2<Complex> w2 = construct_word(pair.second, box.center());

    if (y_power(pair.first) > 0) {
      if (inside_var_nbd_x(w1, box.center())) return variety_center;
      if (tubes_intersect_x(w1, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    }
    if (y_power(pair.second) > 0) {
      if (inside_var_nbd_x(w2, box.center())) return variety_center;
      if (tubes_intersect_x(w2, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    }
    if (x_power(pair.first) > 0) {
      if (inside_var_nbd_y(w2, box.center())) return variety_center;
      if (tubes_intersect_y(w1, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    }
    if (x_power(pair.second) > 0) {
      if (inside_var_nbd_y(w2, box.center())) return variety_center;
      if (tubes_intersect_y(w2, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    }

    if (margulis_smaller_than_xy(w1, w2, box.x_center(), box.y_center())) return bad_margulis_center;
    return open;
}

box_state TestCollection::evaluate_ACJ(word_pair pair, const Box& box, string& aux_word,
                                       vector<string>& new_qrs, unordered_map< string,SL2<ACJ> >& words_cache)
{
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n     Word Pair: %s and %s\n +++++++++++++++++++++++++++++++++++++++++++\n", pair.first.c_str(), pair.second.c_str());
    // we can assume that words pairs use canonical names
    SL2<ACJ> w1 = construct_word(pair.first, box.cover());
    SL2<ACJ> w2 = construct_word(pair.second, box.cover());

    if (y_power(pair.first) > 0) {
      if (inside_var_nbd_x(w1, box.cover())) return variety_nbd;
      if (tubes_intersect_x(w1, box.x_cover(), box.y_cover(), box.cover())) return killed_bad_tubes;
    }
    if (y_power(pair.second) > 0) {
      if (inside_var_nbd_x(w2, box.cover())) return variety_nbd;
      if (tubes_intersect_x(w2, box.x_cover(), box.y_cover(), box.cover())) return killed_bad_tubes;
    }
    if (x_power(pair.first) > 0) {
      if (inside_var_nbd_y(w2, box.cover())) return variety_nbd;
      if (tubes_intersect_y(w1, box.x_cover(), box.y_cover(), box.cover())) return killed_bad_tubes;
    }
    if (x_power(pair.second) > 0) {
      if (inside_var_nbd_y(w2, box.cover())) return variety_nbd;
      if (tubes_intersect_y(w2, box.x_cover(), box.y_cover(), box.cover())) return killed_bad_tubes;
    }

    if (margulis_smaller_than_xy(w1, w2, box.x_cover(), box.y_cover())) return killed_margulis;
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

box_state TestCollection::evaluateCenter(int index, Box& box)
{
  fprintf(stderr, "Evaluating center test index %d\n", index);
	Params<Complex> center = box.center();
	switch(index) {
		case 0:	{ // re_length of both geodesics is short enough
      return check_bounds_center(g_exp_half_margulis_bound <= absLB(center.coshL2 + center.sinhL2));
    } 
		case 1: {
      return check_bounds_center(g_exp_half_margulis_bound <= absLB(center.coshD2 + center.sinhD2));
    }
		case 2: {
      return check_bounds_center(margulis_larger_than_cutoff(box.x_center(), box.y_center(), g_4_cosh_margulis_bound));
    }
		default:
			return evaluate_approx(pairVector[index - num_bound_tests], box);
	}
}

box_state TestCollection::evaluateBox(int index, Box& box, string& aux_word, vector<string>& new_qrs, unordered_map< string,SL2<ACJ> >& words_cache)
{
  fprintf(stderr, "Evaluating box test index %d\n", index);
	Params<ACJ> cover = box.cover();
	switch(index) {
		case 0:	{ // re_length of both geodesics is short enough
      return check_bounds(g_exp_half_margulis_bound <= absLB(cover.coshL2 + cover.sinhL2)); 
    } 
		case 1: {
      return check_bounds(g_exp_half_margulis_bound <= absLB(cover.coshD2 + cover.sinhD2));
    }
		case 2: {
      return check_bounds(margulis_larger_than_cutoff(box.x_cover(), box.y_cover(), g_4_cosh_margulis_bound));
    }
		default:
			return evaluate_ACJ(pairVector[index - num_bound_tests], box, aux_word, new_qrs, words_cache);
	}
}

// Returns the index number for the first basic 2 tests
// or the quasi-relator if the index is 2 or above
const char* TestCollection::getName(int index)
{
	static char buf[500];
	if (index < num_bound_tests) {
		sprintf(buf, "%d", index);
	} else {
    word_pair p = pairVector[index - num_bound_tests];
		sprintf(buf, "(%s,%s)", p.first.c_str(), p.second.c_str());
	}
  return buf;
}

int TestCollection::add(string buf)
{
  size_t start = buf.find('(');   
  size_t comma = buf.find(',');   
  string first;
  string second;
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

int TestCollection::add(const word_pair& p) { 
  map< word_pair,int >::iterator it = pairIndex.find(p);
  if (it == pairIndex.end()) {
//  fprintf(stderr, "adding %lu=%s\n", pairVector.size(), word.c_str());
    pairIndex[p] = pairVector.size();
    pairVector.push_back(p);
    return pairVector.size() + num_bound_tests - 1;
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

void TestCollection::loadImpossibleRelations(const char* fileName)
{
	impossible = ImpossibleRelations::create(fileName);
}

//box_state TestCollection::evaluate_ACJ(string word, Params<ACJ>& params, string& aux_word, vector<string>& new_qrs,
//                                       unordered_map<int,ACJ>& para_cache, unordered_map<string,SL2ACJ>& words_cache)
//{
//    box_state state = open;
//    bool found_qrs = false;
//    aux_word.assign(word);
//    int g_len = g_length(word);
//    double one = 1; // Exact
//	  SL2ACJ w = construct_word(word, params, para_cache, words_cache);
//
//    if (g_len <= g_max_g_len && inside_var_nbd(w)) return variety_nbd;
//
//	if (large_horoball(w,params)) {
//
//        if (not_para_fix_inf(w)) {
//
//			return killed_no_parabolics;
//
//		} else {
//
//			vector<string> mandatory;
//			bool isImpossible = impossible->isAlwaysImpossible(word, mandatory);
//
//			if (isImpossible) return killed_parabolics_impossible;
//            else if (mandatory.size() > 0) {
//                for (vector<string>::iterator it = mandatory.begin(); it != mandatory.end(); ++it) {
//                    SL2ACJ w_sub = construct_word(*it, params, para_cache, words_cache);
//                    if (not_para_fix_inf(w_sub)) {
//                        aux_word.assign(*it);
//                        return killed_elliptic;
//                    }
//                }
//            }
//   
//            // Look for lattice points. We guess at the center
//            // No reason to look for a unique lattice point if w.b has large size
//            ACJ L = params.lattice;
//            if (absLB(L) > 2*w.b.size) {
//                XComplex cL = L.f;
//                XComplex cT = (absUB((w.d.f - one).z) < 2 || absUB((w.a.f - one).z) < 2) ? w.b.f : -w.b.f;
//                // XComplex cT = (w.b.f/w.d.f).z;
//                // We expect T to be near the lattice point M_pow + N_pow*L
//                int N_pow = (int) floor(cT.im / cL.im);
//                int M_pow = (int) floor((cT - (cL * N_pow).z).z.re);
//                // We look over 16 nearby lattice points
//                int s[4] = {0,-1,1,2};
//                SL2ACJ w_k;
//                ACJ T;
//                pair<unordered_map<int,ACJ>::iterator,bool> lookup_para;
//                int N, M;
//                for (int i = 0; i < 4; ++i) {
//                    N = N_pow + s[i];
//                    for (int j = 0; j < 4; ++j) {
//                        state = open;
//                        M = M_pow + s[j];
//                        if (i == 0 && j == 0) {
//                            w_k = w;
//                        } else {
//                            if (abs(M) > 1024 || abs(N) > 1024) { fprintf(stderr, "Error constructing word: huge translation\n"); }
//                            lookup_para = para_cache.emplace(4096*M+N, ACJ());
//                            if (lookup_para.second) {
//                                T = params.lattice*double(N) + double(M);
//                                swap(lookup_para.first->second, T);
//                            }
//                            // Shift to "0"
//                            w_k = SL2ACJ(w.a - lookup_para.first->second * w.c, w.b - lookup_para.first->second * w.d, w.c, w.d); // Cheaper multiplying
//                            // What if we now have a variety word?
//                            if (g_len <= g_max_g_len && inside_var_nbd(w_k)) { // TODO: Test with constucted word!
//                                state = variety_nbd;
//                                break;
//                            }
//                            if (not_para_fix_inf(w_k)) {
//                                state = killed_no_parabolics;
//                                break;
//                            }
//                            if (absUB(w_k.b) < 1) {
//                                isImpossible = impossible->isImpossible(word, M, N, mandatory);
//                                if (isImpossible) {
//                                    state = killed_identity_impossible;
//                                    break;
//                                }
//                                // Mandaotry includes list of things that must be parabolic. If they are not parabolic
//                                // anywhere in the box, we can kill the box
//                                for (vector<string>::iterator it = mandatory.begin(); it != mandatory.end(); ++it) {
//                                    SL2ACJ w_sub = construct_word(*it, params, para_cache, words_cache);
//                                    if (not_para_fix_inf(w_sub)) {
//                                        aux_word.assign(*it);
//                                        return killed_elliptic;
//                                    }
//                                }
//                            }
//                        }
//                        if (absUB(w_k.b) < 1 && absLB(w_k.b) > 0) {
//    //                        string word_k = shifted_word(word, - M, - N);
//    //                        fprintf(stderr, "Killed by Failed qr %s\n", word_k.c_str());
//    //                        SL2ACJ new_w_k = construct_word(word_k, params);
//    //                        SL2ACJ gah_k = constructT(params, - M, - N) * w;
//    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
//    //                                        absLB(w_k.b), absLB(w_k.c), absLB(w_k.a - 1.), absLB(w_k.d - 1.), absLB(w_k.a + 1.), absLB(w_k.d + 1.));
//    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
//    //                                        absLB(new_w_k.b), absLB(new_w_k.c), absLB(new_w_k.a - 1.), absLB(new_w_k.d - 1.), absLB(new_w_k.a + 1.), absLB(new_w_k.d + 1.));
//    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
//    //                                        absLB(gah_k.b), absLB(gah_k.c), absLB(gah_k.a - 1.), absLB(gah_k.d - 1.), absLB(gah_k.a + 1.), absLB(gah_k.d + 1.));
//                            state = killed_failed_qr;
//                            break;
//                        }
//                        if (only_bad_parabolics(w_k, params)) { // TODO: Test with constucted word!
//                            // w_k is a bad parabolic
//                            state = killed_bad_parabolic;
//                            break;
//                        }
//                        if (absUB(w_k.b) < 1) {
//                            // If nothing has worked, at least add w as a quasi relator
//                            state = open_with_qr;
//                            found_qrs = true;
//                            new_qrs.push_back(shifted_word(word, - M, - N));
//    //                        string word_k = shifted_word(word, - M, - N);
//    //                        SL2ACJ new_w_k = construct_word(word_k, params);
//    //                        fprintf(stderr,"Horo Ratio for new QR is %f\n", absUB(w_k.c / params.loxodromic_sqrt));
//    //                        fprintf(stderr,"Reconstucted horo ratio for QR is %f\n", absUB(new_w_k.c / params.loxodromic_sqrt));
//                        }
//                    }
//                }
//                if (state != open && state != open_with_qr) {
//                    aux_word.assign(shifted_word(word, - M, - N));
//                    return state;
//                }
//            }
//        }
//    }
//    if (found_qrs) return open_with_qr; 
//    return open;
//}


//int g_maxWordLength = 2;
//void TestCollection::enumerate(const char* w)
//{
//	static vector<int> maxP;
//	if (!*w && !maxP.empty())
//		return;
//	int pCount=0, lCount=0;
//	while (*w) {
//		if (*w == 'g' || *w == 'G')
//			++lCount;
//		else
//			++pCount;
//		++w;
//	}
//	if (lCount == 0) lCount = 1;
//	if (lCount >= maxP.size()) {
//		maxP.resize(lCount+1, -1);
//	}
//	if (pCount <= maxP[lCount]) {
//		return;
//	}
//	
//	if (lCount + pCount > g_maxWordLength)
//		g_maxWordLength = lCount + pCount;
//	
//	maxP[lCount] = pCount;
//	
//	if (pCount > 2)
//		pCount = 2;
//	if (lCount > 2)
//		lCount = 2;
////	printf("ENUMERATING %d,%d\n", pCount, lCount);
////	enumerateTails("", pCount, lCount);
//}
//
//void TestCollection::enumerateTails(string s, int pCount, int lCount)
//{
//	if (pCount < -1 || lCount < -1 || (pCount == -1 && lCount == -1))
//		return;
//	const char* p = "Gg";
//	if (s.size() > 0) {
//		char last = s[s.size()-1];
//		switch(last) {
//			case 'G': p = "GMmNn"; break;
//			case 'g': p = "gMmNn"; break;
//			case 'M': p = "GgMNn"; break;
//			case 'm': p = "GgmNn"; break;
//			case 'N': p = "GgN"; break;
//			case 'n': p = "Ggn"; break;
//		}
//	}
//	for (; *p; ++p) {
//		string n = s;
//		n.append(1, *p);
//		if (*p == 'g' || *p == 'G') {
//			add(n);
//			enumerateTails(n, pCount, lCount-1);
//		} else {
//			enumerateTails(n, pCount-1, lCount);
//		}
//	}
//}
//
//string checkPower(string word, int x, int y)
//{
//	char buf[200];
//	char *bp = buf;
//	if (abs(x) > 10 || abs(y) > 10)
//		return "";
//	while (x > 0) { *bp++ = 'm'; --x; }
//	while (x < 0) { *bp++ = 'M'; ++x; }
//	while (y > 0) { *bp++ = 'n'; --y; }
//	while (y < 0) { *bp++ = 'N'; ++y; }
//	strcpy(bp, word.c_str());
//	int l = strlen(buf);
//	for (int n = 2; n+n <= l; ++n) {
//		int k;
//		for ( k = 1; k*n < l; ++k) ;
//		if (k*n == l) {
//			for (--k; k > 0; --k) {
//				if (strncmp(buf, buf+k*n, n))
//					break;
//			}
//			if (k == 0) {
////				fprintf(stderr, "id %s x%d\n", buf, n);
//				g_testCollectionFullWord = buf;
//				buf[n] = '\0';
//				return buf;
//			}
//		}
//	}
////	fprintf(stderr, "fullWord = %s\n", buf);
//	return "";
//}
