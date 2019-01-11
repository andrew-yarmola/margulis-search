/*
 *  TestCollection.cpp
 *  mom
 *
 *  Created by Nathaniel Thurston on 27/09/2007.
 *  Copyright 2007 THingith ehf.. All rights reserved.
 *
 */
// The bounds:
// 0. |sl| >= 1 (horoball size)
// 1. Im(sl) >= 0 (only sl^2 matters)
// 2. -1/2 <= Re(l) <= 1/2 (reduction modulo M)
// 3. Im(l >= 0 (negation)
// 4. 0 <= Im(p) <= Im(l)/2 (reduction modulo N, flipping sign of N)
// 5. 0 <= Re(p) <= 1/2 (reduction modulo M, flipping sign of M)
// 6. |sl^2| Im(n) <= 4 (area of fundamental paralleogram)
// 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TestCollection.h"
#include "ImpossibleRelations.h"
#include <algorithm>
using namespace std;
// using namespace __gnu_cxx;

int TestCollection::size()
{
	return 7 + indexString.size();
}

box_state TestCollection::evaluate_approx(word_pair pair, const Box& box)
{
    SL2<Complex> w1 = construct_word(word.first, box.center());
    SL2<Complex> w2 = construct_word(word.second, box.center());
    if (inside_var_nbd(w1, box.center())) return variety_center;
    if (tubes_intersect(w1, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    if (inside_var_nbd(w2, box.center())) return variety_center;
    if (tubes_intersect(w2, box.x_center(), box.y_center(), box.center())) return bad_tubes_center;
    if (margulis_larger_than_cutoff(w1, w2, g_margulis_bound)) return bad_margulis_center;
    return open;
}

box_state TestCollection::evaluate_ACJ(string word, Params<ACJ>& params, string& aux_word, vector<string>& new_qrs,
                                       unordered_map<int,ACJ>& para_cache, unordered_map<string,SL2ACJ>& words_cache)
{
    box_state state = open;
    bool found_qrs = false;
    aux_word.assign(word);
    int g_len = g_length(word);
    double one = 1; // Exact
	  SL2ACJ w = construct_word(word, params, para_cache, words_cache);

    if (g_len <= g_max_g_len && inside_var_nbd(w)) return variety_nbd;

	if (large_horoball(w,params)) {

        if (not_para_fix_inf(w)) {

			return killed_no_parabolics;

		} else {

			vector<string> mandatory;
			bool isImpossible = impossible->isAlwaysImpossible(word, mandatory);

			if (isImpossible) return killed_parabolics_impossible;
            else if (mandatory.size() > 0) {
                for (vector<string>::iterator it = mandatory.begin(); it != mandatory.end(); ++it) {
                    SL2ACJ w_sub = construct_word(*it, params, para_cache, words_cache);
                    if (not_para_fix_inf(w_sub)) {
                        aux_word.assign(*it);
                        return killed_elliptic;
                    }
                }
            }
   
            // Look for lattice points. We guess at the center
            // No reason to look for a unique lattice point if w.b has large size
            ACJ L = params.lattice;
            if (absLB(L) > 2*w.b.size) {
                XComplex cL = L.f;
                XComplex cT = (absUB((w.d.f - one).z) < 2 || absUB((w.a.f - one).z) < 2) ? w.b.f : -w.b.f;
                // XComplex cT = (w.b.f/w.d.f).z;
                // We expect T to be near the lattice point M_pow + N_pow*L
                int N_pow = (int) floor(cT.im / cL.im);
                int M_pow = (int) floor((cT - (cL * N_pow).z).z.re);
                // We look over 16 nearby lattice points
                int s[4] = {0,-1,1,2};
                SL2ACJ w_k;
                ACJ T;
                pair<unordered_map<int,ACJ>::iterator,bool> lookup_para;
                int N, M;
                for (int i = 0; i < 4; ++i) {
                    N = N_pow + s[i];
                    for (int j = 0; j < 4; ++j) {
                        state = open;
                        M = M_pow + s[j];
                        if (i == 0 && j == 0) {
                            w_k = w;
                        } else {
                            if (abs(M) > 1024 || abs(N) > 1024) { fprintf(stderr, "Error constructing word: huge translation\n"); }
                            lookup_para = para_cache.emplace(4096*M+N, ACJ());
                            if (lookup_para.second) {
                                T = params.lattice*double(N) + double(M);
                                swap(lookup_para.first->second, T);
                            }
                            // Shift to "0"
                            w_k = SL2ACJ(w.a - lookup_para.first->second * w.c, w.b - lookup_para.first->second * w.d, w.c, w.d); // Cheaper multiplying
                            // What if we now have a variety word?
                            if (g_len <= g_max_g_len && inside_var_nbd(w_k)) { // TODO: Test with constucted word!
                                state = variety_nbd;
                                break;
                            }
                            if (not_para_fix_inf(w_k)) {
                                state = killed_no_parabolics;
                                break;
                            }
                            if (absUB(w_k.b) < 1) {
                                isImpossible = impossible->isImpossible(word, M, N, mandatory);
                                if (isImpossible) {
                                    state = killed_identity_impossible;
                                    break;
                                }
                                // Mandaotry includes list of things that must be parabolic. If they are not parabolic
                                // anywhere in the box, we can kill the box
                                for (vector<string>::iterator it = mandatory.begin(); it != mandatory.end(); ++it) {
                                    SL2ACJ w_sub = construct_word(*it, params, para_cache, words_cache);
                                    if (not_para_fix_inf(w_sub)) {
                                        aux_word.assign(*it);
                                        return killed_elliptic;
                                    }
                                }
                            }
                        }
                        if (absUB(w_k.b) < 1 && absLB(w_k.b) > 0) {
    //                        string word_k = shifted_word(word, - M, - N);
    //                        fprintf(stderr, "Killed by Failed qr %s\n", word_k.c_str());
    //                        SL2ACJ new_w_k = construct_word(word_k, params);
    //                        SL2ACJ gah_k = constructT(params, - M, - N) * w;
    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
    //                                        absLB(w_k.b), absLB(w_k.c), absLB(w_k.a - 1.), absLB(w_k.d - 1.), absLB(w_k.a + 1.), absLB(w_k.d + 1.));
    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
    //                                        absLB(new_w_k.b), absLB(new_w_k.c), absLB(new_w_k.a - 1.), absLB(new_w_k.d - 1.), absLB(new_w_k.a + 1.), absLB(new_w_k.d + 1.));
    //                        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
    //                                        absLB(gah_k.b), absLB(gah_k.c), absLB(gah_k.a - 1.), absLB(gah_k.d - 1.), absLB(gah_k.a + 1.), absLB(gah_k.d + 1.));
                            state = killed_failed_qr;
                            break;
                        }
                        if (only_bad_parabolics(w_k, params)) { // TODO: Test with constucted word!
                            // w_k is a bad parabolic
                            state = killed_bad_parabolic;
                            break;
                        }
                        if (absUB(w_k.b) < 1) {
                            // If nothing has worked, at least add w as a quasi relator
                            state = open_with_qr;
                            found_qrs = true;
                            new_qrs.push_back(shifted_word(word, - M, - N));
    //                        string word_k = shifted_word(word, - M, - N);
    //                        SL2ACJ new_w_k = construct_word(word_k, params);
    //                        fprintf(stderr,"Horo Ratio for new QR is %f\n", absUB(w_k.c / params.loxodromic_sqrt));
    //                        fprintf(stderr,"Reconstucted horo ratio for QR is %f\n", absUB(new_w_k.c / params.loxodromic_sqrt));
                        }
                    }
                }
                if (state != open && state != open_with_qr) {
                    aux_word.assign(shifted_word(word, - M, - N));
                    return state;
                }
            }
        }
    }
    if (found_qrs) return open_with_qr; 
    return open;
}

box_state check_bounds_center(bool result) {
    if (result) return out_of_bounds_center;
    else return open;
}

double g_latticeArea;
box_state TestCollection::evaluateCenter(int index, Box& box)
{
	Params<XComplex> params = box.center();
	switch(index) {
		case 0:	{
			XComplex sl = params.loxodromic_sqrt;
			return check_bounds_center(sl.re*sl.re + sl.im*sl.im < 1.0);
		}
		case 1: return check_bounds_center( params.loxodromic_sqrt.im < 0.0
                                         || params.lattice.im < 0.0
                                         || params.parabolic.im < 0.0
                                         || params.parabolic.re < 0.0);
		case 2: return check_bounds_center(fabs(params.lattice.re) > 0.5);
		case 3: return check_bounds_center(absUB(params.lattice) < 1);
		case 4: return check_bounds_center(params.parabolic.im > 0.5*params.lattice.im);
		case 5: return check_bounds_center(params.parabolic.re > 0.5);
		case 6: {
			g_latticeArea = pow(absLB(params.loxodromic_sqrt),2)*params.lattice.im;
			return check_bounds_center(g_latticeArea > g_maximumArea);
		}
		default:
			return evaluate_approx(indexString[index-7], params);
	}
}

box_state check_bounds(bool result) {
    if (result) return killed_bounds;
    else return open;
}

box_state TestCollection::evaluateBox(int index, Box& box, string& aux_word, vector<string>& new_qrs,
                                      unordered_map<int,ACJ>& para_cache, unordered_map<string,SL2ACJ>& words_cache)
{
	Params<XComplex> nearer = box.nearer();
	Params<XComplex> further = box.further();
	Params<XComplex> greater = box.greater();
	switch(index) {
		case 0: return check_bounds(absUB(further.loxodromic_sqrt) < 1.0);
		case 1: return check_bounds(greater.loxodromic_sqrt.im < 0.0
                                 || greater.lattice.im < 0.0
                                 || greater.parabolic.im < 0.0
                                 || greater.parabolic.re < 0.0);
		case 2: return check_bounds(fabs(nearer.lattice.re) > 0.5);
		case 3: return check_bounds(absUB(further.lattice) < 1.0);
        // Note: we can exclude the box if and only if the parabolic imag part is
        // bigger than half the lattice imag part over the WHOLE box
        // We assume that case 1 has been tested. Multiplication by 0.5 is EXACT (if no underflow or overflow)
		case 4: return check_bounds(nearer.parabolic.im > 0.5*further.lattice.im);
		case 5: return check_bounds(nearer.parabolic.re > 0.5);
		case 6: {
            // Area is |lox_sqrt|^2*|Im(lattice)|.
            double area = areaLB(nearer);
		    return check_bounds(area > g_maximumArea);
		}
		default: {
			Params<ACJ> cover(box.cover());
			box_state result = evaluate_ACJ(indexString[index-7], cover, aux_word, new_qrs, para_cache, words_cache);
//			if (result) {
//                // TODO: Understand this tail enumeration that adds words based on given word
//				enumerate(indexString[index-7].c_str());
//			}
			return result;
		}
	}
}

// Returns the index number for the first basic 7 tests
// or the quasi-relator if the index is 7 or above
const char* TestCollection::getName(int index)
{
	static char buf[4];
	if (index < 7) {
		sprintf(buf, "%d", index);
		return buf;
	} else {
		return indexString[index-7].c_str();
	}
}


int TestCollection::add(string buf)
{
    size_t start = buf.find('(');   
    string word;
    if (start != string::npos) {
        size_t end = buf.find(')');
        word = buf.substr(start+1,end-start-1);
    } else {
        word = buf;
    }  
	map<string, int>::iterator it = stringIndex.find(word);
	if (it == stringIndex.end()) {
//		fprintf(stderr, "adding %lu=%s\n", indexString.size(), word.c_str());
		stringIndex[word] = indexString.size();
		indexString.push_back(word);
		return indexString.size()+6;
	}
	else
		return it->second+7;
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
