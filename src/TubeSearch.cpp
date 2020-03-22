#include <list>
#include <map>
#include <queue>
#include <algorithm>
#include "TubeSearch.hh"
#include "CanonicalName.hh"
#include "IsomH3.hh"
#include "roundoff.h"

using namespace std;

namespace TubeSearchImpl {
	struct Word {
		string name;
		SL2<Complex> matrix;
	};
	
	Word operator * (const Word& lhs, const Word& rhs); // Multiply/concatiante two words
	
	struct WordPair {
		WordPair(Word& first, Word& second, const Params<Complex>& params);
		Word first;
		Word second;
		double ortho_dist_LB;
		double ortho_dist_UB;
    double jorgensen_ww;
    double jorgensen_xw1 = 100000000000;
    double jorgensen_xw2 = 100000000000;
    double jorgensen_yw1 = 100000000000;
    double jorgensen_yw2 = 100000000000;
    double jorgensen_w1x = 100000000000;
    double jorgensen_w2x = 100000000000;
    double jorgensen_w1y = 100000000000;
    double jorgensen_w2y = 100000000000;
	};
	
	struct PairComp : public binary_function<const WordPair&, const WordPair&, bool>
	{
		bool operator() (const WordPair& a, const WordPair& b)
		{
      if (a.ortho_dist_UB > 0 && b.ortho_dist_UB > 0) {
			  return a.ortho_dist_UB < b.ortho_dist_UB;
      } else {
        return a.jorgensen_ww > b.jorgensen_ww;
      }
		}
	};
	
	struct TubeSearch {
		TubeSearch(Params<Complex> params_) :params(params_),
               x(construct_x(params)), X(inverse(x)),
               y(construct_y(params)), Y(inverse(y)),
               exp_re_perp((params.expdx * params.expdy).real()), m_foundGoodPair(false) {
      x_word.name = "x";
      x_word.matrix = x;
      y_word.name = "y";
      y_word.matrix = y;
    }

		~TubeSearch() {}
		void pushWord(Word word);
		void pushWord(string word);
		WordPair findPair();
		bool foundGoodPair() { return m_foundGoodPair; }
		void findNames(Word& word);
		void addRelator(string word);
		
		Params<Complex> params;
		CanonicalName canonicalName;
		SL2<Complex> x;
		SL2<Complex> y;
		SL2<Complex> X;
		SL2<Complex> Y;
        Word x_word;
        Word y_word;
        double exp_re_perp;
		bool m_foundGoodPair;
		vector<Word> words;
		priority_queue<WordPair, vector<WordPair>, PairComp> pairs;
	};
	
	void TubeSearch::findNames(Word& word)
	{
		word.name = canonicalName.getCanonicalName(word.name);
	}
	
	Word inverse(Word& word)
	{
		Word inv;
		inv.name = word.name;
		reverse(inv.name.begin(), inv.name.end());
		string::size_type pos;
		for (pos = 0; pos < inv.name.size(); ++pos) {
			char c = inv.name[pos];
			if (c >= 'a' && c <= 'z')
				inv.name[pos] += 'A' - 'a';
			else
				inv.name[pos] -= 'A' - 'a';
		}
		inv.matrix = inverse(word.matrix);
		return inv;
	}
	
	Word operator * (const Word& lhs, const Word& rhs)
	{
		Word result;
		result.name = lhs.name + rhs.name;
		result.matrix = lhs.matrix * rhs.matrix;
		return result;
	}

	WordPair::WordPair(Word& first, Word& second, const Params<Complex>& params):first(first), second(second)
	{
//  fprintf(stderr, "Building word pair\n");
//  fprintf(stderr, "a_1 = %f + %f i, c_1 = %f + %f i\n a_2 =  %f + %f i, c_2 = %f + %f i\n",
//                   first.matrix.a.re, first.matrix.a.im, first.matrix.c.re, first.matrix.c.im,
//                   second.matrix.a.re, second.matrix.a.im, second.matrix.c.re, second.matrix.c.im); 
    if (first.name.length() == 0 || second.name.length() == 0) {
      ortho_dist_LB = 1;
      ortho_dist_UB = 1;
    } else {
      ortho_dist_LB = e_re_perp_LB(first.matrix, second.matrix);  
      ortho_dist_UB = e_re_perp_UB(first.matrix, second.matrix);  
    }
    jorgensen_ww = absUB(jorgensen(first.matrix, second.matrix));
   
    if (y_power(first.name) > 0) { 
      jorgensen_xw1 = absUB(jorgensen_xw(first.matrix, params));
      jorgensen_w1x = absUB(jorgensen_wx(first.matrix, params));
    }
    if (y_power(second.name) > 0) { 
      jorgensen_xw2 = absUB(jorgensen_xw(second.matrix, params));
      jorgensen_w2x = absUB(jorgensen_wx(second.matrix, params));
    }
    if (x_power(first.name) > 0) { 
    jorgensen_yw1 = absUB(jorgensen_yw(first.matrix, params));
    jorgensen_w1y = absUB(jorgensen_wy(first.matrix, params));
    }
    if (x_power(second.name) > 0) { 
      jorgensen_yw2 = absUB(jorgensen_yw(second.matrix, params));
      jorgensen_w2y = absUB(jorgensen_wy(second.matrix, params));
    }
    fprintf(stderr, "Built pair (%s,%s) with ortho (%f,%f) and jorgensen:\nw1w2 %f\nxw1 %f\nxw2 %f\nyw1 %f\nyw2 %f\nw1x %f\nw2x %f\nw1y %f\nw2y %f\n",
      first.name.c_str(), second.name.c_str(), ortho_dist_LB,  ortho_dist_UB, jorgensen_ww,
      jorgensen_xw1, jorgensen_xw2, jorgensen_yw1, jorgensen_yw2,  
      jorgensen_w1x, jorgensen_w2x, jorgensen_w1y, jorgensen_w2y);
  }
	
	void TubeSearch::pushWord(Word word)
	{
		vector<Word>::iterator it;
//  first, check to see if this word has already been pushed
		findNames(word);
		for (it = words.begin(); it != words.end(); ++it) {
			if (it->name == word.name || it->name == inverse(word).name) {
				return;
      }
    }
    Complex a = word.matrix.a;
    Complex b = word.matrix.b;
    Complex c = word.matrix.c;
    Complex d = word.matrix.d;
    fprintf(stderr, "Word: [%s]\n", word.name.c_str());
    fprintf(stderr, "At the center is has coords\n");
    fprintf(stderr, "a: %f + I %f\n", a.real(), a.imag());
    fprintf(stderr, "b: %f + I %f\n", b.real(), b.imag());
    fprintf(stderr, "c: %f + I %f\n", c.real(), c.imag());
    fprintf(stderr, "d: %f + I %f\n", d.real(), d.imag());
//   fprintf(stderr,"pushWord(%s)\n", word.name.c_str());
		for (it = words.begin(); it != words.end(); ++it) {
//    fprintf(stderr,"Trying to print\n");
//    fprintf(stderr, "First word %s, Second word %s\n", it->name.c_str(), word.name.c_str());
		  WordPair pair(*it, word, params);
      if (pair.jorgensen_ww < 2) {
        fprintf(stderr,"pushWordPair(%s,%s) with jorgensen %f\n", pair.first.name.c_str(), pair.second.name.c_str(), pair.jorgensen_ww);
				pairs.push(pair);
      }
      if (pair.jorgensen_w1x < 2 || pair.jorgensen_xw1 < 2) {
        WordPair p(pair.first, x_word, params);
        pairs.push(p);
      }
      if (pair.jorgensen_w2x < 2 || pair.jorgensen_xw2 < 2) {
        WordPair p(pair.second, x_word, params);
        pairs.push(p);
      }
      if (pair.jorgensen_w1y < 2 || pair.jorgensen_yw1 < 2) {
        WordPair p(pair.first, y_word, params);
        pairs.push(p);
      }
      if (pair.jorgensen_w2y < 2 || pair.jorgensen_yw2 < 2) {
        WordPair p(pair.second, y_word, params);
        pairs.push(p);
      }
			if (pair.ortho_dist_LB > 1.005 && ( pair.ortho_dist_UB < 2*exp_re_perp || pair.ortho_dist_LB > 100000000 )) { // TODO check if this is sane
        fprintf(stderr,"pushWordPair(%s,%s) with ortho dist %f\n", pair.first.name.c_str(), pair.second.name.c_str(), pair.ortho_dist_UB);
				pairs.push(pair);
      }
		}
		words.push_back(word);
//  fprintf(stderr,"Done with pushWord\n");
	}
	
	void TubeSearch::pushWord(string word)
	{
		Word w;
		w.name = word;
    findNames(w);
    if (w.name == "") return;
    w.matrix = construct_word(word, params);
//  printf("pushWord(%s): distance=%f\n", w.name.c_str(), 0.5 / norm(w.matrix.c));
		pushWord(w);
//	Word wInv = inverse(w);
//	findNames(wInv);
//	pushWord(wInv);
//  fprintf(stderr, "pushed words %s(%s) and %s(%s)\n",
//  w.name.c_str(), w.nameClass.c_str(),
//  wInv.name.c_str(), wInv.nameClass.c_str());
	}
	
	WordPair TubeSearch::findPair()
	{
		while (pairs.size() > 0) {
			WordPair smallest(pairs.top());
      fprintf(stderr, "Smallest (%s, %s)\n", smallest.first.name.c_str(), smallest.second.name.c_str());
//      fprintf(stderr, "adding %s\n", (smallest.first*smallest.second).name.c_str());
      pushWord(inverse(smallest.first) * smallest.second);
//      fprintf(stderr, "adding %s\n", (inverse(smallest.first)*smallest.second).name.c_str());
      pushWord(smallest.first * inverse(smallest.second));
//      fprintf(stderr, "adding %s\n", (smallest.first*inverse(smallest.second)).name.c_str());
      pushWord(smallest.first * smallest.second);
      if (smallest.ortho_dist_UB < exp_re_perp || smallest.ortho_dist_LB > 1000000 || smallest.jorgensen_ww < 1 ||
          smallest.jorgensen_xw1 < 1 || smallest.jorgensen_xw2 < 1 || smallest.jorgensen_yw1 < 1 || smallest.jorgensen_yw2 < 1 ||  
          smallest.jorgensen_w1x < 1 || smallest.jorgensen_w2x < 1 || smallest.jorgensen_w1y < 1 || smallest.jorgensen_w2y < 1 ) {
//        fprintf(stderr, "Pair (%s,%s) with ortho (%f,%f) and jorgensen:\nw1w2 %f\nxw1 %f\nxw2 %f\nyw1 %f\nyw2 %f\nw1x %f\nw2x %f\nw1y %f\nw2y %f\n",
//          smallest.first.name.c_str(), smallest.second.name.c_str(), smallest.ortho_dist_LB,  smallest.ortho_dist_UB, smallest.jorgensen,
//          smallest.jorgensen_xw1, smallest.jorgensen_xw2, smallest.jorgensen_yw1, smallest.jorgensen_yw2,  
//          smallest.jorgensen_w1x, smallest.jorgensen_w2x, smallest.jorgensen_w1y, smallest.jorgensen_w2y);
        m_foundGoodPair = true;
        pairs.pop();
	  	  return smallest;
			}
      return pairs.top();
//      pairs.pop();
//      fprintf(stderr, "rejecting %s^-1 *(%d,%d) %s = %s -> %s (%s)\n",
//        largest.first->name.c_str(),
//        largest.distances[largest.distanceIndex-1].x,
//        largest.distances[largest.distanceIndex-1].y,
//        largest.second->name.c_str(),
//        oldName.c_str(),
//        result.name.c_str(),
//        result.nameClass.c_str()
//        );
//      fprintf(stderr, "duplicate; resorting\n");
		}
    Word w;
    return WordPair(w, w, params);
	}
	
	void TubeSearch::addRelator(string w)
	{
		canonicalName.addRelator(w);
	}
}

vector< word_pair > findPairs(Params<Complex> center, vector<string> seedWords, int numWords, int maxLength,
	vector<string> relators)
{
//fprintf(stderr, "Params:\n");
//fprintf(stderr, "L: %f + I %f \n", center.sinhP.re, center.sinhP.im);
//fprintf(stderr, "S: %f + I %f \n", center.sinhD2.re, center.sinhD2.im);
//fprintf(stderr, "P: %f + I %f \n", center.sinhL2.re, center.sinhL2.im);
  TubeSearchImpl::TubeSearch search(center);
//fprintf(stderr, "Params:\n");
//fprintf(stderr, "L: %f + I %f \n", search.params.sinhP.re, search.params.sinhP.im);
//fprintf(stderr, "S: %f + I %f \n", search.params.sinhD2.re, search.params.sinhD2.im);
//fprintf(stderr, "P: %f + I %f \n", search.params.sinhL2.re, search.params.sinhL2.im);
	for (int i = 0; i < relators.size(); ++i)
		if (relators[i].size() > 0)
			search.addRelator(relators[i]);
	for (int i = 0; i < seedWords.size(); ++i)
		if (seedWords[i].size() > 0)
			search.pushWord(seedWords[i]);
	// search.pushWord("");
	search.pushWord("x");
	search.pushWord("y");
	search.pushWord("YXyyxxyxYYXX");
	set< word_pair > foundPairs;
	while (numWords > int(search.words.size()) || (-numWords > int(search.words.size()) && !search.foundGoodPair())) {
    fprintf(stderr, "While loop: words size %d\n", (int) search.words.size());
		TubeSearchImpl::WordPair p = search.findPair();
		if (p.first.name.length() > maxLength || p.second.name.length() > maxLength) {
			search.m_foundGoodPair = false;
			continue;
		}
    word_pair P(p.first.name, p.second.name);
		foundPairs.insert(P);
	}
  vector< word_pair > v(foundPairs.begin(), foundPairs.end());
	return v;
}

// 	Complex centerDiff = (first.matrix.a / first.matrix.c
// 	 - second.matrix.a / second.matrix.c);
//        if (absUB(centerDiff) < infinity()) {
//            int y0 = int(floor(centerDiff.im / params.sinhP.im));
//            int x0 = int(floor(centerDiff.re - y0*params.sinhP.re));
//            fprintf(stderr, "%d, %d\n", x0, y0); 
//            for (int x = -2-x0; x <= 3-x0; ++x) {
//                for (int y = -1-y0; y <= 2-y0; ++y) {
//                    // this is here to avoid having to reduce first^(-1) * second
//                    if (x == 0 && y == 0 && first.name[0] == second.name[0])
//                        continue;
//                    LatticePoint l;
//                    l.params = params;
//                    l.x = -x;
//                    l.y = -y;
//                    Complex t = (centerDiff + double(x) + double(y)*params.sinhP);
//                    l.distance = pow(absUB(t),2);
//                    fprintf(stderr, "%d, %d\n", -x, -y); 
//                    if (abs(l.x) < 4 && abs(l.y) < 4) 
//                        distances.push_back(l);
//                }
//            }
//        }
// 	sort(distances.begin(), distances.end(), LowerDistance());
// 	setCurrentIndex(0);
//	}
	
// void WordPair::setCurrentIndex(int i)
// {
// 	distanceIndex = i;
// 	if (i >= distances.size()) {
// 		ortho_dist_UB = 0;
// 	} else {
// 		currentCombination = (inverse(first) * distances[i]) * second;
//            // Looks like this is an heuristic ball size 
// 		ortho_dist_UB = 0.5 / pow(absLB(currentCombination.matrix.c),2);
// 	}
//	 fprintf(stderr, "WordPair(%s,%s):setCurrentIndex(%d) : %d,%d %f %s size=%f\n",
//			first->name.c_str(), second->name.c_str(),
//			i, distances[i].x, distances[i].y, distances[i].distance,
//			currentCombination.name.c_str(), ortho_dist_UB);
// }
