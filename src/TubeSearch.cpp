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
    double jorgensen;
    double jorgensen_xw1;
    double jorgensen_xw2;
    double jorgensen_yw1;
    double jorgensen_yw2;
    double jorgensen_w1x;
    double jorgensen_w2x;
    double jorgensen_w1y;
    double jorgensen_w2y;
	};
	
	struct OrthoDistComp : public binary_function<const WordPair&, const WordPair&, bool>
	{
		bool operator() (const WordPair& a, const WordPair& b)
		{
			return a.ortho_dist_UB > b.ortho_dist_UB;
		}
	};
	
	struct TubeSearch {
		TubeSearch(Params<Complex> params_) :params(params_),
               x(construct_x(params)), X(inverse(x)),
               y(construct_y(params)), Y(inverse(y)),
               exp_re_perp(absUB(params.coshP + params.sinhP)), m_foundShortOrtho(false) {}
		~TubeSearch() {}
		void pushWord(Word word);
		void pushWord(string word);
		WordPair findPair();
		bool foundShortOrtho() { return m_foundShortOrtho; }
		void findNames(Word& word);
		void addRelator(string word);
		
		Params<Complex> params;
		CanonicalName canonicalName;
		SL2<Complex> x;
		SL2<Complex> y;
		SL2<Complex> X;
		SL2<Complex> Y;
    double exp_re_perp;
		bool m_foundShortOrtho;
		vector<Word> words;
		priority_queue<WordPair, vector<WordPair>, OrthoDistComp> pairs;
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
//  fprintf(stderr, "Buidling word pair\n");
//  fprintf(stderr, "a_1 = %f + %f i, c_1 = %f + %f i\n a_2 =  %f + %f i, c_2 = %f + %f i\n",
//                   first.matrix.a.re, first.matrix.a.im, first.matrix.c.re, first.matrix.c.im,
//                   second.matrix.a.re, second.matrix.a.im, second.matrix.c.re, second.matrix.c.im); 
    ortho_dist_LB = e_re_perp_LB(first.matrix, second.matrix);  
    ortho_dist_UB = e_re_perp_UB(first.matrix, second.matrix);  
    jorgensen = jorgensen_UB(first.matrix, second.matrix);
    jorgensen_xw1 = jorgensen_xw_UB(first.matrix, params);
    jorgensen_xw2 = jorgensen_xw_UB(second.matrix, params);
    jorgensen_yw1 = jorgensen_yw_UB(first.matrix, params);
    jorgensen_yw2 = jorgensen_yw_UB(second.matrix, params);
    jorgensen_w1x = jorgensen_wx_UB(first.matrix, params);
    jorgensen_w2x = jorgensen_wx_UB(second.matrix, params);
    jorgensen_w1y = jorgensen_wy_UB(first.matrix, params);
    jorgensen_w2y = jorgensen_wy_UB(second.matrix, params);
  }
	
	void TubeSearch::pushWord(Word word)
	{
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
		vector<Word>::iterator it;
//  first, check to see if this word has already been pushed
		findNames(word);
		for (it = words.begin(); it != words.end(); ++it) {
			if (it->name == word.name) {
				return;
      }
    }
    fprintf(stderr,"pushWord(%s)\n", word.name.c_str());
		words.push_back(word);
		Word* wp = &words.back();
		for (it = words.begin(); it != words.end(); ++it) {
//    fprintf(stderr,"Trying to print\n");
      if (it->name == word.name || it->name == inverse(word).name) continue; 
//      fprintf(stderr, "First word %s, Second word %s\n", it->name.c_str(), word.name.c_str());
		  WordPair pair(*it, word, params);
			if (pair.ortho_dist_LB > 1.05 && ( pair.ortho_dist_UB < 2*exp_re_perp || pair.ortho_dist_LB > 100000000 )) { // TODO check if this is sane
        fprintf(stderr,"pushWordPair(%s,%s) with ortho dist %f\n", pair.first.name.c_str(), pair.second.name.c_str(), pair.ortho_dist_UB);
				pairs.push(pair);
      }
      if (pair.jorgensen < 4 ||
          pair.jorgensen_xw1 < 4 || pair.jorgensen_xw2 < 4 || pair.jorgensen_yw1 < 4 || pair.jorgensen_yw2 < 4 ||  
          pair.jorgensen_w1x < 4 || pair.jorgensen_w2x < 4 || pair.jorgensen_w1y < 4 || pair.jorgensen_w2y < 4 ) {
        fprintf(stderr,"pushWordPair(%s,%s) with ortho dist %f\n", pair.first.name.c_str(), pair.second.name.c_str(), pair.ortho_dist_UB);
				pairs.push(pair);
      }
		}
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
      if (smallest.ortho_dist_UB < exp_re_perp || smallest.ortho_dist_LB > 1000000 || smallest.jorgensen < 1.5 ||
          smallest.jorgensen_xw1 < 1.5 || smallest.jorgensen_xw2 < 1.5 || smallest.jorgensen_yw1 < 1.5 || smallest.jorgensen_yw2 < 1.5 ||  
          smallest.jorgensen_w1x < 1.5 || smallest.jorgensen_w2x < 1.5 || smallest.jorgensen_w1y < 1.5 || smallest.jorgensen_w2y < 1.5 ) {
        fprintf(stderr, "Pair (%s,%s) with ortho (%f,%f) and jorgensen:\nw1w2 %f\nxw1 %f\nxw2 %f\nyw1 %f\nyw2 %f\nw1x %f\nw2x %f\nw1y %f\nw2y %f\n",
          smallest.first.name.c_str(), smallest.second.name.c_str(), smallest.ortho_dist_LB,  smallest.ortho_dist_UB, smallest.jorgensen,
          smallest.jorgensen_xw1, smallest.jorgensen_xw2, smallest.jorgensen_yw1, smallest.jorgensen_yw2,  
          smallest.jorgensen_w1x, smallest.jorgensen_w2x, smallest.jorgensen_w1y, smallest.jorgensen_w2y);
        m_foundShortOrtho = true;
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
	set< word_pair > foundPairs;
	while (numWords > int(search.words.size()) || (-numWords > int(search.words.size()) && !search.foundShortOrtho())) {
    fprintf(stderr, "While loop: words size %d\n", (int) search.words.size());
		TubeSearchImpl::WordPair p = search.findPair();
		if (p.first.name.length() > maxLength || p.second.name.length() > maxLength) {
			search.m_foundShortOrtho = false;
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
