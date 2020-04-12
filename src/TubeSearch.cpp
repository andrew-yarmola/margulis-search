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

  inline bool pos_less(double a, double b) {
    return a > 0 && a < b;
  }

  inline bool pos_neq(double a, double b) {
    return a > 0 && b > 0 && a != b;
  }
	
	struct Word {
		string name;
		SL2<Complex> matrix;
    void compute_vals(const Params<Complex>& params);
    double four_cosh_x_x = -1;
    double four_cosh_x_y = -1;
    double four_cosh_y_y = -1;
    double four_cosh_y_x = -1;
    double cosh_move = -1;
    double jorg_x = -1;
    double jorg_y = -1;
    bool operator<(const Word& b) const {
      if (x_power(name) > 0 && y_power(name) > 0 && x_power(b.name) > 0 && y_power(b.name) > 0) {
        if (pos_neq(four_cosh_x_x, b.four_cosh_x_x)) {
          return four_cosh_x_x > b.four_cosh_x_x;
        }
        if (pos_neq(four_cosh_y_y, b.four_cosh_y_y)) {
          return four_cosh_y_y > b.four_cosh_y_y;
        }
        if (pos_neq(four_cosh_x_y, b.four_cosh_x_y)) {
          return four_cosh_x_y > b.four_cosh_x_y;
        }
        if (pos_neq(four_cosh_y_x, b.four_cosh_y_x)) {
          return four_cosh_y_x > b.four_cosh_y_x;
        }
        if (pos_neq(cosh_move, b.cosh_move)) {
          return cosh_move > b.cosh_move;
        }
        if (pos_neq(jorg_x, b.jorg_x)) {
          return jorg_x > b.jorg_x;
        }
        if (pos_neq(jorg_y, b.jorg_y)) {
          return jorg_y > b.jorg_y;
        }
      } else {
        if (x_power(name) > 0 && y_power(name) > 0) {
          return false;
        }
        if (x_power(b.name) > 0 && y_power(b.name) > 0) {
          return true;
        }
      }
      return name > b.name;
    }
	};
	
	struct WordPair {
		WordPair(Word& first, Word& second, const Params<Complex>& params);
		Word first;
		Word second;
		double cosh_two_ortho = -1;
    double four_cosh_marg = -1;
    double jorgensen_words = -1;
	};

	struct PairComp : public binary_function<const WordPair&, const WordPair&, bool>
	{
    bool operator() (const WordPair& a, const WordPair& b)
		{
      if (pos_neq(a.cosh_two_ortho, b.cosh_two_ortho)) {
        return a.cosh_two_ortho < b.cosh_two_ortho;
      }
      if (pos_neq(a.four_cosh_marg, b.four_cosh_marg)) {
        return a.four_cosh_marg < b.four_cosh_marg;
      }
      if (pos_neq(a.jorgensen_words, b.jorgensen_words)) {
        return a.jorgensen_words < b.jorgensen_words;
      }
      if (a.first.name != b.first.name) {
        return a.first < b.first;
      }
      return a.second < b.second;
    }
	};

  string name_inverse(string name) {
    string inv = name;
		reverse(inv.begin(), inv.end());
		string::size_type pos;
		for (pos = 0; pos < inv.size(); ++pos) {
			char c = inv[pos];
			if (c >= 'a' && c <= 'z')
				inv[pos] += 'A' - 'a';
			else
				inv[pos] -= 'A' - 'a';
		}
    return inv;
  }

	Word inverse(Word& word, const Params<Complex> params)
	{
		Word inv;
		inv.name = name_inverse(word.name);
		inv.matrix = inverse(word.matrix);
    inv.compute_vals(params);
		return inv;
	}
	
	struct TubeSearch {
		TubeSearch(Params<Complex> params_) :params(params_), m_found_good_pair(false) {
      x_word.name = "x";
      x_word.matrix = construct_x(params);
      x_word.compute_vals(params);
      y_word.name = "y";
      y_word.matrix = construct_y(params);
      y_word.compute_vals(params);
      X_word = inverse(x_word, params);
      Y_word = inverse(y_word, params);
    }

		~TubeSearch() {}
		bool push_word(Word word);
		bool push_word(string name);
		WordPair find_pair();
		bool found_good_pair() { return m_found_good_pair; }
		void add_relator(string word);
		
		const Params<Complex> params;
		CanonicalName canonical_name;
    Word x_word;
    Word y_word;
    Word X_word;
    Word Y_word;
		bool m_found_good_pair;
		vector<Word> word_lookup;
		priority_queue<Word> words;
		priority_queue<WordPair, vector<WordPair>, PairComp> pairs;
	};
	
	Word product(const Word& lhs, const Word& rhs, const Params<Complex>& params, CanonicalName& canonical_name)
	{
		Word result;
		result.name = canonical_name.get_canonical_name(lhs.name + rhs.name);
    if (result.name == lhs.name + rhs.name) {
		  result.matrix = lhs.matrix * rhs.matrix;
    } else {
      result.matrix = construct_word(result.name, params);
    }
    result.compute_vals(params);
		return result;
	}

	void Word::compute_vals(const Params<Complex>& params)
	{
    if (y_power(name) > 0) {
      four_cosh_x_x = absUB(four_cosh_dist_ax_wax(matrix, params)); 
      four_cosh_y_x = absUB(four_cosh_dist_ay_wax(matrix, params)); 
      jorg_x = min(absUB(jorgensen_xw(matrix, params)),
                         absUB(jorgensen_wx(matrix, params))); 
      if (x_power(name) > 0) {
        cosh_move = absUB(cosh_move_j(matrix));
      }
    }
    if (x_power(name) > 0) {
      four_cosh_y_y = absUB(four_cosh_dist_ay_way(matrix, params)); 
      four_cosh_x_y = absUB(four_cosh_dist_ax_way(matrix, params)); 
      jorg_y = min(absUB(jorgensen_yw(matrix, params)),
                         absUB(jorgensen_wy(matrix, params))); 
    }
  }

	WordPair::WordPair(Word& first, Word& second, const Params<Complex>& params):first(first), second(second)
	{
//  fprintf(stderr, "Building word pair\n");
//  fprintf(stderr, "a_1 = %f + %f i, c_1 = %f + %f i\n a_2 =  %f + %f i, c_2 = %f + %f i\n",
//                   first.matrix.a.re, first.matrix.a.im, first.matrix.c.re, first.matrix.c.im,
//                   second.matrix.a.re, second.matrix.a.im, second.matrix.c.re, second.matrix.c.im); 
    if (first.name.length() != 0 && second.name.length() != 0) {
      cosh_two_ortho = absUB(cosh_2_re_perp(first.matrix, second.matrix));  
      jorgensen_words = absUB(jorgensen(first.matrix, second.matrix));
      four_cosh_marg = absUB(four_cosh_margulis_simple(first.matrix, second.matrix).first);
    }
//    fprintf(stderr, "Built pair (%s,%s)\n",
//      first.name.c_str(), second.name.c_str());
  }
	
	bool TubeSearch::push_word(Word word)
	{
    if (word.name == "") { return false; }
		for (auto it = word_lookup.begin(); it != word_lookup.end(); ++it) {
			if (it->name == word.name) {
				return false;
      }
    }
    // We assume compute_vals has run on the word
    // fprintf(stderr, "Word: %s\n", word.name.c_str());
    bool pushed = false;
  
    if (pos_less(word.four_cosh_x_x, absUB(params.cosh2dx * 6.5))   ||  
        pos_less(word.four_cosh_y_y, absUB(params.cosh2dy * 6.5))   ||
        pos_less(word.four_cosh_x_y, absUB(params.coshdxdy * 6.5))  ||
        pos_less(word.four_cosh_y_x, absUB(params.coshdxdy * 6.5))  ||
        pos_less(word.cosh_move, absUB(params.coshmu * 2.5))    ||
        pos_less(word.jorg_x, 3)                                ||
        pos_less(word.jorg_y, 3))
    {
	  	word_lookup.push_back(word);
	  	words.push(word);
      pushed = true;
      //fprintf(stderr, "Pushed word: %s\n", word.name.c_str());
    }

    if (x_power(word.name) == 0 || y_power(word.name) == 0) {
      return pushed;
    }

//    return pushed;

		for (auto it = word_lookup.begin(); it != word_lookup.end(); ++it) {
//    fprintf(stderr,"Trying to print\n");
//    fprintf(stderr, "First word %s, Second word %s\n", it->name.c_str(), word.name.c_str());
      if (word.name == it->name || word.name == name_inverse(it->name) || x_power(it->name) == 0 || y_power(it->name) == 0) {
        continue;
      }
		  WordPair pair(word, *it, params);
			if (pos_less(pair.jorgensen_words, 1.5)                        ||
          pos_less(pair.four_cosh_marg, absUB(params.coshmu * 4.5))    ||
          pos_less(pair.cosh_two_ortho, absUB(params.coshdxdy * 1.2)))
      {
				pairs.push(pair);
        // pushed = true;
        // fprintf(stderr,"Pushed pair (%s,%s)\n", pair.first.name.c_str(), pair.second.name.c_str());
      }
		}
    return pushed;
//  fprintf(stderr,"Done with push_word\n");
	}
	
	bool TubeSearch::push_word(string name)
	{
	  name = canonical_name.get_canonical_name(name);
    if (name == "") return false;
		Word w;
		w.name = name;
    w.matrix = construct_word(w.name, params);
    w.compute_vals(params);
		return push_word(w);
//	Word wInv = inverse(w);
//	find_names(wInv);
//	push_word(wInv);
//  fprintf(stderr, "pushed words %s(%s) and %s(%s)\n",
//  w.name.c_str(), w.nameClass.c_str(),
//  wInv.name.c_str(), wInv.nameClass.c_str());
	}

#define SEARCH_DEPTH 100
#define POW_DEPTH 6

	WordPair TubeSearch::find_pair()
	{
    int loop_count = 0;
		while ((words.size() > 0 || pairs.size() > 0) && loop_count < SEARCH_DEPTH) {
      if (words.size() > 0) {
        Word w(words.top());
        // fprintf(stderr, "Smallest word %s\n", w.name.c_str());
//        if (w.name == "Xy") {
//          fprintf(stderr, "    4cosh(ax, wax): %f < %f\n    4cosh(ay, way): %f < %f\n    4cosh(ax, way): %f < %f\n    4cosh(ay, wax): %f < %f\n    cosh_move: %f < %f\n    jorg_x: %f < 1.0\n    jorg_y: %f < 1.0\n", w.four_cosh_x_x, absLB(params.cosh2dx * 4), w.four_cosh_y_y, absLB(params.cosh2dy * 4), w.four_cosh_x_y, absLB(params.coshdxdy * 4), w.four_cosh_y_x, absLB(params.coshdxdy * 4), w.cosh_move, absLB(params.coshmu), w.jorg_x, w.jorg_y);
//        } 
        if (y_power(w.name) > 0 && x_power(w.name) > 0 &&
            (pos_less(w.four_cosh_x_x, absLB(params.cosh2dx * 3.999999))   ||  
             pos_less(w.four_cosh_y_y, absLB(params.cosh2dy * 3.999999))   ||
             pos_less(w.four_cosh_x_y, absLB(params.coshdxdy * 3.999999))  ||
             pos_less(w.four_cosh_y_x, absLB(params.coshdxdy * 3.999999))  ||
             pos_less(w.cosh_move, absLB(params.coshmu * 0.999999))        ||
             pos_less(w.jorg_x, 0.9999999)                                 ||
             pos_less(w.jorg_y, 0.9999999)))
        {
          m_found_good_pair = true;
          Word none;
          none.name = "";
          WordPair result(w, none, params);
          // fprintf(stderr, "Found word %s\n", w.name.c_str());
          // fprintf(stderr, "    4cosh(ax, wax): %f < %f\n    4cosh(ay, way): %f < %f\n    4cosh(ax, way): %f < %f\n    4cosh(ay, wax): %f < %f\n    cosh_move: %f < %f\n    jorg_x: %f < 1.0\n    jorg_y: %f < 1.0\n", w.four_cosh_x_x, absLB(params.cosh2dx * 4), w.four_cosh_y_y, absLB(params.cosh2dy * 4), w.four_cosh_x_y, absLB(params.coshdxdy * 4), w.four_cosh_y_x, absLB(params.coshdxdy * 4), w.cosh_move, absLB(params.coshmu), w.jorg_x, w.jorg_y);
          return result;
        }
        words.pop();
        if (w.name[0] == 'y' || w.name[0] == 'Y') {
          int x_pow = 0;
          Word x_mult = product(x_word, w, params, canonical_name);
          while(push_word(x_mult) && x_pow < POW_DEPTH) {
            x_mult = product(x_word, x_mult, params, canonical_name);
            x_pow++;
          }
          int X_pow = 0;
          Word X_mult = product(X_word, w, params, canonical_name);
          while(push_word(X_mult) && X_pow < POW_DEPTH) {
            X_mult = product(X_word, X_mult, params, canonical_name);
            X_pow++;
          }
        } else {
          int y_pow = 0;
          Word y_mult = product(y_word, w, params, canonical_name);
          while(push_word(y_mult) && y_pow < POW_DEPTH) {
            y_mult = product(y_word, y_mult, params, canonical_name);
            y_pow++;
          }
          int Y_pow = 0;
          Word Y_mult = product(Y_word, w, params, canonical_name);
          while(push_word(Y_mult) && Y_pow < POW_DEPTH) {
            Y_mult = product(Y_word, Y_mult, params, canonical_name);
            Y_pow++;
          }
        }
      }
      if (pairs.size() > 0) { 
        WordPair p(pairs.top());
        // fprintf(stderr, "Smallest pair (%s, %s)\n", p.first.name.c_str(), p.second.name.c_str());
        if (pos_less(p.jorgensen_words, 0.999999)                ||
            pos_less(p.four_cosh_marg, absLB(params.coshmu * 3.999999))) 
            //pos_less(p.cosh_two_ortho, absLB(params.coshdxdy)))
        {
          m_found_good_pair = true;
          pairs.pop();
          // fprintf(stderr, "Found pair (%s, %s)\n", p.first.name.c_str(), p.second.name.c_str());
          // fprintf(stderr, "    jorgensen_words: %f < 1.0\n    four_cosh_marg: %f < %f\n    cosh_two_ortho: %f < %f\n", p.jorgensen_words, p.four_cosh_marg, absLB(params.coshmu * 4), p.cosh_two_ortho, absLB(params.coshdxdy));
          return p;
        }
        pairs.pop();
        push_word(product(inverse(p.first, params), p.second, params, canonical_name));
        //push_word(product(inverse(p.second, params), p.first, params, canonical_name));
        //push_word(product(p.first, inverse(p.second, params), params, canonical_name));
        //push_word(product(p.first, p.second, params, canonical_name));
      }
      loop_count++;
		}
    Word w;
    w.name = "";
    return WordPair(w, w, params);
	}
	
	void TubeSearch::add_relator(string w)
	{
		canonical_name.add_relator(w);
	}
}

#define MAX_ATTEMPTS 2

vector< word_pair > find_pairs(Params<Complex> center, vector<string> seed_words,
  int num_words, int max_length, vector<string> relators)
{
  TubeSearchImpl::TubeSearch search(center);
  search.add_relator("xX"); 
  search.add_relator("Xx"); 
  search.add_relator("yY"); 
  search.add_relator("Yy"); 
	for (int i = 0; i < relators.size(); ++i)
		if (relators[i].size() > 0)
			search.add_relator(relators[i]);
	search.push_word("x");
	search.push_word("y");
	search.push_word("X");
	search.push_word("Y");
	for (int i = 0; i < seed_words.size(); ++i)
		if (seed_words[i].size() > 0)
			search.push_word(seed_words[i]);
	set< word_pair > found_pairs;
  int attempts = 0;
	while (attempts < MAX_ATTEMPTS && (num_words > int(found_pairs.size()) || (-num_words > int(found_pairs.size()) && !search.found_good_pair()))) {
    //fprintf(stderr, "While loop: words size %d\n", (int) found_pairs.size());
		TubeSearchImpl::WordPair p = search.find_pair();
    //if (search.found_good_pair()) {
    //   printf("FOUND PAIR\n");
    //}
    attempts++;
		if (p.first.name.length() == 0 || (p.first.name.length() > max_length && p.second.name.length() > max_length)) {
			search.m_found_good_pair = false;
			continue;
		}
    word_pair P(p.first.name, p.second.name);
		found_pairs.insert(P);
	}
  vector< word_pair > v(found_pairs.begin(), found_pairs.end());
	return v;
}
