#ifndef __QUASI_RELATORS_H
#define __QUASI_RELATORS_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include "types.hh"
#include "SL2.hh"
#include "CanonicalName.hh"

class QuasiRelators {
public:
	std::string get_name(std::string w);               // get the canonical name of a quasi-relator
	void add_quasi_relator(std::string w);        // record that this word is a quasi-relator

	std::vector<std::string> all_words();
	std::vector<std::string> word_classes();
	std::string desc();                               // string describing this set of quasi-relators
	std::string min_pow_desc();                               // string describing minimal power set of quasi-relators
	bool is_quasi_relator(std::string w);         // is this word a quasi-relator?

  template<typename T>
  std::string desc(const Params<T>& p);

private:
  CanonicalName canonical_name;
	typedef std::map< std::string, std::string > NameStore;
	NameStore names;
	std::vector<std::string> name_vector;
	std::string inverse(std::string w);
};

#define MAX_ID_SHIFT 5

template<typename T>
std::string likely_identity(std::string word, const Params<T>& p) {
  SL2<T> w = construct_word(word, p);
  if (inside_var_nbd_x(w, p)) {
    SL2<T> x = construct_x(p);
    std::string new_word = x_strip(word);
    for (int i = 0; i < MAX_ID_SHIFT; ++i) {
      new_word = "x" + new_word;
      SL2<T> new_w = construct_word(new_word, p); // order matters
      //if (absUB(jorgensen_wx(new_w, p)) < 0.5) {
      if (absUB(four_cosh_re_length(new_w)) < absLB(four_cosh_re_length(x))) {
        return new_word;
      }      
    }
    new_word = x_strip(word);
    for (int i = 0; i < MAX_ID_SHIFT; ++i) {
      new_word = "X" + new_word;
      SL2<T> new_w = construct_word(new_word, p); // order matters
      // if (absUB(jorgensen_wx(new_w, p)) < 0.5) {
      if (absUB(four_cosh_re_length(new_w)) < absLB(four_cosh_re_length(x))) {
        return new_word;
      }      
    }
  }
  if (inside_var_nbd_y(w, p)) {
    SL2<T> y = construct_y(p);
    std::string new_word = y_strip(word);
    for (int i = 0; i < MAX_ID_SHIFT; ++i) {
      new_word = "y" + new_word;
      SL2<T> new_w = construct_word(new_word, p); // order matters
      // if (absUB(jorgensen_wy(new_w, p)) < 0.5) {
      if (absUB(four_cosh_re_length(new_w)) < absLB(four_cosh_re_length(y))) {
        return new_word;
      }      
    }
    new_word = y_strip(word);
    for (int i = 0; i < MAX_ID_SHIFT; ++i) {
      new_word = "Y" + new_word;
      SL2<T> new_w = construct_word(new_word, p); // order matters
      // if (absUB(jorgensen_wy(new_w, p)) < 0.5) {
      if (absUB(four_cosh_re_length(new_w)) < absLB(four_cosh_re_length(y))) {
        return new_word;
      }      
    }
  }
  return "";
}

template<typename T>
std::string QuasiRelators::desc(const Params<T>& p)
{
  std::string buf;
  std::string word;
  std::set<std::string> words;
	for (std::vector<std::string>::iterator it = name_vector.begin(); it != name_vector.end(); ++it) {
		if (!buf.empty() && buf.back() != ',')
			buf += ",";
    word = likely_identity(*it, p);
    word = canonical_name.get_canonical_name(word);
    if (word.length() > 0) {
      if (words.insert(word).second) {
        buf += word;
      }
    }
	}
	return buf;
}


#endif
