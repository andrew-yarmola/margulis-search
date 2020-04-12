#ifndef __QUASI_RELATORS_H
#define __QUASI_RELATORS_H

#include <string>
#include <map>
#include <vector>

class QuasiRelators {
public:
	std::string get_name(std::string w);               // get the canonical name of a quasi-relator
	void add_quasi_relator(std::string w);        // record that this word is a quasi-relator

	std::vector<std::string> all_words();
	std::vector<std::string> word_classes();
	std::string desc();                               // string describing this set of quasi-relators
	std::string min_pow_desc();                               // string describing minimal power set of quasi-relators
	bool is_quasi_relator(std::string w);         // is this word a quasi-relator?

private:
	typedef std::map< std::string, std::string > NameStore;
	NameStore names;
	std::vector<std::string> name_vector;
	std::string inverse(std::string w);
};

#endif
