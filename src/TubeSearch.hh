#include <vector>
#include <set>
#include <utility>
#include <string>
#include "types.hh"

std::vector< word_pair > find_pairs(
	Params<Complex> center,
	std::vector< std::string > seed_words,
	int num_words,
	int max_length,
	std::vector< std::string > quasi_relators
);
