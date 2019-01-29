#include <vector>
#include <set>
#include <utility>
#include <string>
#include "types.hh"

std::vector< word_pair > findPairs(
	Params<Complex> center,
	std::vector< std::string > seedWords,
	int numWords,
	int maxLength,
	std::vector< std::string > quasiRelators
);
