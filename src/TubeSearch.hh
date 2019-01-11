#include <vector>
#include <set>
#include <utility>
#include <string>
#include "types.hh"

std::set< std::pair<std::string, std::string> > findPairs(
	Params<Complex> center,
	std::vector< std::string > seedWords,
	int numWords,
	int maxLength,
	std::vector< std::string > quasiRelators
);
