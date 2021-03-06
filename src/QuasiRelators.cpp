#include <set>
#include <vector>
#include <algorithm>
#include "QuasiRelators.h"
#include "types.hh"
using namespace std;

string QuasiRelators::getName(string w)
{
	NameStore::iterator it = names.find(w);
	if (it == names.end()) {
		addQuasiRelator(w);
		it = names.find(w);
	}
	return it->second;
}

void QuasiRelators::addQuasiRelator(string w)
{
	if (names.find(w) != names.end())
		return;
//  TODO We kill uniquenes for now as this causes ``large horoball'' failures
//	string W = inverse(w);
//	vector<string> versions;
//	for (string::size_type p = 0; p < w.size(); ++p) {
//		string::size_type prev = (p == 0) ? w.size()-1 : p-1;
//		if (w[prev] == 'g' || w[prev] == 'G')
//			versions.push_back( w.substr(p, string::npos) + w.substr(0, p));
//		if (W[prev] == 'g' || W[prev] == 'G')
//			versions.push_back( W.substr(p, string::npos) + W.substr(0, p));
//	}
//	string first = versions[0];
//	for (vector<string>::iterator it = versions.begin(); it != versions.end(); ++it) {
//		if (*it < first)
//			first = *it;
//	}
//	for (vector<string>::iterator it = versions.begin(); it != versions.end(); ++it)
//		names[*it] = first;
    names[w] = w;
	nameVector.push_back(w);
}

string QuasiRelators::min_pow_desc()
{
    vector<string> qrs(nameVector);
    sort(qrs.begin(), qrs.end(), x_power_sort);
    int min_power = x_power(qrs.front());
	string buf;
	for (vector<string>::iterator it = qrs.begin(); it != qrs.end(); ++it) {
        if (x_power(*it) > min_power) {
            break;
        }
		if (!buf.empty()) {
            buf += ",";
        } 
		buf += *it;
    }
	return buf;
}

string QuasiRelators::desc()
{
	string buf;
	for (vector<string>::iterator it = nameVector.begin(); it != nameVector.end(); ++it) {
		if (!buf.empty())
			buf += ",";
		buf += *it;
	}
	return buf;
}
	std::vector<std::string> allWords();
	std::vector<std::string> wordClasses();

vector<string> QuasiRelators::wordClasses()
{
	return nameVector;
}

vector<string> QuasiRelators::allWords()
{
	vector<string> result;
	for (NameStore::iterator it = names.begin(); it != names.end(); ++it) {
		result.push_back(it->first);
	}
	return result;
}

bool QuasiRelators::isQuasiRelator(std::string w)
{
	return names.find(w) != names.end();
}

string QuasiRelators::inverse(string w)
{
	reverse(w.begin(), w.end()); // reverse order of w
	
	// reverse case of w
	for (string::size_type p = 0; p < w.size(); ++p) {
		switch(w[p]) {
			case 'm': w[p] = 'M'; break;
			case 'M': w[p] = 'm'; break;
			case 'n': w[p] = 'N'; break;
			case 'N': w[p] = 'n'; break;
			case 'g': w[p] = 'G'; break;
			case 'G': w[p] = 'g'; break;
		}
	}

	// commute m and n so that m comes first
	for (bool progress = true; progress; progress = false) {
		for (string::size_type p = 0; p < w.size(); ++p) {
			if ((w[p  ] == 'n' || w[p  ] == 'N')
			 && (w[p+1] == 'm' || w[p+1] == 'M')) {
				char n = w[p];
				w[p] = w[p+1];
				w[p+1] = n;
				progress = true;
			}
		}
	}
	return w;
}
