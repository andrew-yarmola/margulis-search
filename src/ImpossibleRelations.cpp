/*
 *  ImpossibleRelations.cpp
 *  mom
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "ImpossibleRelations.h"
#include <map>
#include <stdio.h>
#include <string.h>

using namespace std;

namespace ImpossibleRelationsImpl {
	
	struct Impl : public ImpossibleRelations {
		bool is_impossible(string word, vector<string>& required_non_identities);
		void load(const char* path);
	private:
		struct PossiblePower {
			bool subword_identity_allowed;
			int power;
			string subword;
		};
		typedef multimap<string, PossiblePower> PossibleStore;
		PossibleStore possible_store;
	};
	
	bool Impl::is_impossible(string word, vector<string>& required_non_identities)
	{
		PossibleStore::iterator it = possible_store.lower_bound(word);
		vector<string> required;
		while (it != possible_store.end() && it->first == word) {
      return true;
      /*PossiblePower possible = it->second;
      if (!possible.subword_identity_allowed) {
        return true;
      } else {
        required.push_back(possible.subword);
      }
			++it;*/
		}
		required_non_identities.swap(required);
		return false;
	}

	void Impl::load(const char* path)
	{
		char buf[1000];
		char word_buf[1000];
		char sub_word_buf[1000];
		int subword_identity_allowed;
		int matchRequired;
		PossiblePower possible;
		FILE* fp = fopen(path, "r");
		while (fp && fgets(buf, sizeof(buf), fp)) {
      buf[strcspn(buf, "\r\n")] = 0;
      possible_store.insert(make_pair(string(buf), possible));
      /*
			int n = sscanf(buf, "PossiblePower %s %d %d %d %d %[gGmMnN]^%d",
				word_buf, &subword_identity_allowed, &matchRequired,
				&possible.matchingMCoeff, &possible.matchingNCoeff,
				sub_word_buf, &possible.power);
			if (n == 7) { // Filled all the values
				possible.subword_identity_allowed = subword_identity_allowed;
				possible.matchRequired = matchRequired;
				possible.subWord = sub_word_buf;
				possible_store.insert(make_pair(string(word_buf), possible));
			} else {
				if (n > 0) fprintf(stderr, "incomplete line %s", buf);
				return;
			}*/
		}
	}
}

ImpossibleRelations* ImpossibleRelations::create(const char* file_path)
{
	ImpossibleRelationsImpl::Impl* impossible = new ImpossibleRelationsImpl::Impl();
	impossible->load(file_path);
	return impossible;
}
