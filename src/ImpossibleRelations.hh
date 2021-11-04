/*
 *  ImpossibleRelations.h
 *  mom
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <vector>

struct ImpossibleRelations
{
public:
	// if the return value is false, sets mandatory to the list of subwords which should not be identies.
	virtual bool is_impossible(std::string word, std::vector<std::string>& required_non_identities) = 0;
	static ImpossibleRelations* create(const char* file);
};
