/*
 *  RelatorTest.h
 *  mom
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <vector>
#include <set>

struct RelatorTest
{
public:
	// if the return value is false, sets mandatory to
  // the list of subwords which should not be identies.
	bool is_impossible(std::string word,
      std::vector<std::string>& required_non_identities);
	bool is_good(std::string word);
	static RelatorTest* create(const char* file);
private:
  void load(const char* path);
  std::set<std::string> always_impossible;
  const std::set<std::string> bad_relators = 
  {
    "XXXYxxYYYYxxY", // m142
    "XXXYXYYxYYXY", // non-realizable
    "YYXXYXyxyx",   // non-realizable
    "xxYYxYXyXy",  // non-realizable
    "YYxyxyXYXX", // non-realizable
    "XXyxyxYXYY", // non-realizable
    "yyyyxYXXYx", // m009
    "yyyXXYXYXX", // m026
    "yyyxxYxYxx", // m026
    "xYYxxYxxYY", // m003 too symmetric
    "YYxxYxxYYx", // m003 too symmetric
    "YYXXYXXYYX", // m003 too symmetric
    "XYYXXYXXYY", // m003 too symmetric
  };
};
