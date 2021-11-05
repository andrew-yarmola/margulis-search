/*
 *  RelatorTest.cpp
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "RelatorTest.hh"
#include "types.hh"
#include <cstdio>

using namespace std;


bool RelatorTest::is_good(string word)
{
  if (word.length() == 0 ||
      bad_relators.find(word) != bad_relators.end()) {
    return false;
  }
  return true;
}

bool RelatorTest::is_impossible(string word,
    vector<string>& required_non_identities)
	{
    required_non_identities.clear();
    if (syllables(word) < 5 ||
        always_impossible.find(word) != always_impossible.end()) {
      return true;
    }
    string cycle(word);
    int rot = 1;
    // could be optimized
    while (rot != word.length()) {
      rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
      if (cycle == word) {
        string sub = word.substr(0, rot);
        if (syllables(sub) < 5) {
          return true;
        } else {
          required_non_identities.push_back(sub);
		      return false;
        }
      }
      ++rot;
    }
		return false;
	}

void RelatorTest::load(const char* path)
{
  char buf[1000];
  char word_buf[1000];
  FILE* fp = fopen(path, "r");
  while (fp && fgets(buf, sizeof(buf), fp)) {
    buf[strcspn(buf, "\r\n")] = 0;
    always_impossible.insert(string(buf));
  }
}

RelatorTest* RelatorTest::create(const char* file_path)
{
	RelatorTest* impossible = new RelatorTest();
	impossible->load(file_path);
	return impossible;
}
