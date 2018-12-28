/*
 *  MomRefine.cpp
 *  mom
 *
 *  Created by Nathaniel Thurston on 10/10/2007.
 *  Copyright 2007 THingith ehf.. All rights reserved.
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Box.h"
#include "TestCollection.h"
#include "BallSearch.h"
#include "QuasiRelators.h"

using namespace std;

struct Options {
	Options() :
        boxName(""), // Binary representation of box
        wordsFile("words"), // Previously generated words
		powersFile("null"), // Output from power parabolic.pl
		maxDepth(18), // Maximum depth for a file
        truncateDepth(6), 
        inventDepth(12),
        maxSize(1000000),
		improveTree(true),
		wordSearchDepth(-1),
		fillHoles(true),
		maxWordLength(20) {}
	const char* boxName;
	const char* wordsFile;
	const char* powersFile;
	int maxDepth;
	int truncateDepth;
	int inventDepth;
	int maxSize;
	bool improveTree;
	int wordSearchDepth;
	bool fillHoles;
	int maxWordLength;
};

Options g_options;
TestCollection g_tests;
typedef vector< vector< box_state > > TestHistory;
static int g_boxesVisited = 0;

struct PartialTree {
	PartialTree() : lChild(NULL), rChild(NULL), testIndex(-1), testResult(open), aux_word(), qr_desc() {}
	PartialTree *lChild;
	PartialTree *rChild;
	int testIndex;
	box_state testResult;
    string aux_word;
    string qr_desc;
};

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree readTree()
{
	PartialTree t;
	char buf[1000];
	if (!fgets(buf, sizeof(buf), stdin)) {
		fprintf(stderr, "unexpected EOF\n");
		exit(1);
	}
	int n = strlen(buf);
	if (buf[n-1] == '\n')
		buf[n-1] = '\0';
	if (buf[0] == 'X') {
		t.lChild = new PartialTree(readTree());
		t.rChild = new PartialTree(readTree());
	} else if (strstr(buf, "HOLE") != NULL) {
		t.testIndex = -2;
	} else {
		if (isdigit(buf[0])) {
			t.testIndex = atoi(buf);
		} else {
            // Add word as eliminator to test collection
            if (strchr("xXyY", buf[0]) != NULL) {
			    t.testIndex = g_tests.add(buf);
            } else {
                if (buf[0] == 'E') {
                    char * comma = strchr(buf,',');
                    comma[0] = '\0'; 
                } else {
                    buf[n-2] = '\0';
                }
                t.testIndex = g_tests.add(&buf[2]);
            }
		}
	}
	return t;
}

void truncateTree(PartialTree& t)
{
	if (t.lChild) {
		truncateTree(*t.lChild);
		delete t.lChild;
		t.lChild = 0;
	}
	if (t.rChild) {
		truncateTree(*t.rChild);
		delete t.rChild;
		t.rChild = 0;
	}
}

int treeSize(PartialTree& t) {
	int size = 1;
	if (t.lChild)
		size += treeSize(*t.lChild);
	if (t.rChild)
		size += treeSize(*t.rChild);
	return size;
}


extern string g_testCollectionFullWord;
extern double g_margulis_bound;

unordered_map<string,SL2ACJ> short_words_cache;

bool refineRecursive(Box box, PartialTree& t, int depth, TestHistory& history, vector< Box >& place, int newDepth, int& searchedDepth)
{
	//fprintf(stderr, "rr: %s depth %d placeSize %lu\n", box.name.c_str(), depth, place.size());
	place.push_back(box);
	int oldTestIndex = t.testIndex;
    vector<string> new_qrs;
    short_words_cache.clear();

    string aux_word;
	if (t.testIndex >= 0) {
        box_state result = g_tests.evaluateBox(t.testIndex, box, aux_word, new_qrs, short_words_cache);
        if (result != open && result != open_with_qr) {
            t.aux_word.assign(aux_word);
            t.testResult = result;
            fprintf(stderr, "Eliminated %s with test %s with result %d\n", box.name.c_str(), g_tests.getName(t.testIndex), result);
            return true;
        } else if (result == open_with_qr) {
            fprintf(stderr,"Retested %d, new qrs len %d\n", result, new_qrs.size());
            for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
                fprintf(stderr, "New QR is %s\n", (*it).c_str());
                box.qr.getName(*it); // Also adds qr to the box's list
            }
            t.qr_desc = box.qr.min_pow_desc();
      } else { 
            fprintf(stderr, "FAILED to eliminate %s with test %s with result %d\n", box.name.c_str(), g_tests.getName(t.testIndex), result);
        }
    }

	if (t.testIndex == -2 && !g_options.fillHoles) {
	    return true;
    }

    // Check if the box is now small enough that some former qrs actually kill it
    Params<ACJ> cover = box.cover();
    vector<string> quasiRelators = box.qr.wordClasses();
    for (vector<string>::iterator it = quasiRelators.begin(); it != quasiRelators.end(); ++it) {
        // So not idenity and absUB(w.b) < 1
        SL2ACJ w = g_tests.construct_word(*it, cover, short_words_cache); 
        if (not_identity(w)) {
//            fprintf(stderr, "Failed qr %s at %s\n", (*it).c_str(), box.name.c_str());
//            fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
//                            absLB(w.b), absLB(w.c), absLB(w.a - 1.), absLB(w.d - 1.), absLB(w.a + 1.), absLB(w.d + 1.));
            t.aux_word.assign(*it);
            t.testResult = killed_failed_qr;
            return true;
        }
    }

	if (g_options.improveTree || !t.lChild) {
		for (int i = 0; i < g_tests.size(); ++i) {
			vector<box_state>& th = history[i];
			while (th.size() <= depth && (th.size() < depth-6 || th.empty() || th.back() == open)) {
				box_state result = g_tests.evaluateCenter(i, place[th.size()]);
				th.push_back(result);
			}
			if (th.back() != open) {
                new_qrs.clear();
				box_state result = g_tests.evaluateBox(i, box, aux_word, new_qrs, para_cache, short_words_cache);

                switch (result) {
                    case variety_nbd : 
                    case killed_no_parabolics :
                    case killed_bad_parabolic :
                    case killed_failed_qr :
                    case killed_identity_impossible :
                    case killed_elliptic : {
                        t.aux_word.assign(aux_word);   
                    }
                    case killed_bounds :
                    case killed_parabolics_impossible : {
                        t.testIndex = i;
                        t.testResult = result;
                        return true;
                    }
                    case open_with_qr : {
//                        fprintf(stderr,"Result %d, new qrs len %d\n", result, new_qrs.size());
                        for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
//                            fprintf(stderr, "New QR is %s\n", (*it).c_str());
                            box.qr.getName(*it); // Also adds qr to the box's list
                        }
                        t.qr_desc = box.qr.min_pow_desc();
                        break;
                    }
                }
            }
		}
	}

	if (g_options.wordSearchDepth >= 0 && (g_options.improveTree || !t.lChild)) {
		while (depth - searchedDepth > g_options.wordSearchDepth) {
			Box& searchPlace = place[++searchedDepth];
			vector<string> searchWords = findWords( searchPlace.center(), vector<string>(), -200, g_options.maxWordLength, box.qr.wordClasses());
            string new_word = searchWords.back();

            int old_size = g_tests.size();
    		int new_index = g_tests.add(new_word);
			history.resize(g_tests.size());

            if (old_size < g_tests.size()) {
                fprintf(stderr, "search (%s) found %s(%s)\n",
                    searchPlace.qr.desc().c_str(), new_word.c_str(), searchPlace.name.c_str());

                new_qrs.clear();
                box_state result = g_tests.evaluateBox(new_index, box, aux_word, new_qrs, para_cache, short_words_cache);

                switch (result) {
                    case variety_nbd : 
                    case killed_no_parabolics :
                    case killed_bad_parabolic :
                    case killed_failed_qr :
                    case killed_identity_impossible :
                    case killed_elliptic : {
                        t.aux_word.assign(aux_word);   
                    }
                    case killed_bounds :
                    case killed_parabolics_impossible : {
                        t.testIndex = new_index;
                        t.testResult = result;
                        return true;
                    }
                    case open_with_qr : {
    //                        fprintf(stderr,"Result %d, new qrs len %d\n", result, new_qrs.size());
                        for (vector<string>::iterator it = new_qrs.begin(); it != new_qrs.end(); ++it) {
    //                            fprintf(stderr, "New QR is %s\n", (*it).c_str());
                            box.qr.getName(*it); // Also adds qr to the box's list
                        }
                        t.qr_desc = box.qr.min_pow_desc();
                        break;
                    }
                }
            }
		}
	}

	t.testIndex = -1;

	if (!t.lChild) {
		if (depth >= g_options.maxDepth || ++g_boxesVisited >= g_options.maxSize || ++newDepth > g_options.inventDepth) {
//            fprintf(stderr,"Deph %d, max depth %d, boxes_visited %d, max size %d, newDepth %d, invent depth %d\n", depth, g_options.maxDepth, g_boxesVisited, g_options.maxSize, newDepth, g_options.inventDepth);
	        Params<XComplex> params = box.center();
	        Params<XComplex> nearer = box.nearer();
            double area = areaLB(nearer);
            fprintf(stderr, "HOLE %s has min area: %f center lat: %f + I %f lox: %f + I %f par: %f + I %f (%s)\n", box.name.c_str(), area, params.lattice.re, params.lattice.im, params.loxodromic_sqrt.re, params.loxodromic_sqrt.im, params.parabolic.re,params.parabolic.im, box.qr.desc().c_str());
            return false;
		}
		t.lChild = new PartialTree();
		t.rChild = new PartialTree();
	}

	bool isComplete = true;

	isComplete = refineRecursive(box.child(0), *t.lChild, depth+1, history, place, newDepth, searchedDepth) && isComplete;
	if (place.size() > depth+1)
		place.resize(depth+1);
	for (int i = 0; i < g_tests.size(); ++i) {
		if (history[i].size() > depth)
			history[i].resize(depth);
	}
	if (searchedDepth > depth)
		searchedDepth = depth;
	if (isComplete || depth < g_options.truncateDepth)
		isComplete = refineRecursive(box.child(1), *t.rChild, depth+1, history, place, newDepth, searchedDepth) && isComplete;
	if (oldTestIndex >= 0 && t.testIndex != oldTestIndex) {
		fprintf(stderr, "invalid box %s(%s) %d %s\n", g_tests.getName(oldTestIndex), box.name.c_str(),
			treeSize(t), isComplete ? "Patched" : "Unpatched");
	}
	if (!isComplete && depth >= g_options.truncateDepth) {
		truncateTree(t);
	}
	return isComplete;
}

void refineTree(Box box, PartialTree& t)
{
	TestHistory history(g_tests.size());
	vector<Box> place;
	int searchedDepth = 0;
	refineRecursive(box, t, 0, history, place, 0, searchedDepth);
}

void printTree(PartialTree& t)
{
    char type = 'F';
//    printf("%d\n", t.testResult);
    switch (t.testResult) {
        case open :
        case open_with_qr : {
            if (t.lChild && t.rChild) {
                printf("X\n");
                printTree(*t.lChild);
                printTree(*t.rChild);
            } else {
		        printf("HOLE (%s)\n", t.qr_desc.c_str());
            }
            return;
        }
        case killed_bounds : {
		    printf("%s\n", g_tests.getName(t.testIndex));
            return;
        }
        case killed_no_parabolics : type = 'K'; break;
        case variety_nbd : type = 'V'; break;
        case killed_parabolics_impossible : type = 'P'; break;
        case killed_identity_impossible : type = 'I'; break;
        case killed_failed_qr : type = 'Q'; break;
        case killed_bad_parabolic : type = 'L'; break;
        case killed_elliptic : type = 'E'; break;
    }
//    printf("%c\n", type);
    string killer;
    if (type == 'P') {
        killer = g_tests.getName(t.testIndex);
    } else if (type == 'E') {
        killer = g_tests.getName(t.testIndex);
        killer += "," + t.aux_word; 
    } else {
        killer = t.aux_word;
    }
    printf("%c(%s)\n", type, killer.c_str());
}

const char* g_programName;


static struct option longOptions[] = {
	{"box",	required_argument, NULL, 'b' },
	{"words", required_argument, NULL, 'w' },
	{"powers", required_argument, NULL, 'p'},
	{"maxDepth", required_argument, NULL, 'd' },
	{"inventDepth", required_argument, NULL, 'i' },
	{"improveTree", no_argument, NULL, 'I'},
	{"truncateDepth", required_argument, NULL, 't' },
	{"maxSize", required_argument, NULL, 's' },
	{"wordSearchDepth", required_argument, NULL, 'B'},
	{"fillHoles", no_argument, NULL, 'f'},
	{"maxArea", required_argument, NULL, 'a'},
	{NULL, 0, NULL, 0}
};

static char optStr[1000] = "";
void setOptStr() {
	char* osp = optStr;
	for (int i = 0; longOptions[i].name; ++i) {
		*osp++ = longOptions[i].val;
		if (longOptions[i].has_arg != no_argument)
			*osp++ = ':';
	}
	*osp = '\0';
}

void usage()
{
	const char* longUsage = "\
       --box <box_id>\n\
       		Box ID for the root of the input and output trees.\n\
		May include characters not in [01], which are ignored.\n\
\n\
Options controlling which relators to use:\n\
	[ --words	<words_file> ]\n\
		File containing starting list of words to try.\n\
	[ --wordSearchDepth <n> ]\n\
		Perform a search for relators when visiting a node at least n levels deep.\n\
\n\
Options controlling which relators eliminate boxes:\n\
	[ --powers <powers_file> ]\n\
		File containing impossible relator definitions.\n\
		See ImpossibleRelators::load(...)\n\
\n\
Options controlling tree manipulation:\n\
	[ --maxDepth <n> ]\n\
		Don't descend more than n levels deeper than the root box.\n\
	[ --inventDepth <n> ]\n\
		Don't descend more than n levels deeper than the terminal node of the input tree.\n\
	[ --maxSize <n> ]\n\
		Don't allow the output tree to have more than n nodes.\n\
	[ --truncateDepth <n> ]\n\
		Don't emit holes more than n levels deeper than the root node.\n\
		Instead, replace the subtree-with-holes with a single hole.\n\
	[ --improveTree ]\n\
		If set, attempt to directly eliminate internal nodes of the input tree.\n\
	[ --fillHoles ]\n\
		If set, attempt to patch holes in the input tree.\n\
";
	fprintf(stderr, "Usage: %s %s\n\n%s", g_programName, optStr, longUsage);
}

void loadWords(set<string>& s, const char* fileName)
{
	FILE* fp = fopen(fileName, "r");
	char buf[10000];
	while (fp && fgets(buf, sizeof(buf), fp)) {
		int n = -1 + strlen(buf);
		if (buf[n] == '\n')
			buf[n] = '\0';
		s.insert(buf);
	}
}

int main(int argc, char** argv)
{
	setOptStr();
	if (argc < 2) {
		usage();
		exit(1);
	}

	int ch;
	while ((ch = getopt_long(argc, argv, optStr, longOptions, NULL)) != -1) {
        fprintf(stderr,"Arg %c, %s\n", ch, optarg);
		switch(ch) {
		case 'b': g_options.boxName = optarg; break;
		case 'w': g_options.wordsFile = optarg; break;
		case 'p': g_options.powersFile = optarg; break;
		case 'd': g_options.maxDepth = atoi(optarg); break;
		case 'i': g_options.inventDepth = atoi(optarg); break;
		case 'I': g_options.improveTree = atoi(optarg); break;
		case 't': g_options.truncateDepth = atoi(optarg); break;
		case 's': g_options.maxSize = atoi(optarg); break;
		case 'B': g_options.wordSearchDepth = atoi(optarg); break;
		case 'f': g_options.fillHoles = atoi(optarg); break;
		case 'm': g_margulis_bound = atof(optarg); break;
		}
	}

	Box box;
	for (const char* boxP = g_options.boxName; *boxP; ++boxP) {
		if (*boxP == '0') {
			box = box.child(0);
		} else if (*boxP == '1') {
			box = box.child(1);
		}
	}
	
	g_tests.load(g_options.wordsFile);
	g_tests.loadImpossibleRelations(g_options.powersFile);
	
	PartialTree t = readTree();
	refineTree(box, t);
	printTree(t);
	fprintf(stderr, "%d nodes added\n", g_boxesVisited);
}
