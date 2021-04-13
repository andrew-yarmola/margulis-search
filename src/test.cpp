#include <stdio.h>
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJ.h"
#include "types.hh"
#include "TubeSearch.hh"
#include "TestCollection.hh"

#define MAX_DEPTH 200

bool g_debug=true;

using namespace std;

int main(int argc,char**argv)
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s position < data\n", argv[0]);
        exit(1);
    }
    char where[MAX_DEPTH];
    size_t depth = 0;
    while (argv[1][depth] != '\0') {
        if (argv[1][depth] != '0' && argv[1][depth] != '1'){
            fprintf(stderr,"bad position %s\n",argv[1]);
            exit(2);
        }
        where[depth] = argv[1][depth];
        depth++;
    }
    where[depth] = '\0';

	Box box;
	for (char* dir = (char *) &where; *dir; ++dir) {
		if (*dir == '0') {
			box = box.child(0);
		} else if (*dir == '1') {
			box = box.child(1);
		}
	}

//  printf("Syllables %d\n", syllables("xxXyyYx"));
/*  printf("x strip xxXyxyYx : %s\n", x_strip("xxXyxyYx").c_str());
  printf("y strip yYyxxXyyYxyyY : %s\n", y_strip("yYyxxXyyYxyyY").c_str());
  printf("y strip yYyxxXyyYx : %s\n", y_strip("yYyxxXyyYx").c_str());
  printf("x rstrip xxXyxyYx : %s\n", x_rstrip("xxXyxyYx").c_str());
  printf("y rstrip yYyxxXyyYxyyY : %s\n", y_rstrip("yYyxxXyyYxyyY").c_str());
  printf("y rstrip yYyxxXyyYx : %s\n", y_rstrip("yYyxxXyyYx").c_str());
*/

  printf("Box: %s", box.desc().c_str());
  SL2<AJ> w = construct_word("yxyXXyx", box.cover());
  string proven = proven_identity("yxyXXyx", box.cover());
  print_SL2(w);
  AJ rl = four_cosh_re_length(w);
  //printf("4 cosh(re(L)) between %f and %f\n", absLB(rl), absUB(rl));
  print_type(rl);
  print_type(w.a + w.d);
  print_SL2(w*w);
  printf("proven %s\n", proven.c_str());
  

//  SL2<AJ> x = box.x_cover(); 
//  SL2<AJ> y = box.y_cover();
//  printf("x is :\n");
//  print_SL2(x);
//  printf("y is :\n");
//  print_SL2(y);
//
//  pair<AJ, AJ> marg_pair = four_cosh_margulis_simple(x,y);  
//
//  printf("exp(2t) between %f and %f\n", absLB(marg_pair.second), absUB(marg_pair.second));  

}
