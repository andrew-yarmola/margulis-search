#include <stdio.h>
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJ.h"
#include "types.hh"
#include "TubeSearch.hh"
#include "TestCollection.hh"

#define MAX_DEPTH 200

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

  printf("Box: %s", box.desc().c_str());

  SL2<AJ> x = box.x_cover(); 
  SL2<AJ> y = box.y_cover();
//  printf("x is :\n");
//  print_SL2(x);
//  printf("y is :\n");
//  print_SL2(y);
//  SL2<Complex> c_x = construct_x(box.center());
//  SL2<Complex> c_y = construct_y(box.center());

//  float_pair up_up = four_cosh_margulis(x,y,true,true);
//  float_pair up_lb = four_cosh_margulis(x,y,true,false);
//  float_pair lb_up = four_cosh_margulis(x,y,false,true);
//  float_pair lb_lb = four_cosh_margulis(x,y,false,false);

//  printf("4 cosh(margulis) between %f (%f) and %f (%f)\n", lb_lb.first, lb_up.first, up_up.first, up_lb.first);
//  printf("exp(2t) between %f (%f) and %f (%f)\n", lb_lb.second, up_lb.second, up_up.second, lb_up.second);

}
