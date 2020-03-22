#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "rapidcsv.h"
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJ.h"
#include "roundoff.h"
#include "types.hh"
#include "TubeSearch.hh"
#include "TestCollection.hh"

#define ERR 0.000001
#define CERR 0.00001
#define DERR 0.00015
#define AJERR 0.000835

using namespace std;

const XComplex to_XComplex(const Complex& z) {
    return XComplex(z.real(), z.imag());
}

const AJ to_AJ(const Complex& z) {
    return AJ(to_XComplex(z));
}

int main(int argc,char**argv)
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s census_csv_file\n", argv[0]);
        exit(1);
    }

    initialize_roundoff();

    rapidcsv::Document doc(argv[1], rapidcsv::LabelParams(0,-1));
    const size_t row_count = doc.GetRowCount();
     for (size_t i = 0; i < row_count; ++i) {
        string name = doc.GetCell<string>("name", i);
        string codes = doc.GetCell<string>("box codes", i);
        double sinhdx = doc.GetCell<double>("sinhdx", i);
        double sinhdy = doc.GetCell<double>("sinhdy", i);
        double coshmu = doc.GetCell<double>("coshmu", i);
        double cosf = doc.GetCell<double>("cosf", i);
        double sintx2 = doc.GetCell<double>("sintx2", i);
        double sinty2 = doc.GetCell<double>("sinty2", i);

        Complex Lx = parse_complex(doc.GetCell<string>("x", i)); 
        Complex Ly = parse_complex(doc.GetCell<string>("y", i)); 
        Lx = shift_imag_around_zero(Lx);
        Ly = shift_imag_around_zero(Ly);

        vector<string> box_codes;
        split_string(codes, "\"[',] ", box_codes);

        for (string box_code : box_codes) {
          Box box = get_box(box_code);
          Params<Complex> center = box.center();
          Params<AJ> cover = box.cover();
          /* printf("%s %f %f %f %f %f %f %s\n", name.c_str(), sinhdx, sinhdy, coshmu, cosf, sintx2, sinty2, box_code.c_str());
          printf("Box: %s", box.desc().c_str());
          printf("Manifold %s\n", name.c_str());
          printf("Lx: %f + i %f\n", Lx.real(), Lx.imag());
          Complex sinhLx2 = sinh(Lx/2);
          printf("sinhLx2: %f + i %f\n", sinhLx2.real(), sinhLx2.imag());
          printf("center sinhLx2: %f + i %f\n", center.sinhLx2.real(), center.sinhLx2.imag());
          */

          // printf("%s\n", name.c_str());
          // Center testing
          assert(absUB(center.sinhdx - sinhdx) < ERR);
          assert(absUB(center.sinhdy - sinhdy) < ERR);
          assert(absUB(center.coshmu - coshmu) < ERR);
          assert(absUB(center.cosf - cosf) < ERR);
          assert(absUB(center.sintx2 - sintx2) < ERR);
          assert(absUB(center.sinty2 - sinty2) < ERR);
          assert(absUB(center.coshLx2 - cosh(Lx/2)) < CERR);
          assert(absUB(center.coshLy2 - cosh(Ly/2)) < CERR);

          // Test margulis compuations
          pair<Complex, Complex> center_pair = four_cosh_margulis_simple(box.x_center(), box.y_center());
          // Complex four_cosh_mu = center.coshmu * 4;
          // printf("center 4coshmu: %f + i %f\n", four_cosh_mu.real(), four_cosh_mu.imag());
          // printf("Box: %s", box.desc().c_str());
          assert(absUB(center_pair.first - center.coshmu * 4) < CERR);
          assert(absUB(center_pair.second - center.expdx * center.expdx) < CERR);

          // Test construct math
          SL2<Complex> x_center = box.x_center();          
          SL2<Complex> y_center = box.y_center();          
          SL2<Complex> xy_center = construct_word("xy", center); 
          SL2<Complex> xyx_center = construct_word("xyx", center); 
          SL2<Complex> yxxYXy_center = construct_word("yxxYXy", center); 
          assert(absUB(dist(xy_center, x_center * y_center)) < ERR);          
          assert(absUB(dist(xyx_center, x_center * y_center * x_center)) < ERR);          
          assert(absUB(dist(yxxYXy_center, y_center * x_center * x_center * inverse(y_center) * inverse(x_center) * y_center)) < ERR);          

          // Test Jorgensen
          Complex j_ww_xy_center = jorgensen(x_center, y_center);
          Complex j_ww_yx_center = jorgensen(y_center, x_center);
          Complex j_xw_xy_center = jorgensen_xw(y_center, center);
          Complex j_wx_yx_center = jorgensen_wx(y_center, center);
          Complex j_wy_xy_center = jorgensen_wy(x_center, center);
          Complex j_yw_yx_center = jorgensen_yw(x_center, center);
          Complex j_xy_center = jorgensen_xy(center);
          Complex j_yx_center = jorgensen_yx(center);
          assert(absUB(j_ww_xy_center - j_xy_center) < ERR); 
          assert(absUB(j_xw_xy_center - j_xy_center) < ERR); 
          assert(absUB(j_wy_xy_center - j_xy_center) < ERR); 
          assert(absUB(j_ww_yx_center - j_yx_center) < ERR); 
          assert(absUB(j_wx_yx_center - j_yx_center) < ERR); 
          assert(absUB(j_yw_yx_center - j_yx_center) < ERR); 
          assert(absLB(j_xy_center) >= 1.0); 
          assert(absLB(j_yx_center) >= 1.0); 

          // Cover testing
          assert(absUB(cover.sinhdx - sinhdx) < ERR);
          assert(absUB(cover.sinhdy - sinhdy) < ERR);
          assert(absUB(cover.coshmu - coshmu) < ERR);
          assert(absUB(cover.cosf - cosf) < ERR);
          assert(absUB(cover.sintx2 - sintx2) < ERR);
          assert(absUB(cover.sinty2 - sinty2) < ERR);
          assert(absUB(cover.coshLx2 - to_AJ(cosh(Lx/2))) < AJERR);
          assert(absUB(cover.coshLy2 - to_AJ(cosh(Ly/2))) < AJERR);

          // Test margulis compuations
          pair<AJ, AJ> cover_pair = four_cosh_margulis_simple(box.x_cover(), box.y_cover());
          // print_type("cover margulis exact:", cover.coshmu * 4);
          // print_type("computed margulis:", cover_pair.first);
          // print_type("computed diff:", cover_pair.first - cover.coshmu * 4);
          assert(absUB(cover_pair.first - cover.coshmu * 4) < AJERR);
          assert(absUB(cover_pair.second - cover.expdx * cover.expdx) < AJERR);

          // Test construct math
          SL2<AJ> x_cover = box.x_cover();          
          SL2<AJ> y_cover = box.y_cover();          
          SL2<AJ> xy_cover = construct_word("xy", cover); 
          SL2<AJ> xyx_cover = construct_word("xyx", cover); 
          SL2<AJ> yxxYXy_cover = construct_word("yxxYXy", cover); 
          assert(absUB(dist(xy_cover, x_cover * y_cover)) < ERR);          
          assert(absUB(dist(xyx_cover, x_cover * y_cover * x_cover)) < ERR);          
          assert(absUB(dist(yxxYXy_cover, y_cover * x_cover * x_cover * inverse(y_cover) * inverse(x_cover) * y_cover)) < ERR);          

          // Test Jorgensen
          AJ j_ww_xy_cover = jorgensen(x_cover, y_cover);
          AJ j_ww_yx_cover = jorgensen(y_cover, x_cover);
          AJ j_xw_xy_cover = jorgensen_xw(y_cover, cover);
          AJ j_wx_yx_cover = jorgensen_wx(y_cover, cover);
          AJ j_wy_xy_cover = jorgensen_wy(x_cover, cover);
          AJ j_yw_yx_cover = jorgensen_yw(x_cover, cover);
          AJ j_xy_cover = jorgensen_xy(cover);
          AJ j_yx_cover = jorgensen_yx(cover);
          // print_type("j_ww_xy_cover:", j_ww_yx_cover);
          // print_type("j_xy_cover:", j_yx_cover);
          // print_type("diffr:", j_ww_yx_cover - j_yx_cover);
          assert(absUB(j_ww_xy_cover - j_xy_cover) < DERR); 
          assert(absUB(j_xw_xy_cover - j_xy_cover) < DERR); 
          assert(absUB(j_wy_xy_cover - j_xy_cover) < DERR); 
          assert(absUB(j_ww_yx_cover - j_yx_cover) < DERR); 
          assert(absUB(j_wx_yx_cover - j_yx_cover) < DERR); 
          assert(absUB(j_yw_yx_cover - j_yx_cover) < DERR); 
          assert(absLB(j_xy_cover) >= 1.0); 
          assert(absLB(j_yx_cover) >= 1.0); 
  
        }
    }

    if(!roundoff_ok()){
        printf("Underflow may have occurred\n");
        exit(1);
    }

    printf("PASSED\n");

}

/*
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
*/