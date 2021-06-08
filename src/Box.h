#ifndef __Box_h
#define __Box_h
#include "types.hh"
#include "SL2.hh"
#include "AJCC.h"
#include "QuasiRelators.h"

#define DIM 4
#define SCL 2   
// Initial box dimensions are therefore
//  2 * (2, 2^(3/4), 2^(2/4), 2^(1/4)). The last is > 2.37

struct Box {
    Box();
    std::string name;
    std::string desc();
    QuasiRelators qr;
    Box child(int dir) const;
    Params<Complex> center() const { return _center; }
    Params<AJCC> cover() const { return _cover; }
    SL2<Complex> x_center() const { return _x_center; }
    SL2<Complex> y_center() const { return _y_center; }
    SL2<AJCC> x_cover() const { return _x_cover; }
    SL2<AJCC> y_cover() const { return _y_cover; }
private:
    int pos;
    double center_digits[DIM];
    double size_digits[DIM];
    double box_center[DIM];
    double box_size[DIM];
    Params<Complex> _center;
    Params<AJCC> _cover;
    void compute_center_and_size();
    void compute_cover();
    SL2<Complex> _x_center;
    SL2<Complex> _y_center;
    SL2<AJCC> _x_cover;
    SL2<AJCC> _y_cover;
};

Box get_box(std::string code);

#endif // __Box_h
