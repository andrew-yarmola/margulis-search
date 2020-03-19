#ifndef __Box_h
#define __Box_h
#include "types.hh"
#include "SL2.hh"
#include "AJ.h"
#include "QuasiRelators.h"

struct Box {
    Box();
    std::string name;
    std::string desc();
    QuasiRelators qr;
    Box child(int dir) const;
    Params<Complex> center() const { return _center; }
    Params<AJ> cover() const { return _cover; }
//    Params<Complex> nearer() const { return _nearer; } // returns all values closer to 0 than in box or 0 if box overlaps
//    Params<Complex> further() const { return _further; } // returns all values futher from 0 that in the box
//    Params<Complex> greater() const { return _greater; } // returns all values greater than in the box
    SL2<Complex> x_center() const { return _x_center; }
    SL2<Complex> y_center() const { return _y_center; }
    SL2<AJ> x_cover() const { return _x_cover; }
    SL2<AJ> y_cover() const { return _y_cover; }
private:
    int pos;
    double center_digits[6];
    double size_digits[6];
    double box_center[6];
    double box_size[6];
    Params<Complex> _center;
    Params<AJ> _cover;
//    Params<Complex> _nearer;
//    Params<Complex> _further;
//    Params<Complex> _greater;
    void compute_center_and_size();
    void compute_cover();
//    void compute_nearer();
//    void compute_further();
//    void compute_greater();
    SL2<Complex> _x_center;
    SL2<Complex> _y_center;
    SL2<AJ> _x_cover;
    SL2<AJ> _y_cover;
};

Box get_box(std::string code);

#endif // __Box_h
