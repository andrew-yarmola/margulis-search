#include "SL2C.h"
#include "SL2ACJ.h"
#include "Params.h"

/*
* Let x and y realize the margulis constant such that re(length(x)) =< re(length(y))
* We will conjugate x to have axis(x) = {0,\infty}
*/

const SL2C construct_x(const Params<XComplex>& params)
{
    // We want the matrix {{ exp(L/2), 0 }, { 0, exp(-L/2) }}
    XComplex sl2 = params.sinhL2;
    XComplex cl2 = sqrt(((sl2 * sl2).z + 1.).z).z;
	return SL2C((cl2 + sl2).z, 0.,
                 0., (cl2 - sl2).z);
}

const SL2C construct_y(const Params<XComplex>& params)
{
    // We want the matrix {{ cosh(D/2) + cosh(P)sinh(D/2), -sinh(D/2)sinh(P) }, 
    //                     { sinh(D/2)sinh(P), cosh(D/2) - cosh(P)sinh(D/2) }}
    XComplex sd2 = params.sinhD2; 
    XComplex sd2 = sqrt(((sd2 * sd2).z + 1.).z).z;
    XComplex sp = params.sinhP;
    XComplex cp = sqrt(((sp * sp).z + 1.).z).z;
	return SL2C((cd2 + (cp * sd2).z).z, -(sd2 * sp).z,
                (sd2 * sp).z, (cd2 - (cp * sd2).z).z);
}

const SL2ACJ construct_x(const Params<ACJ>& params)
{
    ACJ zero = ACJ(0.);
    ACJ sl2 = params.sinhL2;
    ACJ cl2 = sqrt(sl2 * sl2 + 1.);
	return SL2ACJ(cl2 + sl2, zero,
                  zero, cl2 - sl2);
}

const SL2ACJ construct_y(const Params<ACJ>& params)
{
    // We want the matrix {{ cosh(D/2) + cosh(P)sinh(D/2), -sinh(D/2)sinh(P) }, 
    //                     { sinh(D/2)sinh(P), cosh(D/2) - cosh(P)sinh(D/2) }}
    ACJ sd2 = params.sinhD2; 
    ACJ cd2 = sqrt(sd2 * sd2 + 1.);
    ACJ sp = params.sinhP;
    ACJ cp = sqrt(sp * sp + 1.);
	return SL2ACJ(cd2 + cp * sd2, -(sd2 * sp),
                sd2 * sp, cd2 - cp * sd2);
}

int x_power(std::string w) {
    int count = 0;
    for (std::string::size_type p = 0; p < w.size(); ++p) {
        if (w[p] == 'x' || w[p] == 'X') ++count;
    }
    return count;
} 

int y_power(std::string w) {
    int count = 0;
    for (std::string::size_type p = 0; p < w.size(); ++p) {
        if (w[p] == 'y' || w[p] == 'Y') ++count;
    }
    return count;
} 

bool x_power_sort(std::string a, std::string b) { return x_power(a) < x_power(b); }

bool y_power_sort(std::string a, std::string b) { return y_power(a) < y_power(b); }

