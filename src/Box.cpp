#include "Box.h"

// Initial box dimensions are 2^(18/6), 2^(17/6), ..., 2^(13/6). The last is > 4.49

double scale[6];
static bool scale_initialized = false; 
Box::Box() {
	if (!scale_initialized) {
		scale_initialized = true;
		for (int i = 0; i < 6; ++i) {
			scale[i] = pow(2, -i / 6.0);
		}
	}
	for (int i = 0; i < 6; ++i) {
		center_digits[i] = 0;
		size_digits[i] = 8;
	}
	pos = 0;
  compute_center_and_size();
	compute_cover();
  compute_nearer();
	compute_further();
	compute_greater();
}

Box Box::child(int dir) const
{
	Box child(*this);
	child.size_digits[pos] *= 0.5;
	child.center_digits[pos] += (2*dir-1)*child.size_digits[pos];
	++child.pos;
	if (child.pos == 6) { child.pos = 0; }

	child.name = name;
	child.name.append(1, '0'+dir);

	child.qr = qr;

  child.compute_center_and_size();
	child.compute_cover();
  child.compute_nearer();
	child.compute_further();
	child.compute_greater();
	return child;
}

void Box::compute_center_and_size()
{
	for (int i = 0; i < 6; ++i) {
        // GMT paper page 419 of Annals
        // box_size guarantees that :
        // box_center - box_size <= true_center - true_size
        // box_center + box_size >= true_center + true_size
        // where box operations are floating point. 
        box_center[i] = scale[i]*center_digits[i];
        box_size[i]= (1+2*EPS)*(size_digits[i]*scale[i]+HALFEPS*fabs(center_digits[i]));
    }
    _center.sinhP = XComplex(box_center[3], box_center[0]);
    _center.sinhD2 = XComplex(box_center[4], box_center[1]);
    _center.sinhL2 = XComplex(box_center[5], box_center[2]);
    _center.coshP = sqrt(((_center.sinhP * _center.sinhP).z + 1).z).z;
    _center.coshD2 = sqrt(((_center.sinhD2 * _center.sinhD2).z + 1).z).z;
    _center.coshL2 = sqrt(((_center.sinhL2 * _center.sinhL2).z + 1).z).z;
}

void Box::compute_cover()
{
	_cover.sinhP = ACJ(
		XComplex(box_center[3], box_center[0]),
		XComplex(box_size[3], box_size[0]),
		0.,
		0.
	);
	_cover.sinhD2 = ACJ(
		XComplex(box_center[4], box_center[1]),
		0.,
		XComplex(box_size[4], box_size[1]),
		0.
	);
	_cover.sinhL2 = ACJ(
		XComplex(box_center[5], box_center[2]),
		0.,
		0.,
		XComplex(box_size[5], box_size[2])
	);
  _cover.coshP  = sqrt( _cover.sinhP  * _cover.sinhP  + 1);
  _cover.coshD2 = sqrt( _cover.sinhD2 * _cover.sinhD2 + 1);
  _cover.coshL2 = sqrt( _cover.sinhL2 * _cover.sinhL2 + 1);
}

void Box::compute_nearer()
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
        if (center_digits[i] > 0 && // center is positive 
            center_digits[i] > size_digits[i] &&  // true diff is positive
            box_center[i]    > box_size[i]) { // machine diff is >= 0
            // Want lower bound on true_center - true_size.  Assume no overflow or underflow 
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      box_center - box_size <= true_center - true_size
            // Now, in machine arthimetric, by IEEE, if 
            //      box_center > box_size then box_center (-) box_size >= 0.
            // Lemma 7 gives,
            //      (1-EPS)(*)( box_center (-) box_size ) <= box_center - box_size <= true_center - box_size. 
            m[i] = (1-EPS)*(box_center[i] - box_size[i]);
        } else if (center_digits[i] < 0 && // center is negative
                   center_digits[i] < -size_digits[i] && // true sum is negative
                   box_center[i]    < -box_size[i]) {  // machine sum is negative
            // Want upper bound on true_center - true_size.  Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size
            // Now, in machine arthimetric, by IEEE, if 
            //      -box_center > box_size then (-box_center) (-) box_size >= 0.
            // Lemma 7 gives,
            //      (1-EPS)(*)( (-box_center) (-) box_size ) <= -box_center - box_size <= -true_center - true_size.
            // So,
            //      -((1-EPS)(*)( (-box_center) (-) box_size )) >= true_center + true_size.
            // Note, negation is exact for machine numbers
            m[i] = -((1-EPS)*((-box_center[i]) - box_size[i]));
        }
	}
	
	_nearer.sinhP  = XComplex(m[3], m[0]);
	_nearer.sinhD2 = XComplex(m[4], m[1]);
	_nearer.sinhL2 = XComplex(m[5], m[2]);
}

void Box::compute_further()
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
		if (center_digits[i] > -size_digits[i]) { // true sum is positive 
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size
            // By IEEE (+) and (-) resepct <= and >=, so box_center (+) box_size >=0 and
            // Lemma 7 for floating point arithmetic gives and upper bound
            //      (1+EPS)(*)(box_center (+) box_size) >= box_center + box_size >= true_center + true_size
		    m[i] = (1+EPS)*(box_center[i] + box_size[i]);
        } else { // true sum is <= 0
            // Want lower bound of true_center - true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      box_center - box_size <= true_center - true_size
            // By IEEE, (+) and (-) respects <= and >=, and negation is exact.
            // Thus, (-box_center) (+) box_size >=0 and Lemma 7 for floating point arithmetic gives
            //        (1+EPS)(*)( (-box_center) (+) box_size) ) >= (-box_center) + box_size
            // So,
            //      -((1+EPS)(*)( (-box_center) (+) box_size) ))<= box_center - box_size <= true_center - true_size
            m[i] = -((1+EPS)*((-box_center[i]) + box_size[i]));
        }
	}
	
	_further.sinhP  = XComplex(m[3], m[0]);
	_further.sinhD2 = XComplex(m[4], m[1]);
	_further.sinhL2 = XComplex(m[5], m[2]);
}

void Box::compute_greater()
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
		if (center_digits[i] > -size_digits[i]) { // true sum is positive
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size.
            // Notice that box_center + box_size >= true_center + true_size > 0.
            // By IEEE, box_center (+) box_size >=0, as it's guanrateed to evaluate to nearest representable.
            // Lemma 7 for floating point arithmetic gives and upper bound
            //      (1+EPS)(*)(box_center (+) box_size) >= box_center + box_size >= true_center + true_size
		    m[i] = (1+EPS)*(box_center[i] + box_size[i]);
        } else if (center_digits[i] < -size_digits[i] && // true sum is negative
                   box_center[i]    < -box_size[i]) { // machine sum is <= 0
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size.
            // Notice that box_center + box_size < 0.
            // By IEEE, box_center (+) box_size <= 0, as it's guanrateed to evaluate to nearest representable.
            // Lemma 7 for floating point arithmetic gives a bound
            //      (1-EPS)(*)| box_center (+) box_size | < | box_center + box_size |
            // So,
            //      -((1-EPS)(*)(-(box_center (+) box_size))) >= box_center + box_size >= true_center + true_size
            m[i] = -((1-EPS)*(-(box_center[i] + box_size[i])));
        }
	}
	
	_greater.sinhP  = XComplex(m[3], m[0]);
	_greater.sinhD2 = XComplex(m[4], m[1]);
	_greater.sinhL2 = XComplex(m[5], m[2]);
}

//Params<XComplex> Box::offset(const double* offset) const
//{
//	Params<XComplex> result;
//	result.sinhP = XComplex(
//		scale[3]*(offset[3]*size_digits[3] + center_digits[3]),
//		scale[0]*(offset[0]*size_digits[0] + center_digits[0])
//	);
//	result.sinhD2 = XComplex(
//		scale[4]*(offset[4]*size_digits[4] + center_digits[4]),
//		scale[1]*(offset[1]*size_digits[1] + center_digits[1])
//	);
//	result.sinhL2 = XComplex(
//		scale[5]*(offset[5]*size_digits[5] + center_digits[5]),
//		scale[2]*(offset[2]*size_digits[2] + center_digits[2])
//	);
//	return result;
//}
