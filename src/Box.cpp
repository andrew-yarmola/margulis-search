#include "Box.h"
#include "Params.hh"

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

std::string Box::desc() {
  ACJ sinhP = _cover.sinhP;
  ACJ sinhD2 = _cover.sinhD2;
  ACJ sinhL2 = _cover.sinhL2; 
  ACJ coshP = _cover.coshP;
  ACJ coshD2 = _cover.coshD2;
  ACJ coshL2 = _cover.coshL2; 
  Complex c_sinhP = _center.sinhP;
  Complex c_sinhD2 = _center.sinhD2;
  Complex c_sinhL2 = _center.sinhL2; 
  Complex c_coshP = _center.coshP;
  Complex c_coshD2 = _center.coshD2;
  Complex c_coshL2 = _center.coshL2; 

  char _desc[1000];
  sprintf(_desc, "%s\n", name.c_str());
  sprintf(_desc + strlen(_desc), "sinhP %f + %f I with size %f, absLB %f, and absUB %f\n", sinhP.f.re, sinhP.f.im, sinhP.size, absLB(sinhP), absUB(sinhP));
  sprintf(_desc + strlen(_desc), "sinhD2 %f + %f I with size %f, absLB %f, and absUB %f\n", sinhD2.f.re, sinhD2.f.im, sinhD2.size, absLB(sinhD2), absUB(sinhD2));
  sprintf(_desc + strlen(_desc), "sinhL2 %f + %f I with size %f, absLB %f, and absUB %f\n", sinhL2.f.re, sinhL2.f.im, sinhL2.size, absLB(sinhL2), absUB(sinhL2));
  sprintf(_desc + strlen(_desc), "coshP %f + %f I with size %f, absLB %f, and absUB %f\n", coshP.f.re, coshP.f.im, coshP.size, absLB(coshP), absUB(coshP));
  sprintf(_desc + strlen(_desc), "coshD2 %f + %f I with size %f, absLB %f, and absUB %f\n", coshD2.f.re, coshD2.f.im, coshD2.size, absLB(coshD2), absUB(coshD2));
  sprintf(_desc + strlen(_desc),"coshL2 %f + %f I with size %f, absLB %f, and absUB %f\n", coshL2.f.re, coshL2.f.im, coshL2.size, absLB(coshL2), absUB(coshL2));
  sprintf(_desc + strlen(_desc), "Center sinhP %f + %f I, coshP %f + %f I\n", c_sinhP.real(), c_sinhP.imag(), c_coshP.real(), c_coshP.imag());
  sprintf(_desc + strlen(_desc), "Center sinhD2 %f + %f I, coshD2 %f + %f I\n", c_sinhD2.real(), c_sinhD2.imag(), c_coshD2.real(), c_coshD2.imag());
  sprintf(_desc + strlen(_desc), "Center sinhL2 %f + %f I, coshL2 %f + %f I\n", c_sinhL2.real(), c_sinhL2.imag(), c_coshL2.real(), c_coshL2.imag());

  std::string s(_desc);
  return s;
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
  _center.sinhP = Complex(box_center[3], box_center[0]);
  _center.sinhD2 = Complex(box_center[4], box_center[1]);
  _center.sinhL2 = Complex(box_center[5], box_center[2]);
  Complex one = Complex(1,0);
  _center.coshP = sqrt(_center.sinhP * _center.sinhP + one);
  if (absLB( _center.coshP + _center.sinhP) < 1) {
    _center.coshP = - _center.coshP;
  }   
  _center.coshD2 = sqrt(_center.sinhD2 * _center.sinhD2 + one);
  if (absLB( _center.coshD2 + _center.sinhD2) < 1) {
    _center.coshD2 = - _center.coshD2;
  }   
  _center.coshL2 = sqrt(_center.sinhL2 * _center.sinhL2 + one);
  if (absLB( _center.coshL2 + _center.sinhL2) < 1) {
    _center.coshL2 = - _center.coshL2;
  }
  _x_center = construct_x(_center); 
  _y_center = construct_y(_center); 
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
  _cover.coshP = sqrt( _cover.sinhP  * _cover.sinhP  + 1);
  if (absLB( _cover.coshP + _cover.sinhP) < 1) {
    _cover.coshP = - _cover.coshP;
  }   
  _cover.coshD2 = sqrt( _cover.sinhD2 * _cover.sinhD2 + 1);
  if (absLB( _cover.coshD2 + _cover.sinhD2) < 1) {
    _cover.coshD2 = - _cover.coshD2;
  }   
  _cover.coshL2 = sqrt( _cover.sinhL2 * _cover.sinhL2 + 1);
  if (absLB( _cover.coshL2 + _cover.sinhL2) < 1) {
    _cover.coshL2 = - _cover.coshL2;
  }   
  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
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
	
	_nearer.sinhP  = Complex(m[3], m[0]);
	_nearer.sinhD2 = Complex(m[4], m[1]);
	_nearer.sinhL2 = Complex(m[5], m[2]);
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
	
	_further.sinhP  = Complex(m[3], m[0]);
	_further.sinhD2 = Complex(m[4], m[1]);
	_further.sinhL2 = Complex(m[5], m[2]);
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
	
	_greater.sinhP  = Complex(m[3], m[0]);
	_greater.sinhD2 = Complex(m[4], m[1]);
	_greater.sinhL2 = Complex(m[5], m[2]);
}

//Params<Complex> Box::offset(const double* offset) const
//{
//	Params<Complex> result;
//	result.sinhP = Complex(
//		scale[3]*(offset[3]*size_digits[3] + center_digits[3]),
//		scale[0]*(offset[0]*size_digits[0] + center_digits[0])
//	);
//	result.sinhD2 = Complex(
//		scale[4]*(offset[4]*size_digits[4] + center_digits[4]),
//		scale[1]*(offset[1]*size_digits[1] + center_digits[1])
//	);
//	result.sinhL2 = Complex(
//		scale[5]*(offset[5]*size_digits[5] + center_digits[5]),
//		scale[2]*(offset[2]*size_digits[2] + center_digits[2])
//	);
//	return result;
//}
