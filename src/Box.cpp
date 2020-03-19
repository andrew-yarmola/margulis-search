#include "Box.h"
#include "Generators.hh"

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
//  compute_nearer();
//	compute_further();
//	compute_greater();
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
//  child.compute_nearer();
//	child.compute_further();
//	child.compute_greater();
	return child;
}

Box get_box(std::string code) {
	Box box;
	for (char dir : code) {
		if (dir == '0') {
			box = box.child(0);
		} else if (dir == '1') {
			box = box.child(1);
		}
	}
  return box;
}

std::string Box::desc() {
  AJ sinhdx = _cover.sinhdx;
  AJ coshdx = _cover.coshdx;
  AJ expdx = _cover.expdx;
  AJ expmdx = _cover.expmdx;
  AJ sinhdy = _cover.sinhdy;
  AJ coshdy = _cover.coshdy;
  AJ expdy = _cover.expdy;
  AJ expmdy = _cover.expmdy;
  AJ coshmu = _cover.coshmu; 
  AJ cosf   = _cover.cosf;
  AJ sinf   = _cover.sinf;
  AJ expif  = _cover.expif;
  AJ expmif = _cover.expmif;
  AJ sintx2 = _cover.sintx2;
  AJ sinty2 = _cover.sinty2; 
  AJ coshlx = _cover.coshlx;
  AJ coshly = _cover.coshly;
  AJ coshLx2 = _cover.coshLx2;
  AJ sinhLx2 = _cover.sinhLx2;
  AJ coshLy2 = _cover.coshLy2;
  AJ sinhLy2 = _cover.sinhLy2;
  Complex c_sinhdx = _center.sinhdx;
  Complex c_sinhdy = _center.sinhdy;
  Complex c_coshmu = _center.coshmu; 
  Complex c_cosf   = _center.cosf;
  Complex c_sintx2 = _center.sintx2;
  Complex c_sinty2 = _center.sinty2; 
//  Complex c_coshLx2 = _center.coshLx2;
//  Complex c_sinhLx2 = _center.sinhLx2;
//  Complex c_coshLy2 = _center.coshLy2;
//  Complex c_sinhLy2 = _center.sinhLy2;

  char _desc[10000];
  sprintf(_desc, "%s\n", name.c_str());
  sprintf(_desc + strlen(_desc), "sinh(d_x) = %f with size %f, absLB %f, and absUB %f\n",
                                  sinhdx.f.re, sinhdx.size, absLB(sinhdx), absUB(sinhdx));
  sprintf(_desc + strlen(_desc), "cosh(d_x) = %f with size %f, absLB %f, and absUB %f\n",
                                  coshdx.f.re, coshdx.size, absLB(coshdx), absUB(coshdx));
  sprintf(_desc + strlen(_desc), "exp(d_x) = %f with size %f, absLB %f, and absUB %f\n",
                                  expdx.f.re, expdx.size, absLB(expdx), absUB(expdx));
  sprintf(_desc + strlen(_desc), "exp(-d_x) = %f with size %f, absLB %f, and absUB %f\n",
                                  expmdx.f.re, expmdx.size, absLB(expmdx), absUB(expmdx));
  sprintf(_desc + strlen(_desc), "sinh(d_y) =  %f with size %f, absLB %f, and absUB %f\n",
                                  sinhdy.f.re, sinhdy.size, absLB(sinhdy), absUB(sinhdy));
  sprintf(_desc + strlen(_desc), "cosh(d_y) =  %f with size %f, absLB %f, and absUB %f\n",
                                  coshdy.f.re, coshdy.size, absLB(coshdy), absUB(coshdy));
  sprintf(_desc + strlen(_desc), "exp(d_y) = %f with size %f, absLB %f, and absUB %f\n",
                                  expdy.f.re, expdy.size, absLB(expdy), absUB(expdy));
  sprintf(_desc + strlen(_desc), "exp(-d_y) = %f with size %f, absLB %f, and absUB %f\n",
                                  expmdy.f.re, expmdy.size, absLB(expmdy), absUB(expmdy));
  sprintf(_desc + strlen(_desc), "cosh(mu) = %f with size %f, absLB %f, and absUB %f\n",
                                  coshmu.f.re, coshmu.size, absLB(coshmu), absUB(coshmu));
  sprintf(_desc + strlen(_desc), "cos(phi) = %f with size %f, absLB %f, and absUB %f\n",
                                  cosf.f.re, cosf.size, absLB(cosf), absUB(cosf));
  sprintf(_desc + strlen(_desc), "sin(t_x/2) = %f with size %f, absLB %f, and absUB %f\n",
                                  sintx2.f.re, sintx2.size, absLB(sintx2), absUB(sintx2));
  sprintf(_desc + strlen(_desc), "sin(t_y/2) = %f with size %f, absLB %f, and absUB %f\n",
                                  sinty2.f.re, sinty2.size, absLB(sinty2), absUB(sinty2));
  sprintf(_desc + strlen(_desc), "sin(phi) = %f with size %f, absLB %f, and absUB %f\n",
                                  sinf.f.re, sinf.size, absLB(sinf), absUB(sinf));
  sprintf(_desc + strlen(_desc), "exp(i phi) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  expif.f.re, expif.f.im, expif.size, absLB(expif), absUB(expif));
  sprintf(_desc + strlen(_desc), "exp(-i phi) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  expmif.f.re, expmif.f.im, expmif.size, absLB(expmif), absUB(expmif));
  sprintf(_desc + strlen(_desc), "cosh(lx)) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshlx.f.re, coshlx.f.im, coshlx.size, absLB(coshlx), absUB(coshlx));
  sprintf(_desc + strlen(_desc), "cosh(ly)) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshly.f.re, coshly.f.im, coshly.size, absLB(coshly), absUB(coshly));
  sprintf(_desc + strlen(_desc), "sinh(L_x/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhLx2.f.re, sinhLx2.f.im, sinhLx2.size, absLB(sinhLx2), absUB(sinhLx2));
  sprintf(_desc + strlen(_desc), "cosh(L_x/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshLx2.f.re, coshLx2.f.im, coshLx2.size, absLB(coshLx2), absUB(coshLx2));
  sprintf(_desc + strlen(_desc), "sinh(L_y/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhLy2.f.re, sinhLy2.f.im, sinhLy2.size, absLB(sinhLy2), absUB(sinhLy2));
  sprintf(_desc + strlen(_desc), "cosh(L_y/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshLy2.f.re, coshLy2.f.im, coshLy2.size, absLB(coshLy2), absUB(coshLy2));

  sprintf(_desc + strlen(_desc),
  "Center\n    sinh(d_x) %f | sinh(d_y) %f\n    cosh(mu) %f | cos(phi) % f\n    sin(t_x/2) %f | sin(t_y/2) %f\n",
  c_sinhdx.real(), c_sinhdy.real(), c_coshmu.real(), c_cosf.real(), c_sintx2.real(), c_sinty2.real() );

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
  _center.sinhdx = Complex(16 * box_center[0], 0);
  _center.sinhdy = Complex(16 * box_center[1], 0);
  _center.coshmu = Complex(box_center[2], 0);
  _center.cosf = Complex(box_center[3], 0);
  _center.sintx2 = Complex(box_center[4], 0);
  _center.sinty2 = Complex(box_center[5], 0);

  fill_derived(_center);

  _x_center = construct_x(_center); 
  _y_center = construct_y(_center); 
}

void Box::compute_cover()
{
  // Let A = { (z0,z1,z2) \in C^3 | |zi| <= 1 }
  // Our parameters are functions on A with the following defintions
  // sinh(d_x) = 16(s[0]re(z0) + c[0])
  // sinh(d_y) = 16(s[1]im(z0) + c[1])
  // cosh(mu) = s[2]re[z1] + c[2]
  // cos(phi) = s[3]im[z1] + c[3]
  // sin(t_x/2) = s[4]re[z2] + c[4]
  // sin(t_y/2) = s[5]im[z2] + c[5]

  // We need extra scaling as sinh(d_x) and sinh(d_y) get quit large
  // TODO: verify that this multiplication by powers of 2 is valid
  // Here, we have sinhdx(z0,z1,z2) = 16(size * re(z0) + center)) and sinhdy(z0,z1,z2) = 16(size * im(z0) + center))
  _cover.sinhdx = AJ(XComplex(16 * box_center[0], 0), 
                     XComplex(8 * box_size[0], 0), 0, 0,
                     XComplex(8 * box_size[0], 0), 0, 0);
  // Note: d(im(z0))/dz0 = - i / 2 and d(im(z0)/dconj(z0) = i/2
  _cover.sinhdy = AJ(XComplex(16 * box_center[1], 0), 
                     XComplex(0., - 8 * box_size[1]), 0, 0,
                     XComplex(0.,   8 * box_size[1]), 0, 0);

  // TODO: verify that division by 2 is valid
  _cover.coshmu = AJ(XComplex(box_center[2], 0), 
                  0, XComplex(box_size[2]/2, 0), 0,
                  0, XComplex(box_size[2]/2, 0), 0);

  _cover.cosf = AJ(XComplex(box_center[3], 0), 
                0, XComplex(0, -box_size[3]/2), 0,
                0, XComplex(0,  box_size[3]/2), 0);

  _cover.sintx2 = AJ(XComplex(box_center[4], 0), 
                  0, XComplex(box_size[4]/2, 0), 0,
                  0, XComplex(box_size[4]/2, 0), 0);

  _cover.sinty2 = AJ(XComplex(box_center[5], 0), 
                0, XComplex(0, -box_size[5]/2), 0,
                0, XComplex(0,  box_size[5]/2), 0);

  fill_derived(_cover);

  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
}

/*
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
*/
