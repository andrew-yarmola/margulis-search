#include "Box.h"
#include "Generators.hh"


double scale[DIM];
static bool scale_initialized = false; 
Box::Box() {
  if (!scale_initialized) {
    scale_initialized = true;
    for (int i = 0; i < DIM; ++i) {
      scale[i] = pow(2, -i / float(DIM));
    }
  }
  for (int i = 0; i < DIM; ++i) {
    center_digits[i] = 0;
    size_digits[i] = SCL;
  }
  pos = 0;
  compute_center_and_size();
  compute_nearer();
  compute_cover();
}

Box Box::child(int dir) const
{
  Box child(*this);
  child.size_digits[pos] *= 0.5;
  child.center_digits[pos] += (2*dir-1)*child.size_digits[pos];
  ++child.pos;
  if (child.pos == DIM) { child.pos = 0; }

  child.name = name;
  child.name.append(1, '0'+dir);

  child.qr = qr;
  child.short_words_cache.clear();

  child.compute_center_and_size();
  child.compute_nearer();
  child.compute_cover();
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

// This is now special to the box mapping
std::string Box::desc() {
  AJCC sinhL2 = _cover.sinhL2;
  AJCC sinhD2 = _cover.sinhD2;
  AJCC coshL2 = _cover.coshL2;
  AJCC coshD2 = _cover.coshD2;
  AJCC coshmu = _cover.coshmu;
  Complex c_sinhL2 = _center.sinhL2;
  Complex c_sinhD2 = _center.sinhD2;

  char _desc[10000];
  sprintf(_desc, "%s\n", name.c_str());
  sprintf(_desc + strlen(_desc), "sinh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhL2.f.re, sinhL2.f.im, sinhL2.size, absLB(sinhL2), absUB(sinhL2));
  sprintf(_desc + strlen(_desc), "cosh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshL2.f.re, coshL2.f.im, coshL2.size, absLB(coshL2), absUB(coshL2));
  sprintf(_desc + strlen(_desc), "sinh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhD2.f.re, sinhD2.f.im, sinhD2.size, absLB(sinhD2), absUB(sinhD2));
  sprintf(_desc + strlen(_desc), "cosh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshD2.f.re, coshD2.f.im, coshD2.size, absLB(coshD2), absUB(coshD2));
  sprintf(_desc + strlen(_desc), "cosh(mu) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshmu.f.re, coshmu.f.im, coshmu.size, absLB(coshmu), absUB(coshmu));

  sprintf(_desc + strlen(_desc),
      "Center\n    sinh(L/2) %f + i %f | sinh(D/2) %f + i %f\n",
      c_sinhL2.real(), c_sinhL2.imag(),  c_sinhD2.real(), c_sinhD2.imag());

  std::string s(_desc);
  return s;
}

void Box::compute_center_and_size()
{
  for (int i = 0; i < DIM; ++i) {
    // GMT paper page 419 of Annals
    // box_size guarantees that :
    // box_center - box_size <= true_center - true_size
    // box_center + box_size >= true_center + true_size
    // where box operations are floating point. 
    box_center[i] = scale[i]*center_digits[i];
    box_size[i]= (1+2*EPS)*(size_digits[i]*scale[i]+HALFEPS*fabs(center_digits[i]));
  }
  _center.sinhL2 = Complex(box_center[1], box_center[3]);
  _center.sinhD2 = Complex(box_center[0], box_center[2]);

  fill_derived(_center);

  _x_center = construct_x(_center); 
  _y_center = construct_y(_center); 
  _center.coshmu = cosh_move_j(_x_center); 
}

void Box::compute_cover()
{
  // Let A = { (z0,z1) \in C^2 | |zi| <= 1 }
  // Our parameters are functions on A with the following defintions
  // sinh(L/2) = (s[1] + i s[3]) z0 + (c[1] + i c[3]) 
  // sinh(D/2) = (s[0] + i s[2]) z1 + (c[0] + i c[2])

  _cover.sinhL2 = AJCC(XComplex(box_center[1], box_center[3]), 
                     XComplex(box_size[1], box_size[3]), 0, 
                     0, 0,
                     0);

  _cover.sinhD2 = AJCC(XComplex(box_center[0], box_center[2]), 
                     0, XComplex(box_size[0], box_size[2]),
                     0, 0,
                     0);

  fill_derived(_cover);

  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
  _cover.coshmu = cosh_move_j(_x_cover); 
}

void Box::compute_nearer()
{
	double m[DIM];
	for (int i = 0; i < DIM; ++i) {
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
	
  _nearer.sinhL2 = Complex(m[1], m[3]);
  _nearer.sinhD2 = Complex(m[0], m[2]);

  fill_derived(_nearer);
}
