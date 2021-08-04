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

  child.compute_center_and_size();
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
                     XComplex(box_size[1], box_size[3]), 0, 0);

  _cover.sinhD2 = AJCC(XComplex(box_center[0], box_center[2]), 
                     0, XComplex(box_size[0], box_size[2]), 0);

  fill_derived(_cover);

  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
  _cover.coshmu = cosh_move_j(_x_cover); 
}

