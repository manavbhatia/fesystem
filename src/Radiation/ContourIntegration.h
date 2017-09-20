// $Id: ContourIntegration.h,v 1.1.6.1 2008-08-21 00:56:19 manav Exp $

// C/ C++ includes
#include <cmath>

// FESystem includes
#include "FESystem/FESystemNumbers.h"

// libmesh includes
#include "geom/node.h"
#include "geom/elem.h"
#include "fe/fe_base.h"
#include "quadrature/quadrature.h"

inline double	F1(const double &s0,
                 const double &s1,
                 const double &s2,
                 const double &x)	
{
  double b,f;
  b = sqrt(4.0 * s0 * s2 - s1 * s1);
  f = -4.0 * s2 * x;
  f += 2.0 * b * atan2((s1 + 2.0 * s2 * x), b);
  f += (s1 + 2.0 * s2 * x) * log(s0 + x * (s1 + s2 * x));
  f /= (2.0 * s2);
  return f;
}


inline double	F2(const double &s0,
                 const double &s1,
                 const double &s2,
                 const double &x)	
{
  double b,f;
  b = sqrt(4.0 * s0 * s2 - s1 * s1);
  f = 2.0 * (-s1) * b * b * atan2((s1 + 2.0 * s2 * x), b);
  f += b * (2.0 * s2 * x * (s1 - s2 * x) + 
            (2.0 * s2 * (s2 * x * x + s0) - s1 * s1) * log(s0 + x * (s1 + s2 * x)));
  f /= (4.0 * s2 * s2 * b);
  return f;
}


inline double	F3(const double *b,
                 const double &length1,
                 const double &length2)	
{
  double e[5],f,tmp1,tmp2,x;
  bool initialized = false;
  FEType fetype(FIRST, LAGRANGE);
  Node node1(0.,0.,0.,0), node2(0.,0.,0.,1);
  std::auto_ptr<Elem> elem(Elem::build(EDGE2).release());
  std::auto_ptr<FEBase> fe(FEBase::build(1, fetype).release());
  std::auto_ptr<QBase> qbase(QBase::build(QGAUSS, 1, NINTH).release());
  if (!initialized)
    {
    elem->set_node(0) = &node1; elem->set_node(1) = &node2;
    fe->attach_quadrature_rule(qbase.get());
    initialized = true;
    }

  // set the node coordinates and initialize
  node2(0) = length2;
  fe->reinit(elem.get());

  e[0] = 4.0 * b[3] * b[0] - b[1] * b[1];
  e[1] = 4.0 * b[3] * b[2] - 2.0 * b[1] * b[5];
  e[2] = 4.0 * b[3] * b[4] - b[5];
  e[3] = b[1] + 2.0 * b[3] * length1;
  e[4] = b[5];

  // now, iterate over the quad points, and perform the integration
  const std::vector<Point>& xyz = fe->get_xyz();
  const std::vector<double>& JxW = fe->get_JxW();
  
  f = 0.0;
  // this evaluation takes care of both the upper and the lower limits of the 
  // integral. Hence, the atan2 function appears twice
  for (unsigned int i=0; i < xyz.size(); i++)
    {
    x = xyz[i](0);
    tmp1 = sqrtf(e[0] + e[1] * x + e[2] * x * x);
    tmp2 = e[3] + e[4] * x;
    f += JxW[i] * (tmp1 / b[3]) * (atan2(tmp2, tmp1) - atan2((b[1] + b[5] * x), tmp1));
    }
  
  return f;
}



inline double EdgeFactor(const Point& p1,
                         const Point& q1,
                         const Point& p2,
                         const Point& q2)
{
  // calculate the lengths for this pair
  static Point unit_vec1, unit_vec2, vec12;
  double length1, length2, integral, upper_integral, lower_integral, b[6], d[5];
  
  // calculate the lengths
  unit_vec1.assign(q1); unit_vec1 -= p1; length1 = unit_vec1.size();
  unit_vec2.assign(q2); unit_vec2 -= p2; length2 = unit_vec2.size();
  vec12.assign(p2); vec12 -= p1;
  
  // calculate the unit vectors
  unit_vec1 /= length1;
  unit_vec2 /= length2;
  
  // zero the entries
  for (unsigned int i=0; i<6; i++)
    b[i] = 0.0;
  
  for (unsigned int i=0; i<3; i++)
    {
    b[0] += vec12(i) * vec12(i);
    b[1] -= 2.0 * vec12(i) * unit_vec1(i);
    b[2] += 2.0 * vec12(i) * unit_vec2(i);
    b[3] += unit_vec1(i) * unit_vec1(i);
    b[4] += unit_vec2(i) * unit_vec2(i);
    b[5] -= 2.0 * unit_vec1(i) * unit_vec2(i);
    }
  
  // now the d values
  d[0] = b[0] + length1 * (b[1] + length1 * b[3]);
  d[1] = b[2] + length1 * b[5];
  d[2] = b[4];
  d[3] = b[1] + 2.0 * b[3] * length1;
  d[4] = b[5];
  
  
  upper_integral = 
    -2.0 * length1 * length2 + // term1
    (0.5 / b[3]) * (d[3] * F1(d[0], d[1], d[2], length2) +
                    d[4] * F2(d[0], d[1], d[2], length2)) ;
  
  lower_integral = 
    (0.5 / b[3]) * (b[1] * F1(b[0], b[2], b[4], length2) +
                    b[5] * F2(b[0], b[2], b[4], length2));
  
  integral = (0.25 / FESystemNumbers::Pi) * (upper_integral - lower_integral + F3(b, length1, length2));
  return integral;
}


