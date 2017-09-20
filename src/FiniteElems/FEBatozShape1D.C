// $Id: FEBatozShape1D.C,v 1.3.6.1 2007-03-14 22:05:02 manav Exp $

// FESystem inlcludes
//#include "FESystem/FESystemExceptions.h"

// libMesh includes
#include "fe.h"
#include "elem.h"




template <>
Real FE<1,BATOZ>::shape(const ElemType t,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  // params not used here
  (void) t;
  (void) order;
  (void) i;
  (void) p;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,BATOZ>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  // params not used here
  (void) order;
  (void) i;
  (void) p;
  (void) elem;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,BATOZ>::shape_deriv(const ElemType t,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  // params not used here
  (void) t;
  (void) order;
  (void) i;
  (void) j;
  (void) p;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,BATOZ>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  // params not used here
  (void) order;
  (void) i;
  (void) j;
  (void) p;
  (void) elem;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}




template <>
Real FE<1,BATOZ>::shape_second_deriv(const ElemType t,
					const Order order,
					const unsigned int i,
					const unsigned int j,
					const Point& p)
{
  // params not used here
  (void) order;
  (void) i;
  (void) j;
  (void) p;
  (void) t;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,BATOZ>::shape_second_deriv(const Elem* elem,
				        const Order order,
				        const unsigned int i,
				        const unsigned int j,
				        const Point& p)
{
  // params not used here
  (void) order;
  (void) i;
  (void) j;
  (void) p;
  (void) elem;

  // there are no 1-D elements for this family.
  libmesh_error();
  return 0.;
}
