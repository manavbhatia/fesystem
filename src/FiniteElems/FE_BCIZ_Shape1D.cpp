// $Id:$

/*
 *  FE_BCIZ_Shape1D.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 1/30/09.
 *  Copyright 2009 Virginia Tech. All rights reserved.
 *
 */

// FESystem inlcludes
//#include "FESystem/FESystemExceptions.h"

// libMesh includes
#include "fe.h"
#include "elem.h"




template <>
Real FE<1,BCIZ>::shape(const ElemType t,
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
Real FE<1,BCIZ>::shape(const Elem* elem,
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
Real FE<1,BCIZ>::shape_deriv(const ElemType t,
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
Real FE<1,BCIZ>::shape_deriv(const Elem* elem,
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
Real FE<1,BCIZ>::shape_second_deriv(const ElemType t,
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
Real FE<1,BCIZ>::shape_second_deriv(const Elem* elem,
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
