// $Id: FEBatozShape2D.C,v 1.5.6.1 2007-03-14 22:05:02 manav Exp $

// C++ includes
#include <cassert>

// FESystem inlcludes
//#include "FESystem/FESystemExceptions.h"

// Local includes
#include "fe.h"
#include "elem.h"



namespace FEBatoz
{
  /// returns the length of the side defined by vector from node i to j
  double getEdgeLength(const Elem* elem, const unsigned int i, const unsigned int j)
    {
    const Point& point_i = elem->point(i);
    const Point& point_j = elem->point(j);
    
    static Point vec_ji;
    vec_ji.assign(point_j - point_i);
    
    return vec_ji.size();
    }
  
  /// returns the cos of normal to the side defined by vector from node i to j
  double getEdgeNormalCosine(const Elem* elem, const unsigned int i, const unsigned int j)
    {
    static Point normal, vec1, vec2, vec3;
    
    // calculate the normal to the element
    switch (elem->type())
      {
      case TRI3:
        {
          const Point& point0 = elem->point(0);
          const Point& point1 = elem->point(1);
          const Point& point2 = elem->point(2);
          vec1.assign(point1 - point0);
          vec2.assign(point2 - point0);
          vec3.assign(vec1.cross(vec2));
          // this is the unit vector normal to the plane of the triangle
          vec1.assign(vec3.unit());
          
          // cross product of the length vector with the surface normal will 
          // give the needed vector
          const Point& point_i = elem->point(i);
          const Point& point_j = elem->point(j);
          
          vec2.assign(point_j - point_i);
          vec3.assign(vec2.cross(vec1));
          
          // this is the unit vector needed
          normal.assign(vec3.unit());
          
          // cos of angle between this and the x-axis is simply the 
          // 0th component of this vector
          return normal(0);
        }
        break;
        
      default:
        libmesh_error();
      }
    
    }
  
  /// returns the cos of normal to the side defined by vector from node i to j
  double getEdgeNormalSine(const Elem* elem, const unsigned int i, const unsigned int j)
    {
    // for calculation of sine, it is also important to obtain the sign of the
    // angle. The method to obtain it directly from the cross product will only
    // give the absolute value. There are multiple methods to obtain the sign, however, 
    // the one used here is as follows.
    
    static Point normal, vec1, vec2, vec3;
    
    // calculate the normal to the element
    switch (elem->type())
      {
      case TRI3:
        {
          const Point& point0 = elem->point(0);
          const Point& point1 = elem->point(1);
          const Point& point2 = elem->point(2);
          vec1.assign(point1 - point0);
          vec2.assign(point2 - point0);
          vec3.assign(vec1.cross(vec2));
          // this is the unit vector normal to the plane of the triangle
          vec1.assign(vec3.unit());
          
          // cross product of the length vector with the surface normal will 
          // give the needed vector
          const Point& point_i = elem->point(i);
          const Point& point_j = elem->point(j);
          
          vec2.assign(point_j - point_i);
          vec3.assign(vec2.cross(vec1));
          
          // this is the unit vector needed
          normal.assign(vec3.unit());
          
          // cos of angle between this and the x-axis is simply the 
          // 0th component of this vector
          return normal(1);
        }
        break;
        
      default:
        libmesh_error();
      }
    
    }
}


template <>
Real FE<2,BATOZ>::shape(const ElemType type,
                        const Order order,
                        const unsigned int i,
                        const Point& p)
{
  // params not used for now
  (void) type;
  (void) order;
  (void) i;
  (void) p;

  std::cerr << "Batoz elements require the real element\n"
  << "to construct gradient-based degrees of freedom."
  << std::endl;
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BATOZ>::shape(const Elem* elem,
                        const Order order,
                        const unsigned int i,
                        const Point& p)
{
  assert (elem != NULL);
  
  assert(order == SECOND);
  

  // -- keep in mind that the index numbers for the elems start at 0.
  // -- also, the mid side node numbers in the Batoz's paper are different from 
  //   the ones used in this library. Hence, use the following association
  //                   BATOZ TRI6 node #              libMesh TRI6 node #
  //                          1                               0
  //                          2                               1
  //                          3                               2 
  //                          4 (on edge 2,3)                 4 (on edge 1,2)
  //                          5 (on edge 3,1)                 5 (on edge 2,0)
  //                          6 (on edge 1,2)                 3 (on edge 0,1)
  // -- all shape functions for this element are in a variable major format. Hence, 
  //  they follow Hx for w1, w2, w3, thetax1, thetax2, thetax3, thetay1, thetay2, thetay3.
  // And then the same this is followed for Hy (from dofs 9-17)
  
	switch (elem->type())
	  {
    case TRI3:
      {
        // local static variables for shape functions
        static double N1, N2, N3, N4, N5, N6;
        
        // local static variables for edge lengths and sine/cosines
        static double l12, l23, l31, cos4, cos5, cos6, sin4, sin5, sin6;
        
	N1 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 0, p);
	N2 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 1, p);
	N3 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 2, p);
	N4 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 4, p);
	N5 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 5, p);
	N6 = FE<2,LAGRANGE>::shape(TRI6, SECOND, 3, p);
	
	sin4 = FEBatoz::getEdgeNormalSine(elem, 1, 2);
	sin5 = FEBatoz::getEdgeNormalSine(elem, 2, 0);
	sin6 = FEBatoz::getEdgeNormalSine(elem, 0, 1);
	
	cos4 = FEBatoz::getEdgeNormalCosine(elem, 1, 2);
	cos5 = FEBatoz::getEdgeNormalCosine(elem, 2, 0);
	cos6 = FEBatoz::getEdgeNormalCosine(elem, 0, 1);
	
	l12 = FEBatoz::getEdgeLength(elem, 0, 1);
	l23 = FEBatoz::getEdgeLength(elem, 1, 2);
	l31 = FEBatoz::getEdgeLength(elem, 2, 0);
        

        switch (i)
          {
          case 0:
            {
              // Hx, w1
              return 1.5 * (N5 * sin5 / l31 - N6 * sin6 / l12);
            }
            break;
            
          case 1:
            {
              // Hx, w2
              return 1.5 * (N6 * sin6 / l12 - N4 * sin4 / l23);
            }
            break;
            
          case 2:
            {
              // Hx, w3
              return 1.5 * (N4 * sin4 / l23 - N5 * sin5 / l31);
            }
            break;
            
          case 3:
            {
              // Hx, thetax1
              return (-.75) * (N5 * sin5 * cos5 + N6 * sin6 * cos6);
            }
            break;
            
          case 4:
            {
              // Hx, thetax2
              return (-.75) * (N4 * sin4 * cos4 + N6 * sin6 * cos6);
            }
            break;
            
          case 5:
            {
              // Hx, thetax3
              return (-.75) * (N5 * sin5 * cos5 + N4 * sin4 * cos4);
            }
            break;
            
          case 6:
            {
              // Hx, thetay1
              return (N1 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 7:
            {
              // Hx, thetay2
              return (N2 + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 8:
            {
              // Hx, thetay3
              return (N3 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4));
            }
            break;
            
          case 9:
            {
              // Hy, w1
              return 1.5 * (-N5 * cos5 / l31 + N6 * cos6 / l12);
            }
            break;
            
          case 10:
            {
              // Hy, w2
              return 1.5 * (-N6 * cos6 / l12 + N4 * cos4 / l23);
            }
            break;
            
          case 11:
            {
              // Hy, w3
              return 1.5 * (-N4 * cos4 / l23 + N5 * cos5 / l31);
            }
            break;
            
          case 12:
            {
              // Hy, thetax1
              return (-N1 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 13:
            {
              // Hy, thetax2
              return (-N2 + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 14:
            {
              // Hy, thetax3
              return (-N3 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4));
            }
            break;
            
          case 15:
            {
              return (-1.0) * FE<2,BATOZ>::shape(elem, order, 3, p);
            }
            break;
            
          case 16:
            {
              // Hy, thetay2
              return (-1.0) * FE<2,BATOZ>::shape(elem, order, 4, p);
            }
            break;
            
          case 17:
            {
              // Hy, thetay3
              return (-1.0) * FE<2,BATOZ>::shape(elem, order, 5, p);
            }
            break;
            
          default:
            libmesh_error();
          }
      }
      break;
      
	  case QUAD4:
	    {
        libmesh_error();
	    }
      
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << elem->type()
        << std::endl;
	      libmesh_error();
	    }
	  }
  
  
  // should never get here
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BATOZ>::shape_deriv(const ElemType type,
                              const Order order,
                              const unsigned int i,
                              const unsigned int j,
                              const Point& p)
{
  // params not used here
  (void) type;
  (void) order;
  (void) i;
  (void) j;
  (void) p;

  std::cerr << "Batoz elements require the real element\n"
  << "to construct gradient-based degrees of freedom."
  << std::endl;
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BATOZ>::shape_deriv(const Elem* elem,
                              const Order order,
                              const unsigned int i,
                              const unsigned int j,
                              const Point& p)
{
  assert (elem != NULL);
  
  assert(order == SECOND);
  

  // -- keep in mind that the index numbers for the elems start at 0.
  // -- also, the mid side node numbers in the Batoz's paper are different from 
  //   the ones used in this library. Hence, use the following association
  //                   BATOZ TRI6 node #              libMesh TRI6 node #
  //                          1                               0
  //                          2                               1
  //                          3                               2 
  //                          4 (on edge 2,3)                 4 (on edge 1,2)
  //                          5 (on edge 3,1)                 5 (on edge 2,0)
  //                          6 (on edge 1,2)                 3 (on edge 0,1)
  // -- all shape functions for this element are in a variable major format. Hence, 
  //  they follow Hx for w1, w2, w3, thetax1, thetax2, thetax3, thetay1, thetay2, thetay3.
  // And then the same this is followed for Hy (from dofs 9-17)
  
	switch (elem->type())
	  {
    case TRI3:
      {
        // local static variables for shape functions
        static double N1, N2, N3, N4, N5, N6;
        
        // local static variables for edge lengths and sine/cosines
        static double l12, l23, l31, cos4, cos5, cos6, sin4, sin5, sin6;
        
	N1 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 0, j, p);
	N2 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 1, j, p);
	N3 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 2, j, p);
	N4 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 4, j, p);
	N5 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 5, j, p);
	N6 = FE<2,LAGRANGE>::shape_deriv(TRI6, SECOND, 3, j, p);
	
	sin4 = FEBatoz::getEdgeNormalSine(elem, 1, 2);
	sin5 = FEBatoz::getEdgeNormalSine(elem, 2, 0);
	sin6 = FEBatoz::getEdgeNormalSine(elem, 0, 1);
	
	cos4 = FEBatoz::getEdgeNormalCosine(elem, 1, 2);
	cos5 = FEBatoz::getEdgeNormalCosine(elem, 2, 0);
	cos6 = FEBatoz::getEdgeNormalCosine(elem, 0, 1);
	
	l12 = FEBatoz::getEdgeLength(elem, 0, 1);
	l23 = FEBatoz::getEdgeLength(elem, 1, 2);
	l31 = FEBatoz::getEdgeLength(elem, 2, 0);
	
	
        switch (i)
          {
          case 0:
            {
              // Hx, w1
              return 1.5 * (N5 * sin5 / l31 - N6 * sin6 / l12);
            }
            break;
            
          case 1:
            {
              // Hx, w2
              return 1.5 * (N6 * sin6 / l12 - N4 * sin4 / l23);
            }
            break;
            
          case 2:
            {
              // Hx, w3
              return 1.5 * (N4 * sin4 / l23 - N5 * sin5 / l31);
            }
            break;
            
          case 3:
            {
              // Hx, thetax1
              return (-.75) * (N5 * sin5 * cos5 + N6 * sin6 * cos6);
            }
            break;
            
          case 4:
            {
              // Hx, thetax2
              return (-.75) * (N4 * sin4 * cos4 + N6 * sin6 * cos6);
            }
            break;
            
          case 5:
            {
              // Hx, thetax3
              return (-.75) * (N5 * sin5 * cos5 + N4 * sin4 * cos4);
            }
            break;
            
          case 6:
            {
              // Hx, thetay1
              return (N1 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 7:
            {
              // Hx, thetay2
              return (N2 + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 8:
            {
              // Hx, thetay3
              return (N3 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4));
            }
            break;
            
          case 9:
            {
              // Hy, w1
              return 1.5 * (-N5 * cos5 / l31 + N6 * cos6 / l12);
            }
            break;
            
          case 10:
            {
              // Hy, w2
              return 1.5 * (-N6 * cos6 / l12 + N4 * cos4 / l23);
            }
            break;
            
          case 11:
            {
              // Hy, w3
              return 1.5 * (-N4 * cos4 / l23 + N5 * cos5 / l31);
            }
            break;
            
          case 12:
            {
              // Hy, thetax1
              return (-N1 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 13:
            {
              // Hy, thetax2
              return (-N2 + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 14:
            {
              // Hy, thetax3
              return (-N3 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4));
            }
            break;
            
          case 15:
            {
              return (-1.0) * FE<2,BATOZ>::shape_deriv(elem, order, 3, j, p);
            }
            break;
            
          case 16:
            {
              // Hy, thetay2
              return (-1.0) * FE<2,BATOZ>::shape_deriv(elem, order, 4, j, p);
            }
            break;
            
          case 17:
            {
              // Hy, thetay3
              return (-1.0) * FE<2,BATOZ>::shape_deriv(elem, order, 5, j, p);
            }
            break;
            
          default:
            libmesh_error();
          }
      }
      break;
      
	  case QUAD4:
	    {
        libmesh_error();
	    }
      
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << elem->type()
        << std::endl;
	      libmesh_error();
	    }
	  }
  
  
  // should never get here
  libmesh_error();
  return 0.;
}




template <>
Real FE<2,BATOZ>::shape_second_deriv(const ElemType type,
                                     const Order order,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const Point& p) 
{
  // params not used here
  (void) type;
  (void) order;
  (void) i;
  (void) j;
  (void) p;

  std::cerr << "Batoz elements require the real element\n"
  << "to construct gradient-based degrees of freedom."
  << std::endl;
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,BATOZ>::shape_second_deriv(const Elem* elem,
                                     const Order order,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const Point& p)
{
  assert (elem != NULL);
  
  assert(order == SECOND);
  

  // -- keep in mind that the index numbers for the elems start at 0.
  // -- also, the mid side node numbers in the Batoz's paper are different from 
  //   the ones used in this library. Hence, use the following association
  //                   BATOZ TRI6 node #              libMesh TRI6 node #
  //                          1                               0
  //                          2                               1
  //                          3                               2 
  //                          4 (on edge 2,3)                 4 (on edge 1,2)
  //                          5 (on edge 3,1)                 5 (on edge 2,0)
  //                          6 (on edge 1,2)                 3 (on edge 0,1)
  // -- all shape functions for this element are in a variable major format. Hence, 
  //  they follow Hx for w1, w2, w3, thetax1, thetax2, thetax3, thetay1, thetay2, thetay3.
  // And then the same this is followed for Hy (from dofs 9-17)
  
	switch (elem->type())
	  {
    case TRI3:
      {
        // local static variables for shape functions
        static double N1, N2, N3, N4, N5, N6;
        
        // local static variables for edge lengths and sine/cosines
        static double l12, l23, l31, cos4, cos5, cos6, sin4, sin5, sin6;
        
	N1 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 0, j, p);
	N2 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 1, j, p);
	N3 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 2, j, p);
	N4 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 4, j, p);
	N5 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 5, j, p);
	N6 = FE<2,LAGRANGE>::shape_second_deriv(TRI6, SECOND, 3, j, p);

	sin4 = FEBatoz::getEdgeNormalSine(elem, 1, 2);
	sin5 = FEBatoz::getEdgeNormalSine(elem, 2, 0);
	sin6 = FEBatoz::getEdgeNormalSine(elem, 0, 1);

	cos4 = FEBatoz::getEdgeNormalCosine(elem, 1, 2);
	cos5 = FEBatoz::getEdgeNormalCosine(elem, 2, 0);
	cos6 = FEBatoz::getEdgeNormalCosine(elem, 0, 1);

	l12 = FEBatoz::getEdgeLength(elem, 0, 1);
	l23 = FEBatoz::getEdgeLength(elem, 1, 2);
	l31 = FEBatoz::getEdgeLength(elem, 2, 0);


        switch (i)
          {
          case 0:
            {
              // Hx, w1
              return 1.5 * (N5 * sin5 / l31 - N6 * sin6 / l12);
            }
            break;
            
          case 1:
            {
              // Hx, w2
              return 1.5 * (N6 * sin6 / l12 - N4 * sin4 / l23);
            }
            break;
            
          case 2:
            {
              // Hx, w3
              return 1.5 * (N4 * sin4 / l23 - N5 * sin5 / l31);
            }
            break;
            
          case 3:
            {
              // Hx, thetax1
              return (-.75) * (N5 * sin5 * cos5 + N6 * sin6 * cos6);
            }
            break;
            
          case 4:
            {
              // Hx, thetax2
              return (-.75) * (N4 * sin4 * cos4 + N6 * sin6 * cos6);
            }
            break;
            
          case 5:
            {
              // Hx, thetax3
              return (-.75) * (N5 * sin5 * cos5 + N4 * sin4 * cos4);
            }
            break;
            
          case 6:
            {
              // Hx, thetay1
              return (N1 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 7:
            {
              // Hx, thetay2
              return (N2 + N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4) +
                      N6 * (0.5 * cos6 * cos6 - 0.25 * sin6 * sin6));
            }
            break;
            
          case 8:
            {
              // Hx, thetay3
              return (N3 + N5 * (0.5 * cos5 * cos5 - 0.25 * sin5 * sin5) +
                      N4 * (0.5 * cos4 * cos4 - 0.25 * sin4 * sin4));
            }
            break;
            
          case 9:
            {
              // Hy, w1
              return 1.5 * (-N5 * cos5 / l31 + N6 * cos6 / l12);
            }
            break;
            
          case 10:
            {
              // Hy, w2
              return 1.5 * (-N6 * cos6 / l12 + N4 * cos4 / l23);
            }
            break;
            
          case 11:
            {
              // Hy, w3
              return 1.5 * (-N4 * cos4 / l23 + N5 * cos5 / l31);
            }
            break;
            
          case 12:
            {
              // Hy, thetax1
              return (-N1 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 13:
            {
              // Hy, thetax2
              return (-N2 + N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4) +
                      N6 * (0.25 * cos6 * cos6 - 0.5 * sin6 * sin6));
            }
            break;
            
          case 14:
            {
              // Hy, thetax3
              return (-N3 + N5 * (0.25 * cos5 * cos5 - 0.5 * sin5 * sin5) +
                      N4 * (0.25 * cos4 * cos4 - 0.5 * sin4 * sin4));
            }
            break;
            
          case 15:
            {
              return (-1.0) * FE<2,BATOZ>::shape_second_deriv(elem, order, 3, j, p);
            }
            break;
            
          case 16:
            {
              // Hy, thetay2
              return (-1.0) * FE<2,BATOZ>::shape_second_deriv(elem, order, 4, j, p);
            }
            break;
            
          case 17:
            {
              // Hy, thetay3
              return (-1.0) * FE<2,BATOZ>::shape_second_deriv(elem, order, 5, j, p);
            }
            break;
            
          default:
            libmesh_error();
          }
      }
      break;
      
	  case QUAD4:
	    {
        libmesh_error();
	    }
      
	    
	  default:
	    {
	      std::cerr << "ERROR: Unsupported 2D element type!: " << elem->type()
        << std::endl;
	      libmesh_error();
	    }
	  }
  
  
  // should never get here
  libmesh_error();
  return 0.;
}
