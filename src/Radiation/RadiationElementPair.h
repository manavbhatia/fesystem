// $Id: RadiationElementPair.h,v 1.11.6.1 2008-08-21 00:56:19 manav Exp $

#ifndef __fesystem_radiation_element_pair_h__
#define __fesystem_radiation_element_pair_h__

// C++ includes


// FESystem includes
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemNumbers.h"
#include "Radiation/RadiationElement.h"

// libmesh includes
#include "geom/node.h"
#include "geom/elem.h"

#ifdef WITH_MATHEMATICA
// mathematica link
#include "mathlink.h"
#endif

// forward declerations
class RadiationElement;
class Elem;
class Node;

class RadiationElementPair
{
 public:
  /// enum describing the method to be used
  enum ShapeFactorMethod
    {
      AUTOMATIC,
      DOUBLE_AREA_INT,
      CONTOUR_INT,
      MITALAS_CONTOUR_INT,
      ANALYTIC
    };
  
  /// constructor
  RadiationElementPair();

  /// destructor
  ~RadiationElementPair();
	
#ifdef WITH_MATHEMATICA
  /// this method sets the mathematica link ID for using Mathematica to perform some computations
  void setMathLink(MLINK mathlink);
#endif
  
  /// initializes the pair with two elements
  void reinit(RadiationElement* el1, 
              RadiationElement* el2);
  
  /// clears the initialization
  void clear();
  
  /// returns the shape factors between the elements as a 
  /// pair, first element is F1->2, and the second is 
  /// F2->1
  std::pair<double, double> getShapeFactor(ShapeFactorMethod method);
	
 protected:
	
  /// calculates and returns the analytic shape factor
  std::pair<double, double> analyticShapeFactor();
  
  /// calculates and returns the double area integration shape factor
  std::pair<double, double> doubleAreaIntShapeFactor();
  	
  /// calculates and returns the contour integration shape factor
  std::pair<double, double> contourIntShapeFactor();

  /// calculates and returns the shape factors using mitalas' method of intgration
  std::pair<double, double> mitalasContourIntShapeFactor();
  
  /// returns the elem in its local coordinated 
  /// @param elem the elem in global space for which local mesh has to be created
  /// @param nodes nodes for the elem will be returned in this vector of nodes
  std::auto_ptr<Elem> createLocalElem(Elem* elem, 
				      std::vector<Node>& nodes);

  /// boolean to keep track of whether the pair has been initialized
  bool initialized;
  
  // Radiation elems
  RadiationElement* rad_elem1;

  // Radiation elems
  RadiationElement* rad_elem2;
  
#ifdef WITH_MATHEMATICA
  /// mathlink for Mathematica
  MLINK mathematica_link;
#endif
};




inline
void
RadiationElementPair::reinit(RadiationElement* el1, 
                             RadiationElement* el2)
{
  // this element should be initialized
  Assert(!this->initialized, ExcInternalError());
  Assert(el1 != NULL, ExcInternalError());
  Assert(el2 != NULL, ExcInternalError());

  this->rad_elem1 = el1;
  this->rad_elem2 = el2;
  
#ifdef DEBUG
  // make sure that the elements "see" each other 
  // for this, the dot product of the surface normal of an element 
  // to any vector to the other surface should be positive
  const Point& normal1 = this->rad_elem1->getSurfaceNormal();
  const Point& normal2 = this->rad_elem2->getSurfaceNormal();
  
  // calculate the vector from node 0 of element 1 to node 0 of elem 2
  static Point vec;
  vec.assign(this->rad_elem2->getRadiationGeometricElem()->point(0));
  vec -= (this->rad_elem1->getRadiationGeometricElem()->point(0));

  const double eps_num = 1.0e-11;
  if ((vec * normal1 + eps_num) < 0)
    {
    std::cout << "(vec * normal1) >= 0.0 violated: " << (vec * normal1) << std::endl;
    std::cout << "normal 1: " << normal1 << std::endl;
    std::cout << "normal 2: " << normal2 << std::endl;
    std::cout << "rad elem 1 FE_ID: " << this->rad_elem1->FeElementID() << "  rad elem 2 FE_ID: " << this->rad_elem2->FeElementID() << std::endl;
    }
  
  if ((vec * normal2 - eps_num) > 0)
    {
    std::cout << "(vec * normal2) <= 0.0 violated: " << (vec * normal2) << std::endl;
    std::cout << "normal 1: " << normal1 << std::endl;
    std::cout << "normal 2: " << normal2 << std::endl;
    std::cout << "rad elem 1 FE_ID: " << this->rad_elem1->FeElementID() << "  rad elem 2 FE_ID: " << this->rad_elem2->FeElementID() << std::endl;
    }

  AssertThrow((vec * normal1 + eps_num) >= 0 , ExcInternalError());
  AssertThrow((vec * normal2 - eps_num) <= 0 , ExcInternalError());
#endif
  
  this->initialized = true;
}




inline
void
RadiationElementPair::clear()
{
  // this element should be initialized
  this->rad_elem1 = NULL;
  this->rad_elem2 = NULL;
  
  this->initialized = false;
}




#endif // __fesystem_radiation_element_pair_h__

