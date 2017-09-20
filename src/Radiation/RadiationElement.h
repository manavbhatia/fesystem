// $Id: RadiationElement.h,v 1.9.4.1 2007-03-04 03:37:12 manav Exp $

#ifndef __fesystem_radiation_element_h__
#define __fesystem_radiation_element_h__

// C++ includes
#include <iostream>
#include <vector>
#include <memory>

// FESystem includes
#include "Radiation/RadiationCavityFiniteElement.h"

// libMesh includes
#include "geom/point.h"
//#include "numerics/dense_vector.h"

// Forward declerations
class Point;
class Node;
class Elem;

namespace FESystemNumerics
{
  template<typename T> class PetscSeqVector;
}

class RadiationElement
{
 public:
  
  /// constructor
  /// @param rad_cav_elem the radiation cavity element that contains the details of the finite 
  /// element to which this element belongs
  /// @param fe_elem the geometric elem defining the finite element to which this element belongs
  RadiationElement(RadiationCavityFiniteElem& rad_cav_elem, 
		   const Elem* fe_elem);

  /// destructor
  ~RadiationElement();

  /// init this object
  void init();

  /// ID of the elem
  inline unsigned int ID() const;

  /// ID of the finite element that this is attached to
  inline unsigned int FeElementID() const;

  /// face ID for the finite element
  inline int FeElementFaceID() const;

  
  /// this method sets the data for the element.
  /// @param elem_ID is the local radiation elem ID
  /// @param FE_ID is the finite element ID of the elem
  /// @param side_num is the side number of the element
  /// @param x_div_num is the number of the division along the x-axis that
  /// the element will belong to
  /// @param y_div_num is the number of the division along the y-axis that 
  /// the element will belong to
  inline void setData(unsigned int elem_ID,
		      unsigned int x_div_num,
		      unsigned int y_div_num);

	
  /// this function gets the geom elem for the radiation elem
  Elem* getRadiationGeometricElem();

  /// returns the surface normal for this element
  Point& getSurfaceNormal();
  
  /// returns the area for this element
  double getArea();

  /// writes the information for this element 
  void printInfo(std::ostream& output);

  /// returns the element G vector used for interpolation between FE and 
  /// radiation quantities
  /// @param the vector is returned in this arguement
  void getRadElemG_Vector(FESystemNumerics::PetscSeqVector<double>& vector);

  /// returns a reference to the nodes in global ocordinates for this element
  inline const std::vector<Node*>& getGlobalCoordinateNodes();

  /// returns a reference to the nodes in local coorodinates for this element
  inline const std::vector<Node*>& getLocalCoordinateNodes();


  /// method to init the local radiation element in its own coordinates
  void initLocalRadiationElement();


  /// returns the local radiation element
  inline  Elem* getLocalRadiationGeometricElem();

  
 protected:

  /// radiation cavity finite element that defines the finite element to which this radiation 
  /// element is attached
  RadiationCavityFiniteElem& fe_rad_elem;

  /// geometric finite element to which this element belongs
  const Elem* fe_elem;

  /// element ID
  unsigned int elem_ID;

  /// x div number on FE
  unsigned int x_division_number;

  /// y div number on FE
  unsigned int y_division_number;

  /// geometric FE. This is the element in the Finite Element mesh 
  /// to which this radiation element belongs.
  std::auto_ptr<Elem> geometric_rad_elem;

  /// elem in the local coordinate. This is the geometric element that 
  /// defines the radiation element. This is created from the finite element
  /// to which this element belongs, along with the division number. The 
  /// nodes Different sets of nodes are associated with this element.
  /// Upon creation of the element, the nodes associated with the global coordinates
  /// of this element are stored in global_coordinate_nodes. The non-dimensional 
  /// coordinates of these nodes define their location inside the finite element, 
  /// and these are stores in local_coordinate nodes. For the purpose of
  /// integration inside this element, a local axis is created with its origin 
  /// at node 0 and x-axis along node0->node1 vector. These nodes are stored in
  /// local_nodes vector.
  std::auto_ptr<Elem> local_elem;

  /// area
  double area;
  
  /// says if this is initialized 
  bool initialized;
  
  /// normal 
  Point normal;

  /// vector of nodes for the geometric rad elem. These define the global coordinates 
  /// of this radiation element
  std::vector<Node*> global_coordinate_nodes;

  /// vector of local_coordinate nodes. These are the local non-dimensional 
  /// coordinates of the radiation element inside the FE based upon the 
  /// division number for this element.
  std::vector<Node*> local_coordinate_nodes;

  /// vector of ndoes for the local elem. These define the nodes in the 
  /// coordinate system defined by the origin at node 0 and x-axis 
  /// along node0->node1 vector.
  std::vector<Node*> local_nodes;
};


inline
unsigned int RadiationElement::ID() const
{
  return this->elem_ID;
}


inline
unsigned int RadiationElement::FeElementID() const
{
  return this->fe_rad_elem.elem_ID;
}


inline
int RadiationElement::FeElementFaceID() const
{
  return this->fe_rad_elem.side_num;
}



inline 
void RadiationElement::setData(unsigned int el_ID,
			       unsigned int x_div_num, 
			       unsigned int y_div_num)
{
  this->elem_ID = el_ID;
  this->x_division_number = x_div_num;
  this->y_division_number = y_div_num;
}



inline
double RadiationElement::getArea()
{
  if (!this->initialized)
    this->init();
  
  return this->area;
}


inline 
Point& RadiationElement::getSurfaceNormal()
{
  if (!this->initialized)
    this->init();

  return this->normal;
}


inline 
const std::vector<Node*>& RadiationElement::getGlobalCoordinateNodes()
{
  if (!this->initialized)
    this->init();

  return this->global_coordinate_nodes;
}


inline 
const std::vector<Node*>& RadiationElement::getLocalCoordinateNodes()
{
  if (!this->initialized)
    this->init();

  return this->local_coordinate_nodes;
}


inline 
Elem* RadiationElement::getLocalRadiationGeometricElem()
{
  if (!this->initialized)
    this->init();

  return this->local_elem.get();
}

#endif // __fesystem_radiation_element_h__

