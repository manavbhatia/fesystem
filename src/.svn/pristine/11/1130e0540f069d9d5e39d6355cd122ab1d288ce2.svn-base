// $Id: GeometricPoint.C,v 1.4 2006-10-29 05:07:51 manav Exp $

// C++ includes
#include <math.h>
#include <cassert>

// FESystem includes
#include "GeometricPoint.h"

// libMesh includes


GeometricPoint::GeometricPoint(GeometricModel& model):
  GeometricEntity(model)
{

}


GeometricPoint::~GeometricPoint()
{

}


bool GeometricPoint::checkIfDuplicate(const GeometricEntity& point_ref)
{
  // make sure that the two points are not the same object
  assert (&point_ref != this);
  // also, the type should be a point
  if (point_ref.type() != GeometricEntity::POINT)
    return false;

  // convert this point's reference to a geometric point
  const GeometricPoint& point = dynamic_cast<const GeometricPoint&>(point_ref);
  
  // this is the tolerance
  double machine_eps = 1.0e-12;

  // comapre the coordinates, and if they are equal to each other, then the points are duplicate. 
  double x_diff =0.0, y_diff = 0.0, z_diff = 0.0;

  x_diff = fabs(this->point_coords.coord1 - point.point_coords.coord1);
  y_diff = fabs(this->point_coords.coord2 - point.point_coords.coord2);
  z_diff = fabs(this->point_coords.coord3 - point.point_coords.coord3);
 
  if (x_diff <= machine_eps && y_diff <= machine_eps && z_diff <= machine_eps)
    {
      // mark this point as duplicate and set the base entity ID as the ID of the
      // given points base entity ID (since the given entity might be a duplicate
      // of another base entity)
      this->is_duplicate = true;
      this->base_entity_ID = point.baseEntityID();

      return true;
    }
  else 
    return false;
}


void GeometricPoint::readFromInput(std::istream& input)
{
  input >> this->unique_ID;// ID of this point
  input >> this->point_coords.coord1; // coord1 of this point
  input >> this->point_coords.coord2; // coord2
  input >> this->point_coords.coord3; // coord3
}


void GeometricPoint::writeToOutput(std::ostream& output)
{
  output << this->unique_ID << "   "; 
  output << this->point_coords.coord1 << "   ";
  output << this->point_coords.coord2 << "   ";
  output << this->point_coords.coord3 << std::endl;
}
