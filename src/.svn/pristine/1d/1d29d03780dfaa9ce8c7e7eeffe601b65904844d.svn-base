// $Id: GeometricPoint.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_point_h__
#define __geometric_point_h__

// C++ includes


// FESystem includes
#include "GeometricEntity.h"

// libMesh includes

///  this is a struct to store the location of the point
struct Location{double coord1; double coord2; double coord3;};


/// this class defines a geometric point
class GeometricPoint: public GeometricEntity
{
 public:

  /// constructor
  GeometricPoint(GeometricModel& );

  /// destructor
  virtual ~GeometricPoint();

  /// returns the type of this entity
  inline virtual GeometricEntityType type() const;

  /// returns the coordinates of this point
  inline const Location getLocation() const;

  /// method  to check if this object is a duplicate of the given object
  virtual bool checkIfDuplicate(const GeometricEntity& );

 protected:

  /// method to read from input 
  virtual void readFromInput(std::istream& );

  /// method to write to output
  virtual void writeToOutput(std::ostream& );


  /// location of this point
  Location point_coords;
};

inline 
GeometricEntity::GeometricEntityType GeometricPoint::type() const
{
  return GeometricEntity::POINT;
}

inline
const Location GeometricPoint::getLocation() const
{
  return this->point_coords;
}

#endif // __geometric_point_h__
