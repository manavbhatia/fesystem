// $Id: GeometricLine.h,v 1.3 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_line_h__
#define __geometric_line_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Geometry/GeometricEntity.h"
#include "Geometry/GeometricPoint.h"

// libMesh includes


/// this class defines a geometric point
class GeometricLine: public GeometricEntity
{
 public:

  /// constructor
  GeometricLine(GeometricModel& );

  /// destructor
  virtual ~GeometricLine();

  /// returns the type of this entity
  inline virtual GeometricEntityType type() const = 0;

  /// returns the coordinates of this point
  inline const std::vector<GeometricPoint*>& getPoints() const;

  /// method to check if this iobject is a duplicate of the given object
  virtual bool checkIfDuplicate(const GeometricEntity& ) = 0;

 protected:

  /// method to read from input
  virtual void readFromInput(std::istream& ) = 0;

  /// method to write to output
  virtual void writeToOutput(std::ostream& ) = 0;

  /// location of this point
  std::vector<GeometricPoint*> points;
};


inline
const std::vector<GeometricPoint*>& GeometricLine::getPoints() const
{
  return this->points;
}

#endif // __geometric_line_h__
